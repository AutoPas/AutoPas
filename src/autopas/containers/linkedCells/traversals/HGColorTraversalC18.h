/**
 * @file HGColorTraversalC18.h
 * @author atacann
 * @date 02.01.2025
 */

#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC18Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class ParticleCell, class Functor>
class HGColorTraversalC18 : public HGTraversalBase<ParticleCell>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell::ParticleType;

  explicit HGColorTraversalC18(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell>(dataLayout, useNewton3), _functor(functor) {}

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @param traversalInfo traversal info to generate the new traversal
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) {
    const auto traversalInfo = this->getTraversalSelectorInfo(level);
    return std::make_unique<LCC18Traversal<ParticleCell, Functor>>(
        traversalInfo.cellsPerDim, _functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3);
  }

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Not supported with hgrid_color");
    }
    //autopas::utils::Timer crosslevel, intralevel, initTraversal;
    // 4D vector (hierarchy level, 3D cell index) to store next non-empty cell in increasing x direction for each cell
    std::vector<std::vector<std::vector<std::vector<size_t>>>> nextNonEmpty(this->_numLevels);
    // computeInteractions for each level independently first
    std::vector<std::unique_ptr<TraversalInterface>> traversals(this->_numLevels);
    for (size_t level = 0; level < this->_numLevels; level++) {
      // TODO: Do this nextNonEmpty calculation after rebuilds in container code, not every iteration, as it only
      // changes after a rebuild
      const auto dim = this->getTraversalSelectorInfo(level).cellsPerDim;

      nextNonEmpty[level] = std::vector<std::vector<std::vector<size_t>>>(
        dim[2], std::vector<std::vector<size_t>>(dim[1], std::vector<size_t>(dim[0])));
      const auto &cellBlock = this->_levels->at(level)->getCellBlock();
      // calculate next non-empty cell in increasing x direction for each cell
      AUTOPAS_OPENMP(parallel for collapse(2))
      for (size_t z = 0; z < dim[2]; ++z) {
        for (size_t y = 0; y < dim[1]; ++y) {
          // calculate 1d index here to not convert from 3d to 1d every iteration of the loop below
          size_t index1D = utils::ThreeDimensionalMapping::threeToOneD({dim[0] - 1, y, z}, dim);
          std::vector<size_t> &nextRef = nextNonEmpty[level][z][y];
          nextRef[dim[0] - 1] = dim[0];
          for (int x = dim[0] - 2; x >= 0; --x, --index1D) {
            if (cellBlock.getCell(index1D).isEmpty()) {
              // if the next cell is empty, set nextRef[x] to nextRef[x + 1]
              nextRef[x] = nextRef[x + 1];
            } else {
              // if next cell is not empty, set nextRef[x] to x + 1 as it is the next non-empty cell
              nextRef[x] = x + 1;
            }
          }
        }
      }

      // generate new traversal, load cells into it, if dataLayout is SoA also load SoA, but do not store SoA
      // as they can still be later used for cross-level interactions
      traversals[level] = generateNewTraversal(level);
      TraversalInterface *temp = traversals[level].get();
      this->_levels->at(level)->prepareTraversal(temp);
      //initTraversal.start();
      traversals[level]->initTraversal();
      //initTraversal.stop();
      //intralevel.start();
      traversals[level]->traverseParticles();
      //intralevel.stop();
      // do not call endTraversal as SoA's should be stored after cross-level interactions are calculated
      // this->_levels->at(level)->computeInteractions(traversals[level].get());
    }
    //crosslevel.start();
    // computeInteractions across different levels
    for (size_t upperLevel = 0; upperLevel < this->_numLevels; upperLevel++) {
      // only look top-down if newton3 is enabled, both ways otherwise
      const size_t levelLimit = this->_useNewton3 ? upperLevel : this->_numLevels;
      for (size_t lowerLevel = 0; lowerLevel < levelLimit; lowerLevel++) {
        if (lowerLevel == upperLevel) {
          continue;
        }
        // calculate cutoff for level pair
        // if there are multiple different sizes in each level, it might make sense to directly get size of the particle
        // on the upper level (but get _cutoff of lowerLevel as it is the upper bound of the cutoff in a level)
        // because currently here interactions between single particle in upper level <-> cells in lower level are
        // calculated
        const double cutoff = (this->_cutoffs[upperLevel] + this->_cutoffs[lowerLevel]) / 2;
        const std::array<double, 3> upperLength = this->_levels->at(upperLevel)->getTraversalSelectorInfo().cellLength;
        const std::array<double, 3> lowerLength = this->_levels->at(lowerLevel)->getTraversalSelectorInfo().cellLength;

        // We only need to check at most distance cutoff + max displacement of any particle in the container
        // NOTE: if in the future, Hgrid will be used as a base container to a verlet list, interactionLength should be
        // always cutoff + _skin.
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
        const double interactionLength = cutoff + this->_maxDisplacement * 2;
#else
        const double interactionLength = cutoff + this->_skin;
#endif
        const double interactionLengthSquared = interactionLength * interactionLength;

        std::array<unsigned long, 3> stride{};
        for (size_t i = 0; i < 3; i++) {
          if (this->_useNewton3) {
            // find out the stride so that cells we check on lowerLevel do not intersect
            if (this->_dataLayout == DataLayoutOption::soa) {
              stride[i] = 1 + static_cast<unsigned long>(std::ceil(std::ceil(interactionLength / lowerLength[i]) * 2 *
                                                                   lowerLength[i] / upperLength[i]));
            }
            // the stride calculation below can result in less colors, but two different threads can operate on SoA
            // Buffer of the same cell at the same time. They won't update the same value at the same time, but
            // causes race conditions in SoA. (Inbetween loading data into vectors (avx etc.) -> functor calcs -> store again).
            else {
              stride[i] = 1 + static_cast<unsigned long>(std::ceil(interactionLength * 2 / upperLength[i]));
            }
          } else {
            // do c01 traversal if newton3 is disabled
            stride[i] = 1;
          }
        }

        const std::array<double, 3> dir = {interactionLength, interactionLength, interactionLength};

        const auto end = this->getTraversalSelectorInfo(upperLevel).cellsPerDim;
        const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

        // get cellBlocks of upper and lower levels
        const auto &upperLevelCB = this->_levels->at(upperLevel)->getCellBlock();
        auto &lowerLevelLC = *(this->_levels->at(lowerLevel));
        const auto &lowerLevelCB = lowerLevelLC.getCellBlock();

        // lower bound and upper bound of the owned region of the lower level
        std::array<size_t, 3> lowerBound = {0, 0, 0}, upperBound = lowerLevelCB.getCellsPerDimensionWithHalo();
        lowerBound += lowerLevelCB.getCellsPerInteractionLength();
        upperBound -= lowerLevelCB.getCellsPerInteractionLength();

        // do the colored traversal
        const unsigned long numColors = stride[0] * stride[1] * stride[2];
        for (unsigned long col = 0; col < numColors; ++col) {
          const std::array<unsigned long, 3> start(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));

          const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
          const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];

          AUTOPAS_OPENMP(parallel for schedule(dynamic, 1) collapse(3))
          for (unsigned long z = start_z; z < end_z; z += stride_z) {
            for (unsigned long y = start_y; y < end_y; y += stride_y) {
              for (unsigned long x = start_x; x < end_x; x += stride_x) {
                const std::array<unsigned long, 3> upperCellIndex3D{x, y, z};
                auto &cell = upperLevelCB.getCell(upperCellIndex3D);
                // skip if cell is a halo cell and newton3 is disabled
                auto isHalo = cell.getPossibleParticleOwnerships() == OwnershipState::halo;
                if (isHalo && !this->_useNewton3) {
                  continue;
                }
                // variable to determine if we are only interested in owned particles in the lower level
                const bool containToOwnedOnly = isHalo && this->_useNewton3;
                if (this->_dataLayout == DataLayoutOption::aos) {
                  for (auto p1Ptr = cell.begin(); p1Ptr != cell.end(); ++p1Ptr) {
                    const std::array<double, 3> &pos = p1Ptr->getR();
                    auto startIndex3D = lowerLevelCB.get3DIndexOfPosition(pos - dir);
                    auto stopIndex3D = lowerLevelCB.get3DIndexOfPosition(pos + dir);
                    // skip halo cells if we need to consider only owned particles
                    if (containToOwnedOnly) {
                      startIndex3D = this->getMax(startIndex3D, lowerBound);
                      stopIndex3D = this->getMin(stopIndex3D, upperBound);
                    }
                    for (size_t zl = startIndex3D[2]; zl <= stopIndex3D[2]; ++zl) {
                      for (size_t yl = startIndex3D[1]; yl <= stopIndex3D[1]; ++yl) {
                        std::vector<size_t> &nextRef = nextNonEmpty[lowerLevel][zl][yl];
                        for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; xl = nextRef[xl]) {
                          const std::array<unsigned long, 3> lowerCellIndex = {xl, yl, zl};
                          auto &lowerCell = lowerLevelCB.getCell(lowerCellIndex);
                          // skip if cell is farther than interactionLength
                          if (this->getMinDistBetweenCellAndPointSquared(lowerLevelCB, lowerCellIndex, pos) >
                              interactionLengthSquared) {
                            continue;
                          }
                          for (auto &p : lowerCell) {
                            this->_functor->AoSFunctor(*p1Ptr, p, this->_useNewton3);
                          }
                        }
                      }
                    }
                  }
                } else {
                  auto &soa = cell._particleSoABuffer;
                  const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
                  const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
                  const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
                  for (int idx = 0; idx < cell.size(); ++idx) {
                    const std::array<double, 3> pos = {xptr[idx], yptr[idx], zptr[idx]};
                    auto soaSingleParticle = soa.constructView(idx, idx + 1);
                    auto startIndex3D = lowerLevelCB.get3DIndexOfPosition(pos - dir);
                    auto stopIndex3D = lowerLevelCB.get3DIndexOfPosition(pos + dir);
                    if (containToOwnedOnly) {
                      startIndex3D = this->getMax(startIndex3D, lowerBound);
                      stopIndex3D = this->getMin(stopIndex3D, upperBound);
                    }
                    for (size_t zl = startIndex3D[2]; zl <= stopIndex3D[2]; ++zl) {
                      for (size_t yl = startIndex3D[1]; yl <= stopIndex3D[1]; ++yl) {
                        std::vector<size_t> &nextRef = nextNonEmpty[lowerLevel][zl][yl];
                        for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; xl = nextRef[xl]) {
                          const std::array<unsigned long, 3> lowerCellIndex = {xl, yl, zl};
                          // skip if cell is farther than interactionLength
                          if (this->getMinDistBetweenCellAndPointSquared(lowerLevelCB, lowerCellIndex, pos) >
                              interactionLengthSquared) {
                            continue;
                          }
                          // 1 to n SoAFunctorPair
                          this->_functor->SoAFunctorPair(soaSingleParticle,
                                                         lowerLevelCB.getCell(lowerCellIndex)._particleSoABuffer,
                                                         this->_useNewton3);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    //crosslevel.stop();
    // AutoPasLog(INFO, "intralevel: {} interlevel: {} initTraversal: {}",
    //   intralevel.getTotalTime() / 1000000000.0, crosslevel.getTotalTime() / 1000000000.0,
    //   initTraversal.getTotalTime() / 1000000000.0);
    // store SoA now if SoA is used
    for (size_t level = 0; level < this->_numLevels; level++) {
      traversals[level]->endTraversal();
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_color_c18; };

  [[nodiscard]] bool isApplicable() const override { return true; }

  void initTraversal() override {}

  void endTraversal() override {}

 protected:
  Functor *_functor;
  double _storedCutoff = 0;
};
}  // namespace autopas
