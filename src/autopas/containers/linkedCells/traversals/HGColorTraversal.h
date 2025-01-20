/**
 * @file HGColorTraversal.h
 * @author atacann
 * @date 02.01.2025
 */

#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class ParticleCell, class Functor>
class HGColorTraversal : public HGTraversalBase<ParticleCell>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell::ParticleType;

  explicit HGColorTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell>(dataLayout, useNewton3), _functor(functor) {}

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @param traversalInfo traversal info to generate the new traversal
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) {
    const auto traversalInfo = this->getTraversalSelectorInfo(level);
    //_functor->setCutoff(this->_cutoffs[level]);
    return std::make_unique<LCC08Traversal<ParticleCell, Functor>>(
        traversalInfo.cellsPerDim, _functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3);
  }

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Currently only AoS and newton3 is supported on hgrid_color traversal.");
    }
    // computeInteractions for each level independently first
    std::vector<std::unique_ptr<TraversalInterface>> traversals(this->_numLevels);
    for (size_t level = 0; level < this->_numLevels; level++) {
      // generate new traversal, load cells into it, if dataLayout is SoA also load SoA, but do not store SoA
      // as they can still be later used for cross-level interactions
      traversals[level] = generateNewTraversal(level);
      TraversalInterface *temp = traversals[level].get();
      this->_levels->at(level)->prepareTraversal(temp);
      traversals[level]->initTraversal();
      traversals[level]->traverseParticles();
      // do not call endTraversal as it will store SoA
      // this->_levels->at(level)->computeInteractions(traversals[level].get());
    }
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
        //_functor->setCutoff(cutoff);
        // TODO: scale verlet skin by how long passed till last rebuild???
        const double interactionLength = cutoff + this->_skin;

        std::array<unsigned long, 3> stride{};
        for (size_t i = 0; i < 3; i++) {
          if (this->_useNewton3) {
            // find out the stride so that cells we check on lowerLevel do not intersect
            stride[i] = 1 + static_cast<unsigned long>(std::ceil(std::ceil(interactionLength / lowerLength[i]) * 2 *
                                                                 lowerLength[i] / upperLength[i]));
            // stride[i] = 1 + static_cast<unsigned long>(std::ceil(interactionLength  * 2 / upperLength[i]));
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
                auto &cell = upperLevelCB.getCell({x, y, z});
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
                        for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; ++xl) {
                          auto &lowerCell = lowerLevelCB.getCell({xl, yl, zl});
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
                    auto soaView = soa.constructView(idx, idx + 1);
                    auto startIndex3D = lowerLevelCB.get3DIndexOfPosition(pos - dir);
                    auto stopIndex3D = lowerLevelCB.get3DIndexOfPosition(pos + dir);
                    if (containToOwnedOnly) {
                      startIndex3D = this->getMax(startIndex3D, lowerBound);
                      stopIndex3D = this->getMin(stopIndex3D, upperBound);
                    }
                    for (size_t zl = startIndex3D[2]; zl <= stopIndex3D[2]; ++zl) {
                      for (size_t yl = startIndex3D[1]; yl <= stopIndex3D[1]; ++yl) {
                        for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; ++xl) {
                          // TODO: you don't need to check every cell in 3D cubic region, actual region will be like a
                          // sphere
                          // TODO: if (shortest dist between particle and cube) ^ 2 > cutoff ^ 2, can skip this cell
                          // 1 to n SoAFunctorPair
                          this->_functor->SoAFunctorPair(soaView, lowerLevelCB.getCell({xl, yl, zl})._particleSoABuffer,
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
      // store SoA now if SoA is used
      for (size_t level = 0; level < this->_numLevels; level++) {
        traversals[level]->endTraversal();
      }
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_color; };

  [[nodiscard]] bool isApplicable() const override { return true; }

  void initTraversal() override {}

  void endTraversal() override {}

 protected:
  Functor *_functor;
  double _storedCutoff = 0;
};
}  // namespace autopas
