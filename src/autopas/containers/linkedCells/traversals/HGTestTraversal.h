#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class ParticleCell, class Functor>
class HGTestTraversal : public HGTraversalBase<ParticleCell>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell::ParticleType;

  explicit HGTestTraversal(Functor *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell>(dataLayout, useNewton3), _functor(functor) {}

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @param traversalInfo traversal info to generate the new traversal
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) {
    const auto traversalInfo = this->getTraversalSelectorInfo(level);
    return std::make_unique<LCC08Traversal<ParticleCell, Functor>>(
        traversalInfo.cellsPerDim, _functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3);
  }

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Not supported with hgrid_color");
    }
    //autopas::utils::Timer crosslevel, intralevel, initTraversal;
    // computeInteractions for each level independently first
    std::vector<std::unique_ptr<TraversalInterface>> traversals(this->_numLevels);
    for (size_t level = 0; level < this->_numLevels; level++) {

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
        // calculate openmp for loop bounds
        const unsigned long numColors = stride[0] * stride[1] * stride[2];
        std::array<unsigned long, 3> num_index{};
        for (size_t i = 0; i < 3; i++) {
          num_index[i] = (end[i] + stride[i] - 1) / stride[i];
        }

        // We need extra space as in the vector example
        const size_t size0 = num_index[0] + 2;
        const size_t size1 = num_index[1] + 2;
        const size_t size2 = num_index[2] + 2;
        const size_t numColorsPlusOne = numColors + 1;

        // Total number of dependency elements.
        const size_t totalDeps = numColorsPlusOne * size0 * size1 * size2;

        // Create a contiguous dependency array.
        std::vector<char> taskDepend(totalDeps, false);

        // Inline function (lambda) to compute the 1D index from 4D coordinates.
        auto index = [=](size_t col, size_t xi, size_t yi, size_t zi) -> size_t {
          return ((col * size0 + xi) * size1 + yi) * size2 + zi;
        };

        std::vector<std::array<unsigned long, 3>> startIndex(numColors);
        std::vector<std::array<unsigned long, 3>> colorDiff(numColors);
        startIndex = this->oneToThreeDForStart(stride);
        for (int i = 1; i < numColors; i++) {
          colorDiff[i] = startIndex[i] - startIndex[i - 1];
        }

        AUTOPAS_OPENMP(parallel) {
          AUTOPAS_OPENMP(single) {
            for (unsigned long col = 1; col <= numColors; ++col) {
              for (unsigned long zi = 1; zi <= num_index[2]; zi++) {
                for (unsigned long yi = 1; yi <= num_index[1]; yi++) {
                  for (unsigned long xi = 1; xi <= num_index[0]; xi++) {
                    AUTOPAS_OPENMP(task depend(in:
              taskDepend[index(col-1, xi, yi, zi)],
              taskDepend[index(col-1, xi + colorDiff[col - 1][0], yi + colorDiff[col - 1][1], zi + colorDiff[col - 1][2])]
            ) depend(out: taskDepend[index(col, xi, yi, zi)]) ) {
                      unsigned long z = startIndex[col - 1][2] + (zi - 1) * stride_z;
                      unsigned long y = startIndex[col - 1][1] + (yi - 1) * stride_y;
                      unsigned long x = startIndex[col - 1][0] + (xi - 1) * stride_x;
                      if (x < end[0] && y < end[1] && z < end[2]) {
                        const std::array<unsigned long, 3> upperCellIndex3D{x, y, z};
                        auto &cell = upperLevelCB.getCell(upperCellIndex3D);
                        // skip if cell is a halo cell and newton3 is disabled
                        auto isHalo = cell.getPossibleParticleOwnerships() == OwnershipState::halo;
                        if (!(isHalo && !this->_useNewton3)) {
                          // variable to determine if we are only interested in owned particles in the lower level
                          const bool containToOwnedOnly = isHalo && this->_useNewton3;
                          if (this->_dataLayout == DataLayoutOption::aos) {
                          } else {
                            auto &cell = upperLevelCB.getCell({x, y, z});
                            // skip if cell is a halo cell and newton3 is disabled
                            auto isHalo = cell.getPossibleParticleOwnerships() == OwnershipState::halo;
                            // variable to determine if we are only interested in owned particles in the lower level
                            const bool containToOwnedOnly = isHalo && this->_useNewton3;
                            auto &soa = cell._particleSoABuffer;
                            auto [posMin, posMax] = upperLevelCB.getCellBoundingBox({x, y, z});
                            auto startIndex3D = lowerLevelCB.get3DIndexOfPosition(posMin - dir);
                            auto stopIndex3D = lowerLevelCB.get3DIndexOfPosition(posMax + dir);
                            // skip halo cells if we need to consider only owned particles
                            if (containToOwnedOnly) {
                              startIndex3D = this->getMax(startIndex3D, lowerBound);
                              stopIndex3D = this->getMin(stopIndex3D, upperBound);
                            }

                            for (size_t zl = startIndex3D[2]; zl <= stopIndex3D[2]; ++zl) {
                              for (size_t yl = startIndex3D[1]; yl <= stopIndex3D[1]; ++yl) {
                                for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; ++xl) {
                                  // skip if min distance between the two cells is bigger than interactionLength
                                  // n to n SoAFunctorPair
                                  this->_functor->SoAFunctorPair(cell._particleSoABuffer,
                                                                 lowerLevelCB.getCell({xl, yl, zl})._particleSoABuffer,
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
          }
        }
      }
    }
    // crosslevel.stop();
    // AutoPasLog(INFO, "intralevel: {}s interlevel: {}s initTraversal: {}s", intralevel.getTotalTime() / 1000000000.0,
    //            crosslevel.getTotalTime() / 1000000000.0, initTraversal.getTotalTime() / 1000000000.0);
    // store SoA now if SoA is used
    for (size_t level = 0; level < this->_numLevels; level++) {
      traversals[level]->endTraversal();
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_test; };

  [[nodiscard]] bool isApplicable() const override { return this->_dataLayout == DataLayoutOption::soa; }

  void initTraversal() override {}

  void endTraversal() override {}

 protected:
  Functor *_functor;
  double _storedCutoff = 0;
};
}  // namespace autopas
