#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class ParticleCellT, class FunctorT>
class HGTaskSoACellTraversal : public HGTraversalBase<ParticleCellT>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCellT::ParticleType;

  explicit HGTaskSoACellTraversal(FunctorT *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCellT>(dataLayout, useNewton3), _functor(functor) {}

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Not supported with hgrid_task");
    }
    this->computeIntraLevelInteractions();
    // computeInteractions across different levels
    for (size_t upperLevel = 0; upperLevel < this->_numLevels; upperLevel++) {
      // only look top-down if newton3 is enabled, both ways otherwise
      const size_t levelLimit = this->_useNewton3 ? upperLevel : this->_numLevels;
      for (size_t lowerLevel = 0; lowerLevel < levelLimit; lowerLevel++) {
        if (lowerLevel == upperLevel) {
          continue;
        }
        const double interactionLength = this->getInteractionLength(lowerLevel, upperLevel);
        const double interactionLengthSquared = interactionLength * interactionLength;

        std::array<unsigned long, 3> stride = this->computeStride(interactionLength, lowerLevel, upperLevel);

        const std::array<double, 3> dir = {interactionLength, interactionLength, interactionLength};

        const auto end = this->getTraversalSelectorInfo(upperLevel).cellsPerDim;

        // get cellBlocks of upper and lower levels
        const auto &upperLevelCB = this->_levels->at(upperLevel)->getCellBlock();
        auto &lowerLevelLC = *(this->_levels->at(lowerLevel));
        const auto &lowerLevelCB = lowerLevelLC.getCellBlock();

        // lower bound and upper bound of the owned region of the lower level
        std::array<size_t, 3> lowerBound = {0, 0, 0}, upperBound = lowerLevelCB.getCellsPerDimensionWithHalo();
        lowerBound += lowerLevelCB.getCellsPerInteractionLength();
        upperBound -= lowerLevelCB.getCellsPerInteractionLength();

        const long targetBlocksPerColor =
            static_cast<unsigned long>(autopas_get_max_threads() * sqrt(autopas_get_max_threads()));

        unsigned long smallestCriterion = 1e15;
        std::array<unsigned long, 3> bestGroup = {1, 1, 1};
        for (size_t x_group = 1; x_group <= stride[0]; ++x_group)
          for (size_t y_group = 1; y_group <= stride[1]; ++y_group)
            for (size_t z_group = 1; z_group <= stride[2]; ++z_group) {
              std::array<unsigned long, 3> group = {x_group, y_group, z_group};
              std::array<unsigned long, 3> num_index{}, testStride{};
              for (size_t i = 0; i < 3; i++) {
                testStride[i] = 1 + static_cast<unsigned long>(std::ceil(1.0 * (stride[i] - 1) / group[i]));
                num_index[i] = (end[i] + (testStride[i] * group[i]) - 1) / (testStride[i] * group[i]);
              }
              const long numColors = testStride[0] * testStride[1] * testStride[2];
              const long numBlocksPerColor = num_index[0] * num_index[1] * num_index[2];
              const unsigned long criterion = std::abs(numBlocksPerColor * numColors - targetBlocksPerColor);
              if (criterion < smallestCriterion) {
                smallestCriterion = criterion;
                bestGroup = group;
              }
            }
        const std::array<unsigned long, 3> group = bestGroup;
        std::array<unsigned long, 3> num_index{};
        for (size_t i = 0; i < 3; i++) {
          stride[i] = 1 + static_cast<unsigned long>(std::ceil(1.0 * (stride[i] - 1) / group[i]));
          num_index[i] = (end[i] + (stride[i] * group[i]) - 1) / (stride[i] * group[i]);
        }
        const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];
        const unsigned long numColors = stride[0] * stride[1] * stride[2];

        // We need extra space as in the vector example
        const size_t size0 = num_index[0] + 2;
        const size_t size1 = num_index[1] + 2;
        const size_t size2 = num_index[2] + 2;
        const size_t numColorsPlusOne = numColors + 1;

        // Total number of dependency elements.
        const size_t totalDeps = numColorsPlusOne * size0 * size1 * size2;

        // Create a contiguous dependency array.
        std::vector<char> taskDepend(totalDeps, false);

        // Function to compute the 1D index from 4D coordinates.
        auto index = [=](size_t col, size_t xi, size_t yi, size_t zi) -> size_t {
          return ((col * size0 + xi) * size1 + yi) * size2 + zi;
        };

        std::vector<std::array<unsigned long, 3>> startIndex(numColors);
        std::vector<std::array<unsigned long, 3>> colorDiff(numColors);
        startIndex = this->oneToThreeDForStartSnakyPattern(stride);
        for (int i = 1; i < numColors; i++) {
          colorDiff[i] = startIndex[i] - startIndex[i - 1];
        }
        colorDiff[0] = {0, 0, 0};

        // do the colored traversal
        AUTOPAS_OPENMP(parallel) {
          AUTOPAS_OPENMP(single) {
            for (unsigned long col = 1; col <= numColors; ++col) {
              for (unsigned long zi = 1; zi <= num_index[2]; zi++) {
                for (unsigned long yi = 1; yi <= num_index[1]; yi++) {
                  for (unsigned long xi = 1; xi <= num_index[0]; xi++) {
                    unsigned long z_start = startIndex[col - 1][2] * group[2] + ((zi - 1) * stride_z * group[2]);
                    unsigned long y_start = startIndex[col - 1][1] * group[1] + ((yi - 1) * stride_y * group[1]);
                    unsigned long x_start = startIndex[col - 1][0] * group[0] + ((xi - 1) * stride_x * group[0]);
                    if (!(x_start < end[0] && y_start < end[1] && z_start < end[2])) {
                      continue;
                    }
                    AUTOPAS_OPENMP(task depend(
                        in
                        : (taskDepend.data())[index(col - 1, xi, yi, zi)],
                          (taskDepend.data())[index(col - 1, xi + colorDiff[col - 1][0], yi + colorDiff[col - 1][1],
                                                    zi + colorDiff[col - 1][2])])
                                       depend(out
                                              : (taskDepend.data())[index(col, xi, yi, zi)])) {
                      for (unsigned long z = z_start; z < z_start + group[2]; ++z)
                        for (unsigned long y = y_start; y < y_start + group[1]; ++y)
                          for (unsigned long x = x_start; x < x_start + group[0]; ++x) {
                            if (!(x < end[0] && y < end[1] && z < end[2])) {
                              continue;
                            }
                            this->SoATraversalCellToCell(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel,
                                                         interactionLengthSquared, dir, lowerBound, upperBound, true);
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

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_task_soa_cell; };

  [[nodiscard]] bool isApplicable() const override { return this->_dataLayout == DataLayoutOption::soa; }

 protected:
  FunctorT *_functor;

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @param traversalInfo traversal info to generate the new traversal
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) override {
    const auto traversalInfo = this->getTraversalSelectorInfo(level);
    return std::make_unique<LCC08Traversal<ParticleCellT, FunctorT>>(
        traversalInfo.cellsPerDim, _functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3);
  }
};
}  // namespace autopas
