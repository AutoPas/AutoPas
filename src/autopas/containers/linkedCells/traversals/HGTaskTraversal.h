/**
 * @file HGTaskTraversal.h
 * @author atacann
 * @date 22.04.2025
 */

#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class ParticleCell_T, class Functor_T>
class HGTaskTraversal : public HGTraversalBase<ParticleCell_T>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell_T::ParticleType;

  explicit HGTaskTraversal(Functor_T *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell_T>(dataLayout, useNewton3), _functor(functor) {}

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Not supported with hgrid_task");
    }
    this->computeIntraLevelInteractions();
    if (this->_numLevels == 1) {
      return;
    }
    // computeInteractions across different levels
    for (size_t upperLevel = 0; upperLevel < this->_numLevels; upperLevel++) {
      if (this->_useNewton3 && upperLevel == 0) {
        // skip the first level as we go top-down only with newton3
        continue;
      }
      // calculate stride for current level
      std::array<size_t, 3> stride = this->computeStride(upperLevel);

      const auto end = this->getTraversalSelectorInfo(upperLevel).cellsPerDim;

      const auto &upperLevelCB = this->_levels->at(upperLevel)->getCellBlock();
      // only look top-down if newton3 is enabled, both ways otherwise
      const size_t levelLimit = this->_useNewton3 ? upperLevel : this->_numLevels;

      const long targetBlocks = static_cast<unsigned long>(autopas_get_max_threads() * 64);
      const std::array<size_t, 3> group = this->findBestGroupSizeForTargetBlocks(targetBlocks, stride, end);

      std::array<size_t, 3> num_index{};
      for (size_t i = 0; i < 3; i++) {
        stride[i] = 1 + static_cast<size_t>(std::ceil(1.0 * (stride[i] - 1) / group[i]));
        num_index[i] = (end[i] + (stride[i] * group[i]) - 1) / (stride[i] * group[i]);
      }
      const size_t stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];
      const size_t numColors = stride_x * stride_y * stride_z;

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
      auto index = [size0, size1, size2](size_t col, size_t xi, size_t yi, size_t zi) -> size_t {
        return ((col * size0 + xi) * size1 + yi) * size2 + zi;
      };

      std::vector<std::array<long, 3>> startIndex(numColors);
      std::vector<std::array<long, 3>> colorDiff(numColors);
      startIndex = this->oneToThreeDForStartSnakyPattern(stride);
      for (int i = 1; i < numColors; i++) {
        colorDiff[i] = startIndex[i] - startIndex[i - 1];
      }
      colorDiff[0] = {0, 0, 0};

      // do the colored traversal
      AUTOPAS_OPENMP(parallel) {
        AUTOPAS_OPENMP(single) {
          for (size_t col = 1; col <= numColors; ++col) {
            for (long zi = 1; zi <= num_index[2]; zi++) {
              for (long yi = 1; yi <= num_index[1]; yi++) {
                for (long xi = 1; xi <= num_index[0]; xi++) {
                  size_t z_start = startIndex[col - 1][2] * group[2] + ((zi - 1) * stride_z * group[2]);
                  size_t y_start = startIndex[col - 1][1] * group[1] + ((yi - 1) * stride_y * group[1]);
                  size_t x_start = startIndex[col - 1][0] * group[0] + ((xi - 1) * stride_x * group[0]);
                  if (!(x_start < end[0] && y_start < end[1] && z_start < end[2])) {
                    continue;
                  }
                  const auto &sameGroup = taskDepend[index(col - 1, xi, yi, zi)];
                  const auto &diffGroup = taskDepend[index(col - 1, xi + colorDiff[col - 1][0],
                                                           yi + colorDiff[col - 1][1], zi + colorDiff[col - 1][2])];
                  const auto &currentTask = taskDepend[index(col, xi, yi, zi)];
                  // Old clang compilers need firstprivate clause
                  AUTOPAS_OPENMP(task firstprivate(z_start, y_start, x_start) depend(in
                                                                                     : sameGroup, diffGroup)
                                     depend(out
                                            : currentTask)) {
                    for (size_t lowerLevel = 0; lowerLevel < levelLimit; lowerLevel++) {
                      if (lowerLevel == upperLevel) {
                        continue;
                      }
                      // get cellBlocks of upper and lower levels
                      const auto &lowerLevelCB = this->_levels->at(lowerLevel)->getCellBlock();

                      // lower bound and upper bound of the owned region of the lower level
                      std::array<size_t, 3> lowerBound = {0, 0, 0},
                                            upperBound = lowerLevelCB.getCellsPerDimensionWithHalo();
                      lowerBound += lowerLevelCB.getCellsPerInteractionLength();
                      upperBound -= lowerLevelCB.getCellsPerInteractionLength();

                      for (size_t z = z_start; z < z_start + group[2]; ++z)
                        for (size_t y = y_start; y < y_start + group[1]; ++y)
                          for (size_t x = x_start; x < x_start + group[0]; ++x) {
                            if (!(x < end[0] && y < end[1] && z < end[2])) {
                              continue;
                            }
                            if (this->_dataLayout == DataLayoutOption::aos) {
                              this->AoSTraversal(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel,
                                                 lowerBound, upperBound, false);
                            } else {
                              this->SoATraversalParticleToCell(lowerLevelCB, upperLevelCB, {x, y, z}, _functor,
                                                               lowerLevel, lowerBound,
                                                               upperBound, false);
                            }
                          }
                    }
                  }
                }
              }
            }
          }
          AUTOPAS_OPENMP(taskwait)
        }
      }
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_task; };

  [[nodiscard]] bool isApplicable() const override { return true; }

 protected:
  Functor_T *_functor;

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) override {
    const auto traversalInfo = this->getTraversalSelectorInfo(level);
    return std::make_unique<LCC08Traversal<ParticleCell_T, Functor_T>>(
        traversalInfo.cellsPerDim, _functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3);
  }
};
}  // namespace autopas
