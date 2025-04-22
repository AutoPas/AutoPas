/**
 * @file HGBlockTraversal.h
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

template <class ParticleCell_T, class Functor_T>
class HGBlockTraversal : public HGTraversalBase<ParticleCell_T>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell_T::ParticleType;

  explicit HGBlockTraversal(Functor_T *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell_T>(dataLayout, useNewton3), _functor(functor) {}

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Not supported with hgrid_block");
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
        const auto &lowerLevelCB = this->_levels->at(lowerLevel)->getCellBlock();

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
              const unsigned long numColors = testStride[0] * testStride[1] * testStride[2];
              const long numBlocksPerColor = num_index[0] * num_index[1] * num_index[2];
              const unsigned long criterion = numColors * std::abs(numBlocksPerColor - targetBlocksPerColor);
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
        AutoPasLog(INFO, "block NumColors {}, numCells {}, num_blocks_per_color {}, groupx {}, groupy {}, groupz {}",
                   numColors, end[0] * end[1] * end[2], num_index[0] * num_index[1] * num_index[2], group[0], group[1],
                   group[2]);

        AUTOPAS_OPENMP(parallel)
        for (unsigned long col = 0; col < numColors; ++col) {
          const std::array<unsigned long, 3> startIndex(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));
          AUTOPAS_OPENMP(for schedule(dynamic, 1) collapse(3))
          for (unsigned long zi = 0; zi < num_index[2]; ++zi) {
            for (unsigned long yi = 0; yi < num_index[1]; ++yi) {
              for (unsigned long xi = 0; xi < num_index[0]; ++xi) {
                unsigned long start_z = startIndex[2] * group[2] + (zi * stride_z * group[2]);
                unsigned long start_y = startIndex[1] * group[1] + (yi * stride_y * group[1]);
                unsigned long start_x = startIndex[0] * group[0] + (xi * stride_x * group[0]);
                for (unsigned long z = start_z; z < start_z + group[2]; ++z) {
                  for (unsigned long y = start_y; y < start_y + group[1]; ++y) {
                    for (unsigned long x = start_x; x < start_x + group[0]; ++x) {
                      if (!(x < end[0] && y < end[1] && z < end[2])) {
                        continue;
                      }
                      if (this->_dataLayout == DataLayoutOption::aos) {
                        this->AoSTraversal(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel,
                                           interactionLengthSquared, dir, lowerBound, upperBound, true);
                      } else {
                        this->SoATraversalParticleToCell(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel,
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

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_block; };

  [[nodiscard]] bool isApplicable() const override { return true; }

 protected:
  Functor_T *_functor;

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @param traversalInfo traversal info to generate the new traversal
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