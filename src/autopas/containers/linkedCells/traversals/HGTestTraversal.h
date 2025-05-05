#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class ParticleCell_T, class Functor_T>
class HGTestTraversal : public HGTraversalBase<ParticleCell_T>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell_T::ParticleType;

  explicit HGTestTraversal(Functor_T *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell_T>(dataLayout, useNewton3), _functor(functor) {}

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Not supported with hgrid_block");
    }
    this->computeIntraLevelInteractions();
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

      const long targetBlocksPerColor =
          static_cast<size_t>(autopas_get_max_threads() * sqrt(autopas_get_max_threads()));
      const std::array<size_t, 3> group =
          this->findBestGroupSizeForTargetBlocksPerColor(targetBlocksPerColor, stride, end);

      std::array<size_t, 3> num_index{};
      // calculate openmp for loop bound and actual strides with new block length
      for (size_t i = 0; i < 3; i++) {
        stride[i] = 1 + static_cast<size_t>(std::ceil(1.0 * (stride[i] - 1) / group[i]));
        num_index[i] = (end[i] + (stride[i] * group[i]) - 1) / (stride[i] * group[i]);
      }
      const size_t stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

      // do the colored traversal
      const size_t numColors = stride_x * stride_y * stride_z;
      AUTOPAS_OPENMP(parallel)
      for (size_t col = 0; col < numColors; ++col) {
        const std::array<size_t, 3> startIndex(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));

        AUTOPAS_OPENMP(for schedule(dynamic, 1) collapse(3))
        for (size_t zi = 0; zi < num_index[2]; ++zi) {
          for (size_t yi = 0; yi < num_index[1]; ++yi) {
            for (size_t xi = 0; xi < num_index[0]; ++xi) {
              size_t start_z = startIndex[2] * group[2] + (zi * stride_z * group[2]);
              size_t start_y = startIndex[1] * group[1] + (yi * stride_y * group[1]);
              size_t start_x = startIndex[0] * group[0] + (xi * stride_x * group[0]);
              for (size_t lowerLevel = 0; lowerLevel < levelLimit; lowerLevel++) {
                if (lowerLevel == upperLevel) {
                  continue;
                }
                const double interactionLength = this->getInteractionLength(lowerLevel, upperLevel);
                const double interactionLengthSquared = interactionLength * interactionLength;

                const std::array<double, 3> dir = {interactionLength, interactionLength, interactionLength};

                // get cellBlocks of upper and lower levels
                const auto &lowerLevelCB = this->_levels->at(lowerLevel)->getCellBlock();

                // lower bound and upper bound of the owned region of the lower level
                std::array<size_t, 3> lowerBound = {0, 0, 0}, upperBound = lowerLevelCB.getCellsPerDimensionWithHalo();
                lowerBound += lowerLevelCB.getCellsPerInteractionLength();
                upperBound -= lowerLevelCB.getCellsPerInteractionLength();

                for (size_t z = start_z; z < start_z + group[2]; ++z) {
                  for (size_t y = start_y; y < start_y + group[1]; ++y) {
                    for (size_t x = start_x; x < start_x + group[0]; ++x) {
                      if (!(x < end[0] && y < end[1] && z < end[2])) {
                        continue;
                      }
                      if (this->_dataLayout == DataLayoutOption::aos) {
                        this->AoSTraversal(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel,
                                           interactionLengthSquared, lowerBound, upperBound, false);
                      } else {
                        this->SoATraversalParticleToCell(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel,
                                                         interactionLengthSquared, lowerBound, upperBound, false);
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

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_test; };

  [[nodiscard]] bool isApplicable() const override { return this->_dataLayout == DataLayoutOption::soa; }

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