#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class ParticleCell_T, class Functor_T>
class HGTestTraversal3 : public HGTraversalBase<ParticleCell_T>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell_T::ParticleType;

  explicit HGTestTraversal3(Functor_T *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCell_T>(dataLayout, useNewton3), _functor(functor) {}

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Not supported with hgrid_color");
    }
    this->computeIntraLevelInteractions();
    if (this->_numLevels == 1) {
      return;
    }
    // computeInteractions across different levels
    for (size_t upperLevel = 0; upperLevel < this->_numLevels; upperLevel++) {
      // calculate stride for current level
      std::array<size_t, 3> stride = this->computeStride(upperLevel);
      const size_t stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

      const auto end = this->getTraversalSelectorInfo(upperLevel).cellsPerDim;
      const size_t end_x = end[0], end_y = end[1], end_z = end[2];

      const auto &upperLevelCB = this->_levels->at(upperLevel)->getCellBlock();
      // only look top-down if newton3 is enabled, both ways otherwise
      const size_t levelLimit = this->_useNewton3 ? upperLevel : this->_numLevels;

      // do the colored traversal
      const size_t numColors = stride_x * stride_y * stride_z;
      AUTOPAS_OPENMP(parallel)
      for (size_t col = 0; col < numColors; ++col) {
        const std::array<size_t, 3> start(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));

        const size_t start_x = start[0], start_y = start[1], start_z = start[2];

        AUTOPAS_OPENMP(for schedule(dynamic, 1) collapse(3))
        for (size_t z = start_z; z < end_z; z += stride_z) {
          for (size_t y = start_y; y < end_y; y += stride_y) {
            for (size_t x = start_x; x < end_x; x += stride_x) {
              for (size_t lowerLevel = 0; lowerLevel < levelLimit; lowerLevel++) {
                if (lowerLevel == upperLevel) {
                  continue;
                }

                // get cellBlocks of upper and lower levels
                const auto &lowerLevelCB = this->_levels->at(lowerLevel)->getCellBlock();

                // lower bound and upper bound of the owned region of the lower level
                std::array<size_t, 3> lowerBound = {0, 0, 0}, upperBound = lowerLevelCB.getCellsPerDimensionWithHalo();
                lowerBound += lowerLevelCB.getCellsPerInteractionLength();
                upperBound -= lowerLevelCB.getCellsPerInteractionLength();
                if (this->_dataLayout == DataLayoutOption::aos) {
                  this->AoSTraversal(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel,
                                     lowerBound, upperBound, false);
                } else {
                  this->SoATraversalParticleToCellUsingCellRange(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel,
                                                   lowerBound, upperBound);
                }
              }
            }
          }
        }
      }
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_test3; };

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