/**
 * @file HGColorSoACellToCell.h
 * @author atacann
 * @date 17.01.2025
 */

#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {

template <class ParticleCellT, class FunctorT>
class HGColorSoACellToCell : public HGTraversalBase<ParticleCellT>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCellT::ParticleType;

  explicit HGColorSoACellToCell(FunctorT *functor, DataLayoutOption dataLayout, bool useNewton3)
      : HGTraversalBase<ParticleCellT>(dataLayout, useNewton3), _functor(functor) {}

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Not supported with hgrid_color_soa_cell");
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
        const unsigned long stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

        // get cellBlocks of upper and lower levels
        const auto &upperLevelCB = this->_levels->at(upperLevel)->getCellBlock();
        const auto &lowerLevelCB = this->_levels->at(lowerLevel)->getCellBlock();

        // lower bound and upper bound of the owned region of the lower level
        std::array<size_t, 3> lowerBound = {0, 0, 0}, upperBound = lowerLevelCB.getCellsPerDimensionWithHalo();
        lowerBound += lowerLevelCB.getCellsPerInteractionLength();
        upperBound -= lowerLevelCB.getCellsPerInteractionLength();

        // do the colored traversal
        const unsigned long numColors = stride[0] * stride[1] * stride[2];
        AUTOPAS_OPENMP(parallel)
        for (unsigned long col = 0; col < numColors; ++col) {
          const std::array<unsigned long, 3> start(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));

          const unsigned long start_x = start[0], start_y = start[1], start_z = start[2];
          const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];

          AUTOPAS_OPENMP(for schedule(dynamic, 1) collapse(3))
          for (unsigned long z = start_z; z < end_z; z += stride_z) {
            for (unsigned long y = start_y; y < end_y; y += stride_y) {
              for (unsigned long x = start_x; x < end_x; x += stride_x) {
                this->SoATraversalCellToCell(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel,
                                             interactionLengthSquared, dir, lowerBound, upperBound, true);
              }
            }
          }
        }
      }
    }
  }

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::hgrid_color_soa_cell; };

  [[nodiscard]] bool isApplicable() const override { return this->_dataLayout == DataLayoutOption::soa; }

 protected:
  FunctorT *_functor;

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
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
