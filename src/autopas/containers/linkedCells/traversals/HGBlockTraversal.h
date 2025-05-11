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
/**
 * For each level, LCC08Traversal is used. For the cross-level interactions, for each level x only smaller levels
 * are iterated (newton3 on only). The cells on level x are iterated with colors (dynamic color count based on ratio
 * of cell lengths between level x and y) so that the cells on the lower level y
 * that are considered for each cell on level x do not intersect.
 * To reduce number of colors and increase memory efficiency, instead of only 1 upper
 * level cell a block of cells is assigned to a thread at a time. The size of block is calculated dynamically
 * by considering upper and lower cell lengths and number of threads. The number of blocks per color is at least
 * num_threads * 4 or 8, depending on the option.
 * @tparam ParticleCell_T type of Particle cell
 * @tparam Functor_T type of Functor
 */
template <class ParticleCell_T, class Functor_T>
class HGBlockTraversal : public HGTraversalBase<ParticleCell_T>, public HGTraversalInterface {
 public:
  using Particle = typename ParticleCell_T::ParticleType;

  explicit HGBlockTraversal(Functor_T *functor, DataLayoutOption dataLayout, bool useNewton3, int blockMultiplier)
      : HGTraversalBase<ParticleCell_T>(dataLayout, useNewton3), _blockMultiplier(blockMultiplier), _functor(functor) {}

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
    if (not this->isApplicable()) {
      utils::ExceptionHandler::exception("Not supported with hgrid_block");
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

      const int targetBlocksPerColor = autopas_get_max_threads() * _blockMultiplier;
      const std::array<size_t, 3> group =
          this->findBestGroupSizeForTargetBlocksPerColor(targetBlocksPerColor, stride, end);

      std::array<size_t, 3> blocksPerColorPerDim{};
      // calculate openmp for loop bound and actual strides with new block length
      for (size_t i = 0; i < 3; i++) {
        stride[i] = 1 + static_cast<size_t>(std::ceil(1.0 * (stride[i] - 1) / group[i]));
        blocksPerColorPerDim[i] = (end[i] + (stride[i] * group[i]) - 1) / (stride[i] * group[i]);
      }
      const size_t stride_x = stride[0], stride_y = stride[1], stride_z = stride[2];

      // do the colored traversal
      const size_t numColors = stride_x * stride_y * stride_z;
      AutoPasLog(INFO, "HGBlockTraversal: numColors: {}, group: {} {} {}, numBlocksPerColor: {}", numColors, group[0],
                 group[1], group[2], blocksPerColorPerDim[0] * blocksPerColorPerDim[1] * blocksPerColorPerDim[2]);
      AUTOPAS_OPENMP(parallel)
      for (size_t col = 0; col < numColors; ++col) {
        const std::array<size_t, 3> startIndex(utils::ThreeDimensionalMapping::oneToThreeD(col, stride));

        AUTOPAS_OPENMP(for schedule(dynamic, 1) collapse(3))
        for (size_t zi = 0; zi < blocksPerColorPerDim[2]; ++zi) {
          for (size_t yi = 0; yi < blocksPerColorPerDim[1]; ++yi) {
            for (size_t xi = 0; xi < blocksPerColorPerDim[0]; ++xi) {
              size_t start_z = startIndex[2] * group[2] + (zi * stride_z * group[2]);
              size_t start_y = startIndex[1] * group[1] + (yi * stride_y * group[1]);
              size_t start_x = startIndex[0] * group[0] + (xi * stride_x * group[0]);
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

                for (size_t z = start_z; z < start_z + group[2]; ++z) {
                  for (size_t y = start_y; y < start_y + group[1]; ++y) {
                    for (size_t x = start_x; x < start_x + group[0]; ++x) {
                      if (!(x < end[0] && y < end[1] && z < end[2])) {
                        continue;
                      }
                      if (this->_dataLayout == DataLayoutOption::aos) {
                        this->AoSTraversal(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel, upperLevel,
                                           lowerBound, upperBound);
                      } else {
                        this->SoATraversalParticleToCell(lowerLevelCB, upperLevelCB, {x, y, z}, _functor, lowerLevel,
                                                         upperLevel, lowerBound, upperBound);
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

  [[nodiscard]] TraversalOption getTraversalType() const override {
    if (_blockMultiplier == 4) {
      return TraversalOption::hgrid_block4;
    } else if (_blockMultiplier == 8) {
      return TraversalOption::hgrid_block8;
    } else {
      autopas::utils::ExceptionHandler::exception("hgrid_block multiplier is not 4 or 8, blockMultiplier: {}",
                                                  _blockMultiplier);
      return TraversalOption::hgrid_block4;
    }
  };

  [[nodiscard]] bool isApplicable() const override { return true; }

 protected:
  int _blockMultiplier;
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