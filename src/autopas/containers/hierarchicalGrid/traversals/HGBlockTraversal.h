/**
 * @file HGBlockTraversal.h
 * @author atacann
 * @date 02.01.2025
 */

#pragma once

#include "HGTraversalBase.h"
#include "HGTraversalInterface.h"
#include "autopas/containers/linkedCells/traversals/LCC08Traversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {
/**
 * For each level, LCC08Traversal is used. For the cross-level interactions, for each level x, only smaller levels
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
  /**
   * Type of particles stored in the traversed cells.
   */
  using Particle = typename ParticleCell_T::ParticleType;

  /**
   * Constructor.
   * @param functor Pairwise functor used for interactions.
   * @param numLevels Number of levels in the hierarchical grid.
   * @param dataLayout Data layout used by the traversal.
   * @param useNewton3 Whether Newton3 optimization is enabled.
   * @param blockMultiplier Multiplier controlling target blocks per color.
   */
  explicit HGBlockTraversal(Functor_T *functor, size_t numLevels, DataLayoutOption dataLayout, bool useNewton3,
                            int blockMultiplier)
      : HGTraversalBase<ParticleCell_T>(numLevels, dataLayout, useNewton3),
        _blockMultiplier(blockMultiplier),
        _functor(*functor) {}

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
                        this->AoSTraversal(lowerLevelCB, upperLevelCB, {x, y, z}, &_functor, lowerLevel, upperLevel,
                                           lowerBound, upperBound);
                      } else {
                        this->SoATraversalParticleToCell(lowerLevelCB, upperLevelCB, {x, y, z}, &_functor, lowerLevel,
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

 protected:
  /**
   * Multiplier for target blocks per color used for dynamic load balancing.
   */
  int _blockMultiplier;
  /**
   * Functor instance used by this traversal.
   */
  Functor_T &_functor;

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) override {
    const auto traversalInfo = this->getTraversalSelectorInfo(level);
    return std::make_unique<LCC08Traversal<ParticleCell_T, Functor_T>>(
        traversalInfo.cellsPerDim, &_functor, traversalInfo.interactionLength, traversalInfo.cellLength,
        this->_dataLayout, this->_useNewton3);
  }

  [[nodiscard]] bool isApplicable() const override { return this->_numLevels > 1; }

  /**
   * Finds the best group size for a given target number of blocks per color.
   * The block size is selected so that the number of blocks per color is at least the target number, and the number of
   * colors is the smallest possible. If there are two configurations with the same number of colors, the one with the
   * highest number of blocks per color is selected.
   * A block is a unit of work that is assigned to a thread in an OpenMP loop.
   * A group of blocks of 3d size x,y,z will be assigned to a thread per openmp loop iteration.
   * @param targetBlocksPerColor The target number of blocks per color.
   * @param stride The stride for each dimension.
   * @param end The end coordinates for each dimension.
   * @return The best group size for the given target number of blocks per color.
   */
  static std::array<size_t, 3> findBestGroupSizeForTargetBlocksPerColor(int targetBlocksPerColor,
                                                                        const std::array<size_t, 3> &stride,
                                                                        const std::array<unsigned long, 3> &end) {
    unsigned long smallestNumColors = std::numeric_limits<unsigned long>::max();
    unsigned long largestBlocksPerColor = std::numeric_limits<unsigned long>::min();
    std::array<size_t, 3> bestGroup = {1, 1, 1};
    for (size_t x_group = 1; x_group <= std::max(stride[0] - 1, 1ul); ++x_group)
      for (size_t y_group = 1; y_group <= std::max(stride[1] - 1, 1ul); ++y_group)
        for (size_t z_group = 1; z_group <= std::max(stride[2] - 1, 1ul); ++z_group) {
          std::array<size_t, 3> group = {x_group, y_group, z_group};
          std::array<size_t, 3> num_index{}, testStride{};
          for (size_t i = 0; i < 3; i++) {
            testStride[i] = 1 + static_cast<size_t>(
                                    std::ceil((static_cast<double>(stride[i]) - 1.0) / static_cast<double>(group[i])));
            num_index[i] = (end[i] + (testStride[i] * group[i]) - 1) / (testStride[i] * group[i]);
          }
          const size_t numColors = testStride[0] * testStride[1] * testStride[2];
          const unsigned long numBlocksPerColor = num_index[0] * num_index[1] * num_index[2];
          if (numBlocksPerColor >= targetBlocksPerColor and
              (numColors < smallestNumColors or
               (numColors == smallestNumColors and numBlocksPerColor > largestBlocksPerColor))) {
            smallestNumColors = numColors;
            largestBlocksPerColor = numBlocksPerColor;
            bestGroup = group;
          }
        }
    return bestGroup;
  }
};
}  // namespace autopas