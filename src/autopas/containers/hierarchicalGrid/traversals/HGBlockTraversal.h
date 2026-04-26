/**
 * @file HGBlockTraversal.h
 * @author atacann
 * @date 02.01.2025
 */

#pragma once

#include "HGFitC08Traversal.h"
#include "HGTraversalBase.h"
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
class HGBlockTraversal : public HGTraversalBase<ParticleCell_T>, public TraversalInterface {
 public:
  /**
   * Type of particles stored in the traversed cells.
   */
  using Particle = typename ParticleCell_T::ParticleType;
  /**
   * Type alias for cell blocks used by this traversal.
   */
  using CellBlock = internal::CellBlock3D<FullParticleCell<Particle>>;

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
        TraversalInterface(dataLayout, useNewton3),
        _blockMultiplier(blockMultiplier),
        _functor(*functor) {}

  void traverseParticles() override {
    using namespace autopas::utils::ArrayMath::literals;
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

  [[nodiscard]] bool isApplicable() const override { return true; }

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

  void initTraversal() override {
    this->_traversals.reserve(this->_numLevels);
    for (size_t level = 0; level < this->_numLevels; ++level) {
      this->_traversals.emplace_back(std::make_unique<HGFitC08Traversal<ParticleCell_T, Functor_T>>(
          &_functor, this->_numLevels, this->_dataLayout, this->_useNewton3, level));
      this->_traversals[level]->setLevels(this->_levels, this->_maxCutoffPerLevel, this->_skin);
    }
    // HGFitC08Traversal is an Hgrid traversal which initializes the data layout for all levels
    this->_traversals[0]->initTraversal();
  }

  void endTraversal() override {
    // HGFitC08Traversal is an Hgrid traversal which stores the data layout for all levels
    this->_traversals[0]->endTraversal();
  }

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
   * Traversals used for intra-level interactions for each level.
   */
  std::vector<std::unique_ptr<HGFitC08Traversal>> _traversals;

  /**
   * Compute intra-level interactions for all levels.
   */
  void computeIntraLevelInteractions() {
    for (size_t level = 0; level < this->_numLevels; level++) {
      // We do not simply call computeInteractions() here as we want to store SoA after inter-level traversals are
      // computed. They will be loaded in HGTraversalBase::endTraversal().
      this->_traversals[level]->traverseParticles();
    }
  }

  /**
   * Compute the needed stride to avoid race conditions.
   * For SoA consider that two cells interacting on the same lower level cell is problematic, as two threads can
   * operate on the same buffer.
   * @param level The hierarchy level to compute the stride for.
   * @param topDown If true compute stride for top-down traversal, otherwise bottom-up.
   * @return The stride for each dimension.
   */
  std::array<size_t, 3> computeStride(const size_t level, const bool topDown = true) {
    const std::array<double, 3> levelLength = this->_levels->at(level)->getTraversalSelectorInfo().cellLength;
    std::array<size_t, 3> stride{1, 1, 1};
    for (size_t otherLevel = 0; otherLevel < this->_numLevels; ++otherLevel) {
      if (otherLevel == level) {
        continue;
      }
      if (this->_useNewton3 and ((otherLevel >= level and topDown) or (otherLevel <= level and !topDown))) {
        continue;
      }
      const double interactionLength = this->getInteractionLength(level, otherLevel);
      const std::array<double, 3> otherLevelLength =
          this->_levels->at(otherLevel)->getTraversalSelectorInfo().cellLength;
      std::array<size_t, 3> tempStride{};
      for (size_t i = 0; i < 3; i++) {
        if (this->_useNewton3) {
          // find out the stride so that cells we check on lowerLevel do not intersect
          if (this->_dataLayout == DataLayoutOption::soa) {
            tempStride[i] = 1 + static_cast<size_t>(std::ceil(std::ceil(interactionLength / otherLevelLength[i]) * 2 *
                                                              otherLevelLength[i] / levelLength[i]));
            // the stride calculation below, for AoS, can result in less colors, but two different threads can operate
            // on SoA buffer of the same cell at the same time. They won't update the same value at the same time, but
            // causes race conditions in SoA. (Inbetween loading data into vectors (avx etc.) -> functor calcs ->
            // store again).
          } else {
            tempStride[i] = 1 + static_cast<size_t>(std::ceil(interactionLength * 2 / levelLength[i]));
          }
        } else {
          // do c01 traversal if newton3 is disabled
          tempStride[i] = 1;
        }
        stride[i] = std::max(stride[i], tempStride[i]);
      }
    }
    return stride;
  }

  /**
   * Traverses a single upper-level cell and lower-level cells that are in the interaction range using AoS.
   * @param lowerCB lower level cell block
   * @param upperCB upper level cell block
   * @param upperCellCoords 3d cell index of cell belonging to the upper level
   * @param functor functor to apply
   * @param lowerLevel lower level index
   * @param upperLevel upper level index
   * @param lowerBound inclusive lower coordinate bound for lower-level owned cells when only interacting with owned
   * particles
   * @param upperBound inclusive upper coordinate bound for lower-level owned cells when only interacting with owned
   * particles
   */
  void AoSTraversal(const CellBlock &lowerCB, const CellBlock &upperCB, const std::array<size_t, 3> upperCellCoords,
                    Functor_T *functor, size_t lowerLevel, size_t upperLevel, const std::array<size_t, 3> &lowerBound,
                    const std::array<size_t, 3> &upperBound) {
    // can also use sorted cell optimization? need to be careful to sort only once per upper cell, to not sort for
    // each particle in the upper cell
    using namespace autopas::utils::ArrayMath::literals;

    const bool enableCellDistanceChecking = this->checkDistance(upperLevel, lowerLevel);
    const auto &lowerLevelDims = lowerCB.getCellsPerDimensionWithHalo();
    auto &upperCell = upperCB.getCell(upperCellCoords);
    if (upperCell.isEmpty()) {
      return;
    }
    // skip if cell is a halo cell and newton3 is disabled
    auto isHalo = upperCell.getPossibleParticleOwnerships() == OwnershipState::halo;
    if (isHalo && !this->_useNewton3) {
      return;
    }
    const double lowerInteractionLength = this->_skin + this->_maxCutoffPerLevel[lowerLevel] / 2;
    const std::array<double, 3> lowerDir{lowerInteractionLength, lowerInteractionLength, lowerInteractionLength};
    // variable to determine if we are only interested in owned particles in the lower level
    const bool containToOwnedOnly = isHalo && this->_useNewton3;
    for (auto p1Ptr = upperCell.begin(); p1Ptr != upperCell.end(); ++p1Ptr) {
      const std::array<double, 3> &pos = p1Ptr->getR();
      const double radius = p1Ptr->getSize() / 2;
      const std::array<double, 3> interactionLength = lowerDir + radius;
      const double interactionLengthSquared = interactionLength[0] * interactionLength[0];
      auto startIndex3D = lowerCB.get3DIndexOfPosition(pos - interactionLength);
      auto stopIndex3D = lowerCB.get3DIndexOfPosition(pos + interactionLength);
      // skip halo cells if we need to consider only owned particles
      if (containToOwnedOnly) {
        startIndex3D = utils::ArrayMath::max(startIndex3D, lowerBound);
        stopIndex3D = utils::ArrayMath::min(stopIndex3D, upperBound);
      }
      for (size_t zl = startIndex3D[2]; zl <= stopIndex3D[2]; ++zl) {
        for (size_t yl = startIndex3D[1]; yl <= stopIndex3D[1]; ++yl) {
          auto cellIndex1D = autopas::utils::ThreeDimensionalMapping::threeToOneD(
              {static_cast<size_t>(startIndex3D[0]), yl, zl}, lowerLevelDims);
          for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; ++xl, ++cellIndex1D) {
            auto &lowerCell = lowerCB.getCell(cellIndex1D);
            if (lowerCell.isEmpty()) {
              continue;
            }
            // skip if cell is farther than interactionLength
            if (enableCellDistanceChecking and
                this->getMinDistBetweenCellAndPointSquared(lowerCB, cellIndex1D, pos) > interactionLengthSquared) {
              continue;
            }
            for (auto &p : lowerCell) {
              functor->AoSFunctor(*p1Ptr, p, this->_useNewton3);
            }
          }
        }
      }
    }
  }

  /**
   * Traverses a single upper level cell and lower level cells that are in the interaction range using SoA.
   * SoA is calculated between a single upper level particle and a full lower level cell.
   * @param lowerCB lower level cell block
   * @param upperCB upper level cell block
   * @param upperCellCoords 3d cell index of cell belonging to the upper level
   * @param functor functor to apply
   * @param lowerLevel lower level index
   * @param upperLevel upper level index
   * @param lowerBound lower bound of the coordinates of lower level cells that contain owned particles
   * @param upperBound upper bound of the coordinates of lower level cells that contain owned particles
   */
  void SoATraversalParticleToCell(const CellBlock &lowerCB, const CellBlock &upperCB,
                                  const std::array<size_t, 3> upperCellCoords, Functor_T *functor, size_t lowerLevel,
                                  size_t upperLevel, const std::array<size_t, 3> &lowerBound,
                                  const std::array<size_t, 3> &upperBound) {
    using namespace autopas::utils::ArrayMath::literals;
    const bool enableCellDistanceChecking = this->checkDistance(upperLevel, lowerLevel);
    const auto &lowerLevelDims = lowerCB.getCellsPerDimensionWithHalo();
    auto &upperCell = upperCB.getCell(upperCellCoords);
    if (upperCell.isEmpty()) {
      return;
    }
    // skip if cell is a halo cell and newton3 is disabled
    auto isHalo = upperCell.getPossibleParticleOwnerships() == OwnershipState::halo;
    if (isHalo && !this->_useNewton3) {
      return;
    }
    // variable to determine if we are only interested in owned particles in the lower level
    const bool containToOwnedOnly = isHalo && this->_useNewton3;
    auto &soa = upperCell._particleSoABuffer;
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const double lowerInteractionLength = this->_skin + this->_maxCutoffPerLevel[lowerLevel] / 2;
    const std::array<double, 3> lowerDir{lowerInteractionLength, lowerInteractionLength, lowerInteractionLength};

    for (int idx = 0; idx < upperCell.size(); ++idx) {
      const std::array<double, 3> pos = {xptr[idx], yptr[idx], zptr[idx]};
      const double radius = upperCell[idx].getSize() / 2;
      const std::array<double, 3> interactionLength = lowerDir + radius;
      const double interactionLengthSquared = interactionLength[0] * interactionLength[0];
      auto soaSingleParticle = soa.constructView(idx, idx + 1);
      auto startIndex3D = lowerCB.get3DIndexOfPosition(pos - interactionLength);
      auto stopIndex3D = lowerCB.get3DIndexOfPosition(pos + interactionLength);
      if (containToOwnedOnly) {
        startIndex3D = utils::ArrayMath::max(startIndex3D, lowerBound);
        stopIndex3D = utils::ArrayMath::min(stopIndex3D, upperBound);
      }

      for (size_t zl = startIndex3D[2]; zl <= stopIndex3D[2]; ++zl) {
        for (size_t yl = startIndex3D[1]; yl <= stopIndex3D[1]; ++yl) {
          auto cellIndex1D = autopas::utils::ThreeDimensionalMapping::threeToOneD(
              {static_cast<size_t>(startIndex3D[0]), yl, zl}, lowerLevelDims);
          for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; ++xl, ++cellIndex1D) {
            auto &lowerCell = lowerCB.getCell(cellIndex1D);
            if (lowerCell.isEmpty()) {
              continue;
            }
            // skip if cell is farther than interactionLength
            if (enableCellDistanceChecking and
                this->getMinDistBetweenCellAndPointSquared(lowerCB, cellIndex1D, pos) > interactionLengthSquared) {
              continue;
            }
            // 1 to n SoAFunctorPair
            functor->SoAFunctorPair(soaSingleParticle, lowerCell._particleSoABuffer, this->_useNewton3);
          }
        }
      }
    }
  }

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