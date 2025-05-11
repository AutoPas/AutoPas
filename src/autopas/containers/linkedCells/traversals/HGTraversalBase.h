/**
 * @file HGTraversalBase.h
 * @author atacann
 * @date 09.12.2024
 */

#pragma once

#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {
/**
 * Class for common operations of HGridTraversal classes. It should be a parent class for all HGrid traversals.
 * @tparam ParticleCell_T type of Particle cell
 */
template <class ParticleCell_T>
class HGTraversalBase : public TraversalInterface {
 public:
  using Particle = typename ParticleCell_T::ParticleType;

  explicit HGTraversalBase(DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3),
        _numLevels(0),
        _levels(nullptr),
        _skin(0),
        _maxDisplacement(0),
        _stepsSinceLastRebuild(0),
        _rebuildFrequency(1) {}

  /**
   * Store HGrid data
   * @param levels LinkedCells for each HGrid level
   * @param cutoffs cutoffs of each HGrid level
   * @param skin verlet skin of HGrid container
   * @param maxDisplacement maximum displacement of any particle in the container, ignored if dynamic containers is not
   * used
   * @param stepsSinceLastRebuild number of time-steps since last rebuild
   * @param rebuildFrequency frequency of rebuild
   */
  void setLevels(std::vector<std::unique_ptr<LinkedCells<Particle>>> *levels, std::vector<double> &cutoffs, double skin,
                 double maxDisplacement, unsigned int stepsSinceLastRebuild, const unsigned int rebuildFrequency) {
    _numLevels = cutoffs.size();
    _cutoffs = cutoffs;
    _levels = levels;
    _skin = skin;
    _maxDisplacement = maxDisplacement;
    _stepsSinceLastRebuild = stepsSinceLastRebuild;
    _rebuildFrequency = rebuildFrequency;
  }

 protected:
  size_t _numLevels;
  std::vector<std::unique_ptr<LinkedCells<Particle>>> *_levels;
  std::vector<double> _cutoffs;
  double _skin;
  double _maxDisplacement;
  unsigned int _stepsSinceLastRebuild;
  unsigned int _rebuildFrequency;
  const double distCheckRatio = 0.15;
  // Intralevel traversals used for each level for this iteration
  std::vector<std::unique_ptr<TraversalInterface>> _traversals;

  using CellBlock = internal::CellBlock3D<FullParticleCell<Particle>>;

  /**
   * Generate a new Traversal from the given data, needed as each level of HGrid has different cell sizes
   * @param level which HGrid level to generate a traversal for
   * @return A new traversal that is applicable to a specific LinkedCells level
   */
  virtual std::unique_ptr<TraversalInterface> generateNewTraversal(const size_t level) = 0;

  /**
   * Get the interaction length for a given level pair.
   * This is simply the average of the cutoffs of the two levels plus the verlet skin.
   * @param lowerLevel The lower level index.
   * @param upperLevel The upper level index.
   * @return The interaction length for the given level pair.
   */
  [[nodiscard]] double getInteractionLength(const size_t lowerLevel, const size_t upperLevel) const {
    const double cutoff = (this->_cutoffs[upperLevel] + this->_cutoffs[lowerLevel]) / 2;
    return cutoff + currentSkin();
  }

  /**
   * Check if ratio between two levels is small enough to consider checking distance of each cell to early break
   * and not call functors.
   * @param upperLevel The upper level index.
   * @param lowerLevel The lower level index.
   * @return True if the distance is small enough, false otherwise.
   */
  [[nodiscard]] bool checkDistance(const size_t upperLevel, const size_t lowerLevel) const {
    // size of the lower level cell divided by interactionLength
    // We don't use cellLength because cellSizeFactor might have influenced it
    return (_cutoffs[lowerLevel] + _skin) / getInteractionLength(lowerLevel, upperLevel) <= distCheckRatio;
  }

  [[nodiscard]] double currentSkin() const {
    // We only need to check at most distance cutoff + max displacement of any particle in the container if dynamic
    // containers are used.
    // NOTE: if in the future, if Hgrid will be used as a base container to a verlet list, interactionLength should be
    // always cutoff + _skin.
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    return std::min(this->_maxDisplacement * 2 + 1e-9, this->_skin);
#else
    return std::min((this->_skin / _rebuildFrequency) * _stepsSinceLastRebuild + 1e-9, this->_skin);
#endif
  }

  /**
   * Compute the color for the given level.
   * Computed the needed stride by checking all lower levels and takes the maximum for each dimension.
   * @param level The hierarchy level to compute the stride for.
   * @return The stride for each dimension.
   */
  std::array<size_t, 3> computeStride(const size_t level) {
    const std::array<double, 3> levelLength = this->_levels->at(level)->getTraversalSelectorInfo().cellLength;
    std::array<size_t, 3> stride{1, 1, 1};
    for (size_t otherLevel = 0; otherLevel < this->_numLevels; ++otherLevel) {
      if (otherLevel == level) {
        continue;
      }
      if (this->_useNewton3 && otherLevel >= level) {
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
            // the stride calculation below can result in less colors, but two different threads can operate on SoA
            // buffer of the same cell at the same time. They won't update the same value at the same time, but
            // causes race conditions in SoA. (Inbetween loading data into vectors (avx etc.) -> functor calcs -> store
            // again).
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
   *
   * @tparam T Coordinate integral type
   * @param x x coordinate
   * @param y y coordinate
   * @param z z coordinate
   * @param dim0 x dimension size
   * @param dim1 y dimension size
   * @return 1D index of the 3D coordinates. The 3d manhattan distance between consecutive 1d indexes is 1.
   */
  template <typename T>
  static T computeIndex(T x, T y, T z, T dim0, T dim1) {
    int term = y * dim0 + ((y % 2 == 0) ? x : (dim0 - 1 - x));
    if (z % 2 == 0) {
      return z * dim0 * dim1 + term;
    } else {
      return (z + 1) * dim0 * dim1 - 1 - term;
    }
  }

  /**
   * Convert a 1d index (or color) to a 3d index in a snaky pattern so that between consecutive indexes,
   * the manhattan distance of the 3D coordinates is 1. It is used in TaskTraversal.
   * @tparam T Type of the indices.
   * @param dims The total dimensions of the index space.
   * @return A vector that contains the 3d indexes for each 1d index.
   */
  template <typename T>
  std::vector<std::array<long, 3>> oneToThreeDForStartSnakyPattern(const std::array<T, 3> &dims) {
    std::vector<std::array<long, 3>> coords(dims[0] * dims[1] * dims[2]);
    for (T z = 0; z < dims[2]; ++z) {
      for (T y = 0; y < dims[1]; ++y) {
        for (T x = 0; x < dims[0]; ++x) {
          int idx = computeIndex(x, y, z, dims[0], dims[1]);
          coords[idx] = {static_cast<long>(x), static_cast<long>(y), static_cast<long>(z)};
        }
      }
    }
    return coords;
  }

  /**
   *
   * @param upperCB first CellBlock3D
   * @param upperCell 3d cell index of cell belonging to first CellBlock3D
   * @param lowerCB second CellBlock3D
   * @param lowerCell 3d cell index of cell belonging to second CellBlock3D
   * @return minimum distance between two cells squared
   */
  double getMinDistBetweenCellsSquared(const CellBlock &upperCB, const std::array<size_t, 3> &upperCell,
                                       const CellBlock &lowerCB, const std::array<size_t, 3> &lowerCell) {
    const auto cellPos1 = upperCB.getCellBoundingBox(upperCell);
    const auto cellPos2 = lowerCB.getCellBoundingBox(lowerCell);
    double totalDist = 0;
    for (size_t i = 0; i < 3; ++i) {
      const auto dist =
          std::max(0.0, std::max(cellPos1.first[i] - cellPos2.second[i], cellPos2.first[i] - cellPos1.second[i]));
      totalDist += dist * dist;
    }
    return totalDist;
  }

  /**
   *
   * @param CB cell block that contains the cell
   * @param cellIndex 1d cell index of cell belonging to first CellBlock3D
   * @param point
   * @return minimum distance between a cell and a point squared
   */
  double getMinDistBetweenCellAndPointSquared(const CellBlock &CB, const size_t cellIndex,
                                              const std::array<double, 3> &point) {
    const auto cellPos = CB.getCellBoundingBox(cellIndex);
    double totalDist = 0;
    for (size_t i = 0; i < 3; ++i) {
      const auto dist = std::max(0.0, std::max(cellPos.first[i] - point[i], point[i] - cellPos.second[i]));
      totalDist += dist * dist;
    }
    return totalDist;
  }

  /**
   *
   * @tparam T type of elements in arrays
   * @tparam SIZE size of arrays
   * @param a first array
   * @param b second array
   * @return max of two arrays element-wise
   */
  template <class T, std::size_t SIZE>
  constexpr std::array<T, SIZE> getMax(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
    std::array<T, SIZE> ret{};
    for (size_t i = 0; i < SIZE; ++i) {
      ret[i] = std::max(a[i], b[i]);
    }
    return ret;
  }

  /**
   *
   * @tparam T type of elements in arrays
   * @tparam SIZE size of arrays
   * @param a first array
   * @param b second array
   * @return min of two arrays element-wise
   */
  template <class T, std::size_t SIZE>
  constexpr std::array<T, SIZE> getMin(const std::array<T, SIZE> &a, const std::array<T, SIZE> &b) {
    std::array<T, SIZE> ret{};
    for (size_t i = 0; i < SIZE; ++i) {
      ret[i] = std::min(a[i], b[i]);
    }
    return ret;
  }

  /**
   * Returns traversalSelectorInfo for the specific level
   * @param level Which level to get info from
   * @return TraversalSelectorInfo of level
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo(int level) const {
    if (level < 0 || level >= _numLevels) {
      autopas::utils::ExceptionHandler::exception("Hierarchical grid level out of range: {}", level);
    }
    // return traversal info of hierarchy with biggest interactionLength
    const auto temp = _levels->at(level)->getTraversalSelectorInfo();
    // adjust interactionLength to actual value
    // _overlap will be smaller than needed, because smaller levels have big halo regions compared to their actual
    // interactionlength so cells past _overlap can also be halo cells, resulting in unnecessary iteration of those halo
    // cells
    const TraversalSelectorInfo ret{temp.cellsPerDim, this->getInteractionLength(level, level), temp.cellLength,
                                    temp.clusterSize};
    return ret;
  }

  /**
   * Compute intra-level interactions for all levels.
   */
  void computeIntraLevelInteractions() {
    for (size_t level = 0; level < this->_numLevels; level++) {
      // We do not simply call computeInteractions() here as we want to store SoA after inter-level traversals are
      // computed. They will be loaded in HGTraversalBase::endTraversal().
      TraversalInterface *temp = this->_traversals[level].get();
      // set actual halo region length of the traversal
      // this will let traversal skip going through unnecessary halo cells
      const double haloRegionLength = _cutoffs.back() + _skin;
      temp->setHaloRegionLength(haloRegionLength);
      this->_levels->at(level)->prepareTraversal(temp);
      this->_traversals[level]->initTraversal();
      this->_traversals[level]->traverseParticles();
    }
  }

  void initTraversal() override {
    this->_traversals.resize(this->_numLevels);
    for (size_t i = 0; i < this->_numLevels; ++i) {
      this->_traversals[i] = this->generateNewTraversal(i);
    }
  }

  void endTraversal() override {
    // stores SoA if SoA is used
    for (size_t level = 0; level < this->_numLevels; level++) {
      this->_traversals[level]->endTraversal();
    }
  }

  /**
   * Traverses a single upper level cell and a lower level cells that are in the interaction range using AoS.
   * @tparam Functor type of the functor
   * @param lowerCB lower level cell block
   * @param upperCB upper level cell block
   * @param upperCellCoords 3d cell index of cell belonging to the upper level
   * @param functor functor to apply
   * @param lowerLevel lower level index
   * @param lowerBound lower bound of the coordinates of lower level cells that contain owned particles
   * @param upperBound upper bound of the coordinates of lower level cells that contain owned particles
   * @param distanceCheck if true, checks if the minimum distance between upper level particle and lower level cell is
   * less than the interaction length, otherwise skips the lower level cell
   */
  template <class Functor>
  void AoSTraversal(const CellBlock &lowerCB, const CellBlock &upperCB, const std::array<size_t, 3> upperCellCoords,
                    Functor *functor, size_t lowerLevel, size_t upperLevel, const std::array<size_t, 3> &lowerBound,
                    const std::array<size_t, 3> &upperBound) {
    // can also use sorted cell optimization? need to be careful to sort only once per upper cell, to not sort for each
    // particle in the upper cell
    using namespace autopas::utils::ArrayMath::literals;
    using utils::ArrayUtils::operator<<;

    bool distanceCheck = checkDistance(upperLevel, lowerLevel);
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
    const double lowerInteractionLength = currentSkin() + _cutoffs[lowerLevel] / 2;
    const std::array<double, 3> lowerDir{lowerInteractionLength, lowerInteractionLength, lowerInteractionLength};
    size_t cellIndex1D;
    std::array<size_t, 3> particleCellIndex3D;
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
        startIndex3D = this->getMax(startIndex3D, lowerBound);
        stopIndex3D = this->getMin(stopIndex3D, upperBound);
      }
      if (distanceCheck) {
        particleCellIndex3D = lowerCB.get3DIndexOfPosition(pos);
      }
      for (size_t zl = startIndex3D[2]; zl <= stopIndex3D[2]; ++zl) {
        for (size_t yl = startIndex3D[1]; yl <= stopIndex3D[1]; ++yl) {
          cellIndex1D = autopas::utils::ThreeDimensionalMapping::threeToOneD(
              {static_cast<size_t>(startIndex3D[0]), yl, zl}, lowerLevelDims);
          for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; ++xl, ++cellIndex1D) {
            auto &lowerCell = lowerCB.getCell(cellIndex1D);
            if (lowerCell.isEmpty()) {
              continue;
            }
            // skip if cell is farther than interactionLength
            if (distanceCheck and
                this->getMinDistBetweenCellAndPointSquared(lowerCB, cellIndex1D, pos) > interactionLengthSquared) {
              continue;
            } else if (distanceCheck and xl >= particleCellIndex3D[0]) {
              break;
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
   * Traverses a single upper level cell and a lower level cells that are in the interaction range using SoA.
   * SoA is calculated between a single upper level particle and a full lower level cell.
   * @tparam Functor type of the functor
   * @param lowerCB lower level cell block
   * @param upperCB upper level cell block
   * @param upperCellCoords 3d cell index of cell belonging to the upper level
   * @param functor functor to apply
   * @param lowerLevel lower level index
   * @param lowerBound lower bound of the coordinates of lower level cells that contain owned particles
   * @param upperBound upper bound of the coordinates of lower level cells that contain owned particles
   * @param distanceCheck if true, checks if the minimum distance between upper level particle and lower level cell is
   * less than the interaction length, otherwise skips the lower level cell
   */
  template <class Functor>
  void SoATraversalParticleToCell(const CellBlock &lowerCB, const CellBlock &upperCB,
                                  const std::array<size_t, 3> upperCellCoords, Functor *functor, size_t lowerLevel,
                                  size_t upperLevel, const std::array<size_t, 3> &lowerBound,
                                  const std::array<size_t, 3> &upperBound) {
    using namespace autopas::utils::ArrayMath::literals;
    bool distanceCheck = checkDistance(upperLevel, lowerLevel);
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

    const double lowerInteractionLength = currentSkin() + _cutoffs[lowerLevel] / 2;
    const std::array<double, 3> lowerDir{lowerInteractionLength, lowerInteractionLength, lowerInteractionLength};
    size_t cellIndex1D;

    for (int idx = 0; idx < upperCell.size(); ++idx) {
      const std::array<double, 3> pos = {xptr[idx], yptr[idx], zptr[idx]};
      const double radius = upperCell[idx].getSize() / 2;
      const std::array<double, 3> interactionLength = lowerDir + radius;
      const double interactionLengthSquared = interactionLength[0] * interactionLength[0];
      auto soaSingleParticle = soa.constructView(idx, idx + 1);
      auto startIndex3D = lowerCB.get3DIndexOfPosition(pos - interactionLength);
      auto stopIndex3D = lowerCB.get3DIndexOfPosition(pos + interactionLength);
      if (containToOwnedOnly) {
        startIndex3D = this->getMax(startIndex3D, lowerBound);
        stopIndex3D = this->getMin(stopIndex3D, upperBound);
      }
      for (size_t zl = startIndex3D[2]; zl <= stopIndex3D[2]; ++zl) {
        for (size_t yl = startIndex3D[1]; yl <= stopIndex3D[1]; ++yl) {
          cellIndex1D = autopas::utils::ThreeDimensionalMapping::threeToOneD(
              {static_cast<size_t>(startIndex3D[0]), yl, zl}, lowerLevelDims);
          for (size_t xl = startIndex3D[0]; xl <= stopIndex3D[0]; ++xl, ++cellIndex1D) {
            auto &lowerCell = lowerCB.getCell(cellIndex1D);
            if (lowerCell.isEmpty()) {
              continue;
            }
            // skip if cell is farther than interactionLength
            if (distanceCheck and
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
    unsigned long smallestNumColors = 1e15;
    unsigned long largestBlocksPerColor = 0;
    std::array<size_t, 3> bestGroup = {1, 1, 1};
    for (size_t x_group = 1; x_group <= std::max(stride[0] - 1, static_cast<size_t>(1)); ++x_group)
      for (size_t y_group = 1; y_group <= std::max(stride[1] - 1, static_cast<size_t>(1)); ++y_group)
        for (size_t z_group = 1; z_group <= std::max(stride[2] - 1, static_cast<size_t>(1)); ++z_group) {
          std::array<size_t, 3> group = {x_group, y_group, z_group};
          std::array<size_t, 3> num_index{}, testStride{};
          for (size_t i = 0; i < 3; i++) {
            testStride[i] = 1 + static_cast<size_t>(std::ceil(1.0 * (stride[i] - 1) / group[i]));
            num_index[i] = (end[i] + (testStride[i] * group[i]) - 1) / (testStride[i] * group[i]);
          }
          const size_t numColors = testStride[0] * testStride[1] * testStride[2];
          const long numBlocksPerColor = num_index[0] * num_index[1] * num_index[2];
          if (numBlocksPerColor >= targetBlocksPerColor &&
              (numColors < smallestNumColors ||
               (numColors == smallestNumColors && numBlocksPerColor > largestBlocksPerColor))) {
            smallestNumColors = numColors;
            largestBlocksPerColor = numBlocksPerColor;
            bestGroup = group;
          }
        }
    return bestGroup;
  }

  /**
   * Finds the best group size for a given target number of blocks across all colors.
   * A block is a unit of work that is assigned to a thread in an OpenMP loop.
   * A group of blocks of 3d size x,y,z will be assigned to a thread per openmp loop iteration.
   * @param targetBlocks The target number of blocks.
   * @param stride The stride for each dimension.
   * @param end The end coordinates for each dimension.
   * @return The best group size for the given target number of blocks per color.
   */
  static std::array<size_t, 3> findBestGroupSizeForTargetBlocks(int targetBlocks, const std::array<size_t, 3> &stride,
                                                                const std::array<unsigned long, 3> &end) {
    unsigned long smallestCriterion = 1e15;
    std::array<size_t, 3> bestGroup = {1, 1, 1};
    for (size_t x_group = 1; x_group <= stride[0] * 2; ++x_group)
      for (size_t y_group = 1; y_group <= stride[1] * 2; ++y_group)
        for (size_t z_group = 1; z_group <= stride[2] * 2; ++z_group) {
          std::array<size_t, 3> group = {x_group, y_group, z_group};
          std::array<size_t, 3> num_index{}, testStride{};
          for (size_t i = 0; i < 3; i++) {
            testStride[i] = 1 + static_cast<size_t>(std::ceil(1.0 * (stride[i] - 1) / group[i]));
            num_index[i] = (end[i] + (testStride[i] * group[i]) - 1) / (testStride[i] * group[i]);
          }
          const long numColors = testStride[0] * testStride[1] * testStride[2];
          const long numBlocksPerColor = num_index[0] * num_index[1] * num_index[2];
          const unsigned long criterion = std::abs(numBlocksPerColor * numColors - targetBlocks);
          if (criterion < smallestCriterion) {
            smallestCriterion = criterion;
            bestGroup = group;
          }
        }
    return bestGroup;
  }
};
}  // namespace autopas
