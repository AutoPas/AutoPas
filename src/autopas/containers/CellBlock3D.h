/**
 * @file CellBlock3D.h
 *
 * @date 19 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>

#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace autopas::internal {
/**
 * Class that manages a block of ParticleCells.
 * It is used to resize the cell block and to handle the conversion of 3d to 1d indices.
 * @tparam ParticleCell Type of the handled ParticleCells.
 */
template <class ParticleCell>
class CellBlock3D : public CellBorderAndFlagManager {
 public:
  /**
   * The index type to access the particle cells
   */
  using index_t = std::size_t;
  /**
   * Constructor of CellBlock3D
   * @param vec Vector of ParticleCells that this class manages.
   * @param bMin Lower corner of the cell block.
   * @param bMax Higher corner of the cell block.
   * @param interactionLength Max. radius of interaction between particles.
   * @param cellSizeFactor Cell size factor relative to interactionLength.
   */
  CellBlock3D(std::vector<ParticleCell> &vec, const std::array<double, 3> &bMin, const std::array<double, 3> &bMax,
              double interactionLength, double cellSizeFactor = 1.0)
      : _cells(&vec), _boxMin(bMin), _boxMax(bMax), _interactionLength(interactionLength) {
    rebuild(vec, bMin, bMax, interactionLength, cellSizeFactor);

    for (int i = 0; i < 3; ++i) {
      if (bMax[i] < bMin[i] + interactionLength) {
        AutoPasLog(ERROR, "Interaction length too large is {}, bmin {}, bmax {}", interactionLength, bMin[i], bMax[i]);
        utils::ExceptionHandler::exception("Error in CellBlock3D: interaction Length too large!");
      }
    }
  }

  /**
   * Deleted copy constructor.
   */
  CellBlock3D(const CellBlock3D &) = delete;

  /**
   * Deleted assignment operator.
   * @return
   */
  CellBlock3D &operator=(const CellBlock3D) = delete;

  [[nodiscard]] bool cellCanContainHaloParticles(index_t index1d) const override {
    if (index1d < _firstOwnedCellIndex or index1d > _lastOwnedCellIndex) {
      return true;
    }
    const auto index3d = oneToThreeD(index1d);
    bool isHaloCell = false;
    for (size_t i = 0; i < index3d.size(); ++i) {
      if (index3d[i] < _cellsPerInteractionLength or
          index3d[i] >= _cellsPerDimensionWithHalo[i] - _cellsPerInteractionLength) {
        isHaloCell = true;
        break;
      }
    }
    return isHaloCell;
  }

  [[nodiscard]] bool cellCanContainOwnedParticles(index_t index1d) const override {
    return not cellCanContainHaloParticles(index1d);
  }

  /**
   * get the ParticleCell of a specified 1d index.
   * @param index1d The 1D index of the cell.
   * @return Non-const cell reference.
   */
  ParticleCell &getCell(index_t index1d) const;

  /**
   * Get the ParticleCell of a specified 3d index.
   * @param index3d The 3D index of the cell.
   * @return Non-const cell reference.
   */
  ParticleCell &getCell(const std::array<index_t, 3> &index3d) const;

  /**
   * Reserve memory for a given number of particles.
   * @param numParticles Particles incl. halo particles.
   */
  void reserve(size_t numParticles);

  /**
   * Rebuild the cell block. This resizes the cell block and sets the appropriate internal variables.
   * @param vec New vector of ParticleCells to which the internal pointer is set.
   * @param bMin New lower corner of the cell block.
   * @param bMax New higher corner of the cell block.
   * @param interactionLength Max. radius of interaction between particles.
   * @param cellSizeFactor Cell size factor relative to interactionLength.
   */
  void rebuild(std::vector<ParticleCell> &vec, const std::array<double, 3> &bMin, const std::array<double, 3> &bMax,
               double interactionLength, double cellSizeFactor);

  /**
   * Get the containing cell of a specified position.
   * @note If pos is outside the domain vector->operator[]() will trigger undefined behavior.
   *
   * @param pos The position for which the cell is needed.
   * @return Cell at the given position.
   */
  ParticleCell &getContainingCell(const std::array<double, 3> &pos) const;

  /**
   * Get the lower and upper corner of the cell at the 1d index index1d
   * @param index1d The 1d cell index.
   * @return std::pair of boxMin (lower corner) and boxMax (upper corner) of the box.
   */
  std::pair<std::array<double, 3>, std::array<double, 3>> getCellBoundingBox(index_t index1d) const;

  /**
   * Get the lower and upper corner of the cell at the 3d index index3d.
   * @param index3d The 3d cell index.
   * @return std::pair of boxMin (lower corner) and boxMax (upper corner) of the box.
   */
  std::pair<std::array<double, 3>, std::array<double, 3>> getCellBoundingBox(
      const std::array<index_t, 3> &index3d) const;

  /**
   * Get the 3d index of the cell block for a given position.
   * @note If pos is outside the domain the returned 3d index is also outside the domain.
   *
   * @param pos The position of interest.
   * @return The 3d cell index,
   */
  [[nodiscard]] std::array<index_t, 3> get3DIndexOfPosition(const std::array<double, 3> &pos) const;

  /**
   * Get the 1d index of the cell block for a given position.
   * @param pos the position of interest.
   * @return The 1d cell index.
   */
  [[nodiscard]] index_t get1DIndexOfPosition(const std::array<double, 3> &pos) const;

  /**
   * Get the dimension of the cell block including the halo boxes.
   * @return The dimensions of the cell block.
   */
  [[nodiscard]] const std::array<index_t, 3> &getCellsPerDimensionWithHalo() const {
    return _cellsPerDimensionWithHalo;
  }

  /**
   * Checks whether a given position is inside the halo region of the managed cell block.
   * @param position The given position.
   * @return true if the position is inside the halo region.
   */
  [[nodiscard]] bool checkInHalo(const std::array<double, 3> &position) const;

  /**
   * Deletes all particles in the halo cells of the managed cell block.
   */
  void clearHaloCells();

  /**
   * Invokes a callback on each halo cell that is at least partially within allowedDistance of the given position.
   * Always checks the cell containing position first (if it is a halo cell) before computing neighboring cell bounds.
   * Stops early and returns true as soon as func(cell) returns true.
   * @note The 1-norm (per-dimension Manhattan bounds) is used to determine nearby cells.
   *
   * @tparam Func Callable with signature bool(ParticleCell &).
   * @param position Cells close to this position are visited.
   * @param allowedDistance The maximal distance to the position.
   * @param func Callback invoked for each nearby halo cell; returning true terminates iteration early.
   * @return True iff func returned true for one of the visited halo cells.
   */
  template <typename Func>
  bool forEachNearbyHaloCell(const std::array<double, 3> &position, const double allowedDistance, Func &&func) const {
    using namespace autopas::utils::ArrayMath::literals;

    const auto index3D = get3DIndexOfPosition(position);
    const auto index1D = threeToOneD(index3D);

    // First, check the cell the particle is currently in if it is a halo cell.
    if ((*_cells)[index1D].getPossibleParticleOwnerships() == OwnershipState::halo) {
      if (func((*_cells)[index1D])) {
        return true;
      }
    }

    // Second, check neighboring halo cells within allowedDistance.
    // These indices are (already) at least 0 and at most _cellsPerDimensionWithHalo[i]-1.
    const auto lowIndex3D = get3DIndexOfPosition(position - allowedDistance);
    const auto highIndex3D = get3DIndexOfPosition(position + allowedDistance);

    auto currentIndex3D = lowIndex3D;
    for (currentIndex3D[0] = lowIndex3D[0]; currentIndex3D[0] <= highIndex3D[0]; ++currentIndex3D[0]) {
      for (currentIndex3D[1] = lowIndex3D[1]; currentIndex3D[1] <= highIndex3D[1]; ++currentIndex3D[1]) {
        for (currentIndex3D[2] = lowIndex3D[2]; currentIndex3D[2] <= highIndex3D[2]; ++currentIndex3D[2]) {
          // We have already checked the cell containing position, so skip it here.
          if (currentIndex3D == index3D) {
            continue;
          }
          const auto currentIndex1D = threeToOneD(currentIndex3D);
          if ((*_cells)[currentIndex1D].getPossibleParticleOwnerships() == OwnershipState::halo) {
            if (func((*_cells)[currentIndex1D])) {
              return true;
            }
          }
        }
      }
    }

    return false;
  }

  /**
   * Finds the next particle matching matchesParticle starting at (cellIndex, particleIndex).
   * When regionIter is true, restricts cell traversal to the 3D cell sub-box overlapping
   * [boxMin - verletSkin, boxMax + verletSkin].
   *
   * @tparam regionIter Whether iteration is restricted to a bounding box.
   * @tparam Predicate Callable with signature bool(const ParticleType &).
   * @param cellIndex Current 1D cell index (0 with particleIndex == 0 indicates start of iteration).
   * @param particleIndex Current particle index within the cell.
   * @param iteratorBehavior Ownership and sequentiality flags.
   * @param boxMin Lower corner of the search region (unused when regionIter is false).
   * @param boxMax Upper corner of the search region (unused when regionIter is false).
   * @param verletSkin Safety margin added around [boxMin, boxMax] when regionIter is true.
   * @param matchesParticle Predicate deciding whether a candidate particle fulfills iterator requirements.
   * @return Tuple of {particle pointer, cell index, particle index}, or {nullptr, 0, 0} if none remain.
   */
  template <bool regionIter, typename Predicate>
  std::tuple<const typename ParticleCell::ParticleType *, size_t, size_t> getNextParticle(
      size_t cellIndex, size_t particleIndex, IteratorBehavior iteratorBehavior, const std::array<double, 3> &boxMin,
      const std::array<double, 3> &boxMax, double verletSkin, Predicate &&matchesParticle) const {
    if (cellIndex >= _cells->size()) {
      return {nullptr, 0, 0};
    }

    const bool isStart = (cellIndex == 0 and particleIndex == 0);

    // Fast path: when advancing within an already-validated cell, scan its remaining particles directly without
    // recomputing cell bounds or 3D indices.
    if (not isStart) {
      if (const auto result = findMatchingParticleInCell(cellIndex, particleIndex, matchesParticle);
          std::get<0>(result) != nullptr) {
        return result;
      }
    }

    const bool includeOwned = static_cast<bool>(iteratorBehavior & IteratorBehavior::owned);
    const bool includeHalo = static_cast<bool>(iteratorBehavior & IteratorBehavior::halo);
    if (not includeOwned and not includeHalo) {
      return {nullptr, 0, 0};
    }

    const size_t stride = (iteratorBehavior & IteratorBehavior::forceSequential) ? 1 : autopas_get_num_threads();
    const size_t threadOffset = (iteratorBehavior & IteratorBehavior::forceSequential) ? 0 : autopas_get_thread_num();

    if constexpr (regionIter) {
      return findNextParticleInRegion(cellIndex, isStart, includeOwned, includeHalo, boxMin, boxMax, verletSkin, stride,
                                      threadOffset, matchesParticle);
    } else {
      return findNextParticleInBlock(cellIndex, isStart, includeOwned, includeHalo, stride, threadOffset,
                                     matchesParticle);
    }
  }

  /**
   * Get the lower corner of the halo region.
   * @return Coordinates of the lower corner.
   */
  [[nodiscard]] const std::array<double, 3> &getHaloBoxMin() const { return _haloBoxMin; }

  /**
   * Get the upper corner of the halo region.
   * @return Coordinates of the upper corner.
   */
  [[nodiscard]] const std::array<double, 3> &getHaloBoxMax() const { return _haloBoxMax; }

  /**
   * Get the number of cells per interaction length.
   * @return Cells per interaction length.
   */
  [[nodiscard]] unsigned long getCellsPerInteractionLength() const { return _cellsPerInteractionLength; }

  /**
   * Get the cell lengths.
   * @return Array of cell lengths.
   */
  [[nodiscard]] const std::array<double, 3> &getCellLength() const { return _cellLength; }

  /**
   * Get the number of cells in this block.
   * @return
   */
  index_t getNumCells() const;

  /**
   * 1D id of the first cell that is not in the halo.
   * @return
   */
  index_t getFirstOwnedCellIndex() const;

  /**
   * 1D id of the last cell before there are only halo cells left.
   * @return
   */
  index_t getLastOwnedCellIndex() const;

 private:
  /**
   * Scans a single cell starting at startParticleIndex for the first particle satisfying matchesParticle.
   * @tparam Predicate Callable with signature bool(const ParticleType &).
   * @param cellIndex 1D index of the cell to scan.
   * @param startParticleIndex First particle index to inspect.
   * @param matchesParticle Predicate deciding whether a particle matches.
   * @return Tuple of {particle pointer, cellIndex, particleIndex}, or {nullptr, 0, 0} if no particle matches.
   */
  template <typename Predicate>
  std::tuple<const typename ParticleCell::ParticleType *, size_t, size_t> findMatchingParticleInCell(
      size_t cellIndex, size_t startParticleIndex, Predicate &&matchesParticle) const {
    const auto &cell = (*_cells)[cellIndex];
    for (size_t pIdx = startParticleIndex; pIdx < cell.size(); ++pIdx) {
      if (matchesParticle(cell[pIdx])) {
        return {&cell[pIdx], cellIndex, pIdx};
      }
    }
    return {nullptr, 0, 0};
  }

  /**
   * Advances a region iterator to the next matching particle by stepping directly over the 3D cell sub-box overlapping
   * [boxMin - verletSkin, boxMax + verletSkin].
   */
  template <typename Predicate>
  std::tuple<const typename ParticleCell::ParticleType *, size_t, size_t> findNextParticleInRegion(
      size_t cellIndex, bool isStart, bool includeOwned, bool includeHalo, const std::array<double, 3> &boxMin,
      const std::array<double, 3> &boxMax, double verletSkin, size_t stride, size_t threadOffset,
      Predicate &&matchesParticle) const {
    using namespace autopas::utils::ArrayMath::literals;

    auto minIdx3D = get3DIndexOfPosition(boxMin - verletSkin);
    auto maxIdx3D = get3DIndexOfPosition(boxMax + verletSkin);

    // If only owned particles are requested, clamp the 3D sub-box directly to owned cells.
    if (not includeHalo) {
      for (size_t d = 0; d < 3; ++d) {
        minIdx3D[d] = std::max(minIdx3D[d], _cellsPerInteractionLength);
        maxIdx3D[d] = std::min(maxIdx3D[d], _cellsPerDimensionWithHalo[d] - _cellsPerInteractionLength - 1);
        if (minIdx3D[d] > maxIdx3D[d]) {
          return {nullptr, 0, 0};
        }
      }
    }

    const std::array<size_t, 3> subDims{
        maxIdx3D[0] - minIdx3D[0] + 1,
        maxIdx3D[1] - minIdx3D[1] + 1,
        maxIdx3D[2] - minIdx3D[2] + 1,
    };
    const size_t totalSubCells = subDims[0] * subDims[1] * subDims[2];

    size_t subIdx{};
    if (isStart) {
      subIdx = threadOffset;
    } else {
      const auto currentIdx3D = oneToThreeD(cellIndex);
      subIdx = ((currentIdx3D[2] - minIdx3D[2]) * subDims[1] + (currentIdx3D[1] - minIdx3D[1])) * subDims[0] +
               (currentIdx3D[0] - minIdx3D[0]) + stride;
    }

    const bool haloOnly = includeHalo and not includeOwned;
    for (; subIdx < totalSubCells; subIdx += stride) {
      const auto local3D = utils::ThreeDimensionalMapping::oneToThreeD(subIdx, subDims);
      const std::array<size_t, 3> cell3D{
          minIdx3D[0] + local3D[0],
          minIdx3D[1] + local3D[1],
          minIdx3D[2] + local3D[2],
      };
      if (haloOnly) {
        const bool isHaloCell = cell3D[0] < _cellsPerInteractionLength or
                                cell3D[0] >= _cellsPerDimensionWithHalo[0] - _cellsPerInteractionLength or
                                cell3D[1] < _cellsPerInteractionLength or
                                cell3D[1] >= _cellsPerDimensionWithHalo[1] - _cellsPerInteractionLength or
                                cell3D[2] < _cellsPerInteractionLength or
                                cell3D[2] >= _cellsPerDimensionWithHalo[2] - _cellsPerInteractionLength;
        if (not isHaloCell) {
          continue;
        }
      }
      const size_t nextCellIndex = threeToOneD(cell3D);
      if (const auto result = findMatchingParticleInCell(nextCellIndex, 0, matchesParticle);
          std::get<0>(result) != nullptr) {
        return result;
      }
    }
    return {nullptr, 0, 0};
  }

  /**
   * Advances a non-region iterator to the next matching particle across the cell block.
   */
  template <typename Predicate>
  std::tuple<const typename ParticleCell::ParticleType *, size_t, size_t> findNextParticleInBlock(
      size_t cellIndex, bool isStart, bool includeOwned, bool includeHalo, size_t stride, size_t threadOffset,
      Predicate &&matchesParticle) const {
    const size_t startCellIndex = includeHalo ? 0 : _firstOwnedCellIndex;
    const size_t endCellIndex = includeHalo ? (_cells->size() - 1) : _lastOwnedCellIndex;

    cellIndex = isStart ? (startCellIndex + threadOffset) : (cellIndex + stride);
    const bool checkCellOwnership = not(includeOwned and includeHalo);
    for (; cellIndex <= endCellIndex; cellIndex += stride) {
      if (checkCellOwnership) {
        const bool isRelevant = (includeOwned and cellCanContainOwnedParticles(cellIndex)) or
                                (includeHalo and cellCanContainHaloParticles(cellIndex));
        if (not isRelevant) {
          continue;
        }
      }
      if (const auto result = findMatchingParticleInCell(cellIndex, 0, matchesParticle);
          std::get<0>(result) != nullptr) {
        return result;
      }
    }
    return {nullptr, 0, 0};
  }

  /**
   * Transform a 1D index to a 3D index in the context of this cell block.
   * @param index1D
   * @return 3D cell index.
   */
  [[nodiscard]] std::array<index_t, 3> oneToThreeD(index_t index1D) const;

  /**
   * Transform a 3D index to a 1D index in the context of this cell block.
   * @param index3D
   * @return 1D cell index.
   */
  [[nodiscard]] index_t threeToOneD(const std::array<index_t, 3> &index3D) const;

  /**
   * Number of cells to be checked in each direction where we can find valid interaction partners.
   * This is also the number of Halo cells per direction (always symmetric in each dimension).
   */
  std::array<index_t, 3> _cellsPerDimensionWithHalo{};
  index_t _firstOwnedCellIndex{};
  index_t _lastOwnedCellIndex{};
  std::vector<ParticleCell> *_cells;

  std::array<double, 3> _boxMin, _boxMax;
  std::array<double, 3> _haloBoxMin{}, _haloBoxMax{};

  double _interactionLength{};

  unsigned long _cellsPerInteractionLength{};

  std::array<double, 3> _cellLength{};

  /**
   * 1.0 / _cellLength
   * Since this value is often needed for sorting particles in cells, we precompute it.
   */
  std::array<double, 3> _cellLengthReciprocal{};
};

template <class ParticleCell>
inline typename CellBlock3D<ParticleCell>::index_t CellBlock3D<ParticleCell>::get1DIndexOfPosition(
    const std::array<double, 3> &pos) const {
  return threeToOneD(get3DIndexOfPosition(pos));
}

template <class ParticleCell>
inline std::array<typename CellBlock3D<ParticleCell>::index_t, 3> CellBlock3D<ParticleCell>::get3DIndexOfPosition(
    const std::array<double, 3> &pos) const {
  std::array<typename CellBlock3D<ParticleCell>::index_t, 3> cellIndex{};

  for (size_t dim = 0; dim < 3; dim++) {
    // 0 <= cellIndex < cellsPerDimWithHalo
    const long int value = (static_cast<long int>(std::floor((pos[dim] - _boxMin[dim]) * _cellLengthReciprocal[dim]))) +
                           _cellsPerInteractionLength;
    const index_t nonNegativeValue = static_cast<index_t>(std::max(value, 0l));
    const index_t nonLargerValue = std::min(nonNegativeValue, _cellsPerDimensionWithHalo[dim] - 1);
    cellIndex[dim] = nonLargerValue;

    // sanity checks and precautions for numerical instabilities around the edge of the container box,
    // which would lead to doubling of particles
    if (pos[dim] >= _boxMax[dim]) {
      // pos is located outside the box. Make sure that cellIndex is at least in the positive halo range.
      cellIndex[dim] = std::max(cellIndex[dim], _cellsPerDimensionWithHalo[dim] - _cellsPerInteractionLength);
    } else if (pos[dim] < _boxMin[dim] and cellIndex[dim] == _cellsPerInteractionLength) {
      // pos is located outside but the index resolved to be the last IN the box. Fix by decrementing the index.
      --cellIndex[dim];
    } else if (pos[dim] < _boxMax[dim] and
               cellIndex[dim] == _cellsPerDimensionWithHalo[dim] - _cellsPerInteractionLength) {
      // pos is located inside the box, but cellIndex resolves to be the first halo cell. Fix by decrementing the index.
      --cellIndex[dim];
    }
  }

  return cellIndex;
}

template <class ParticleCell>
void CellBlock3D<ParticleCell>::reserve(size_t numParticles) {
  const auto particlesPerCell = numParticles / _cells->size();
  for (auto &cell : *_cells) {
    cell.reserve(particlesPerCell);
  }
}

template <class ParticleCell>
inline void CellBlock3D<ParticleCell>::rebuild(std::vector<ParticleCell> &vec, const std::array<double, 3> &bMin,
                                               const std::array<double, 3> &bMax, double interactionLength,
                                               double cellSizeFactor) {
  using namespace autopas::utils::ArrayMath::literals;

  _cells = &vec;
  _boxMin = bMin;
  _boxMax = bMax;
  _interactionLength = interactionLength;

  if (cellSizeFactor >= 1.0) {
    _cellsPerInteractionLength = 1;
  } else {
    _cellsPerInteractionLength = ceil(1.0 / cellSizeFactor);
  }
  // compute cell length and number of cells
  index_t numCells = 1;
  for (int d = 0; d < 3; ++d) {
    const double boxLength = _boxMax[d] - _boxMin[d];
    // The number of cells is rounded down because the cells will be stretched to fit.
    // std::max to ensure there is at least one cell.
    const auto cellsPerDim =
        std::max(static_cast<index_t>(std::floor(boxLength / (_interactionLength * cellSizeFactor))), 1ul);

    _cellsPerDimensionWithHalo[d] = cellsPerDim + 2 * _cellsPerInteractionLength;

    _cellLength[d] = boxLength / static_cast<double>(cellsPerDim);

    _cellLengthReciprocal[d] = static_cast<double>(cellsPerDim) / boxLength;  // compute with least rounding possible

    _haloBoxMin[d] = _boxMin[d] - _cellsPerInteractionLength * _cellLength[d];
    _haloBoxMax[d] = _boxMax[d] + _cellsPerInteractionLength * _cellLength[d];

    numCells *= _cellsPerDimensionWithHalo[d];
  }
  AutoPasLog(TRACE, "Box Length incl Halo : {}", autopas::utils::ArrayUtils::to_string(_haloBoxMax - _haloBoxMin));
  AutoPasLog(TRACE, "Cells/Dim  incl Halo : {}", autopas::utils::ArrayUtils::to_string(_cellsPerDimensionWithHalo));
  AutoPasLog(TRACE, "Cell Length          : {}", autopas::utils::ArrayUtils::to_string(_cellLength));
  AutoPasLog(TRACE, "Interaction Length   : {}", interactionLength);
  AutoPasLog(TRACE, "Cell Size Factor     : {}", cellSizeFactor);

  _firstOwnedCellIndex = _cellsPerDimensionWithHalo[0] * _cellsPerDimensionWithHalo[1] * _cellsPerInteractionLength +
                         _cellsPerDimensionWithHalo[0] * _cellsPerInteractionLength + _cellsPerInteractionLength;
  _lastOwnedCellIndex = numCells - 1 - _firstOwnedCellIndex;

  // initialize cells
  _cells->resize(numCells);

  for (auto &cell : *_cells) {
    cell.setCellLength(_cellLength);
  }

  // determine the OwnershipStates, each cell can contain. This is later used in the CellFunctor to skip calculations
  for (int i = 0; i < numCells; i++) {
    const bool canHaveHalos = cellCanContainHaloParticles(i);
    const bool canHaveOwned = cellCanContainOwnedParticles(i);
    if (canHaveHalos and canHaveOwned) {
      (*_cells)[i].setPossibleParticleOwnerships(OwnershipState::owned | OwnershipState::halo);
    } else if (canHaveHalos) {
      (*_cells)[i].setPossibleParticleOwnerships(OwnershipState::halo);
    } else if (canHaveOwned) {
      (*_cells)[i].setPossibleParticleOwnerships(OwnershipState::owned);
    } else {
      (*_cells)[i].setPossibleParticleOwnerships(OwnershipState::dummy);
    }
  }
}

template <class ParticleCell>
typename CellBlock3D<ParticleCell>::index_t CellBlock3D<ParticleCell>::getNumCells() const {
  return _cells->size();
}

template <class ParticleCell>
typename CellBlock3D<ParticleCell>::index_t CellBlock3D<ParticleCell>::getFirstOwnedCellIndex() const {
  return _firstOwnedCellIndex;
}

template <class ParticleCell>
typename CellBlock3D<ParticleCell>::index_t CellBlock3D<ParticleCell>::getLastOwnedCellIndex() const {
  return _lastOwnedCellIndex;
}

template <class ParticleCell>
inline std::pair<std::array<double, 3>, std::array<double, 3>> CellBlock3D<ParticleCell>::getCellBoundingBox(
    const index_t index1d) const {
  return this->getCellBoundingBox(this->oneToThreeD(index1d));
}

template <class ParticleCell>
inline std::pair<std::array<double, 3>, std::array<double, 3>> CellBlock3D<ParticleCell>::getCellBoundingBox(
    const std::array<index_t, 3> &index3d) const {
  std::array<double, 3> boxmin{}, boxmax{};
  for (int d = 0; d < 3; d++) {
    // defaults
    boxmin[d] = index3d[d] * this->_cellLength[d] + _haloBoxMin[d];

    // stupid rounding errors. Snap values to the exact box values.
    if (index3d[d] == 0) {
      boxmin[d] = _haloBoxMin[d];
    } else if (index3d[d] == _cellsPerInteractionLength) {
      // Case: we are at the lower boundary of the non-halo box
      boxmin[d] = _boxMin[d];
    }

    boxmax[d] = boxmin[d] + this->_cellLength[d];

    // This must not be an else to the if block above as both cases can be true.
    // e.g. if there is only one cell
    if (index3d[d] == this->_cellsPerDimensionWithHalo[d] - _cellsPerInteractionLength - 1) {
      boxmax[d] = _boxMax[d];
    } else if (index3d[d] == this->_cellsPerDimensionWithHalo[d] - 1) {
      boxmin[d] = _haloBoxMax[d] - this->_cellLength[d];
      boxmax[d] = _haloBoxMax[d];
    }
  }
  return {boxmin, boxmax};
}

template <class ParticleCell>
inline ParticleCell &CellBlock3D<ParticleCell>::getContainingCell(const std::array<double, 3> &pos) const {
  auto ind = get1DIndexOfPosition(pos);
  return getCell(ind);
}

template <class ParticleCell>
inline ParticleCell &CellBlock3D<ParticleCell>::getCell(index_t index1d) const {
  return (*_cells)[index1d];
}

template <class ParticleCell>
inline ParticleCell &CellBlock3D<ParticleCell>::getCell(const std::array<index_t, 3> &index3d) const {
  return (*_cells)[threeToOneD(index3d)];
}

template <class ParticleCell>
inline std::array<typename CellBlock3D<ParticleCell>::index_t, 3> CellBlock3D<ParticleCell>::oneToThreeD(
    index_t index1D) const {
  return utils::ThreeDimensionalMapping::oneToThreeD(index1D, _cellsPerDimensionWithHalo);
}

template <class ParticleCell>
inline typename CellBlock3D<ParticleCell>::index_t CellBlock3D<ParticleCell>::threeToOneD(
    const std::array<index_t, 3> &index3D) const {
  return utils::ThreeDimensionalMapping::threeToOneD(index3D, _cellsPerDimensionWithHalo);
}

template <class ParticleCell>
bool CellBlock3D<ParticleCell>::checkInHalo(const std::array<double, 3> &position) const {
  return autopas::utils::inBox(position, _haloBoxMin, _haloBoxMax) &&
         autopas::utils::notInBox(position, _boxMin, _boxMax);
}

template <class ParticleCell>
void CellBlock3D<ParticleCell>::clearHaloCells() {
  std::vector<index_t> haloSlices(2 * _cellsPerInteractionLength);
  auto mid = haloSlices.begin() + _cellsPerInteractionLength;
  std::iota(haloSlices.begin(), mid, 0);

  // x: min and max of x
  std::iota(mid, haloSlices.end(), _cellsPerDimensionWithHalo[0] - _cellsPerInteractionLength);
  for (index_t i : haloSlices) {
    for (index_t j = 0; j < _cellsPerDimensionWithHalo[1]; j++) {
      for (index_t k = 0; k < _cellsPerDimensionWithHalo[2]; k++) {
        index_t index = threeToOneD({i, j, k});
        (*_cells)[index].clear();
      }
    }
  }
  // y: min and max of y
  std::iota(mid, haloSlices.end(), _cellsPerDimensionWithHalo[1] - _cellsPerInteractionLength);
  for (index_t i = 1; i < _cellsPerDimensionWithHalo[0] - 1; i++) {  // 0 and cells-1 already done in previous loop
    for (index_t j : haloSlices) {
      for (index_t k = 0; k < _cellsPerDimensionWithHalo[2]; k++) {
        index_t index = threeToOneD({i, j, k});
        (*_cells)[index].clear();
      }
    }
  }
  // z: min and max of z
  std::iota(mid, haloSlices.end(), _cellsPerDimensionWithHalo[2] - _cellsPerInteractionLength);
  for (index_t i = 1; i < _cellsPerDimensionWithHalo[0] - 1; i++) {    // 0 and cells-1 already done in previous loop
    for (index_t j = 1; j < _cellsPerDimensionWithHalo[1] - 1; j++) {  // 0 and cells-1 already done in previous loop
      for (index_t k : haloSlices) {
        index_t index = threeToOneD({i, j, k});
        (*_cells)[index].clear();
      }
    }
  }
}
}  // namespace autopas::internal
