/**
 * @file CellBlock3D.h
 *
 * @date 19 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <cmath>
#include <numeric>
#include <vector>

#include "autopas/containers/CellBorderAndFlagManager.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/inBox.h"

namespace autopas::internal {
/**
 * Class that manages a block of ParticleCells.
 * It is used to resize the cellblock and to handle the conversion of 3d to 1d
 * indices
 * @tparam ParticleCell type of the handled ParticleCells
 */
template <class ParticleCell>
class CellBlock3D : public CellBorderAndFlagManager {
 public:
  /**
   * the index type to access the particle cells
   */
  using index_t = std::size_t;
  /**
   * Constructor of CellBlock3D
   * @param vec vector of ParticleCells that this class manages
   * @param bMin lower corner of the cellblock
   * @param bMax higher corner of the cellblock
   * @param interactionLength max. radius of interaction between particles
   * @param cellSizeFactor cell size factor relative to interactionLength
   */
  CellBlock3D(std::vector<ParticleCell> &vec, const std::array<double, 3> bMin, const std::array<double, 3> bMax,
              double interactionLength, double cellSizeFactor = 1.0) {
    rebuild(vec, bMin, bMax, interactionLength, cellSizeFactor);

    for (int i = 0; i < 3; ++i) {
      if (bMax[i] < bMin[i] + interactionLength) {
        AutoPasLog(error, "Interaction length too large is {}, bmin {}, bmax {}", interactionLength, bMin[i], bMax[i]);
        utils::ExceptionHandler::exception("Error in CellBlock3D: interaction Length too large!");
      }
    }
  }

  /**
   * deleted copy constructor
   */
  CellBlock3D(const CellBlock3D &) = delete;

  /**
   * delete assignment operator
   * @return deleted
   */
  CellBlock3D &operator=(const CellBlock3D) = delete;

  [[nodiscard]] bool cellCanContainHaloParticles(index_t index1d) const override {
    auto index3d = index3D(index1d);
    bool isHaloCell = false;
    for (size_t i = 0; i < 3; i++) {
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
   * Checks if cell with index1d can be ignored for iteration with currently selected behavior.
   * @param index1d 1d index of checked cell
   * @param behavior @see IteratorBehavior
   * @return false if this cell can contain particles that would be affected by current behavior
   */
  [[nodiscard]] bool ignoreCellForIteration(index_t index1d, IteratorBehavior behavior) const {
    if ((behavior & IteratorBehavior::halo) and cellCanContainHaloParticles(index1d)) {
      return false;
    }
    if ((behavior & IteratorBehavior::owned) and cellCanContainOwnedParticles(index1d)) {
      return false;
    }
    if (behavior & IteratorBehavior::dummy) {
      return false;
    }
    return true;
  }

  /**
   * get the ParticleCell of a specified 1d index
   * @param index1d the index of the cell
   * @return the specified cell
   */
  ParticleCell &getCell(index_t index1d) const;

  /**
   * get the ParticleCell of a specified 3d index
   * @param index3d the index of the cell
   * @return the specified cell
   */
  ParticleCell &getCell(const std::array<index_t, 3> &index3d) const;

  /**
   * rebuild the cellblock. This resizes the cellblock and sets the appropriate
   * internal variables
   * @param vec new vector of ParticleCells to which the internal pointer is set
   * @param bMin new lower corner of the cellblock
   * @param bMax new higher corner of the cellblock
   * @param interactionLength max. radius of interaction between particles
   * @param cellSizeFactor cell size factor relative to interactionLength
   */
  void rebuild(std::vector<ParticleCell> &vec, const std::array<double, 3> &bMin, const std::array<double, 3> &bMax,
               double interactionLength, double cellSizeFactor);

  // this class doesn't actually know about particles
  /**
   * Get the containing cell of a specified position
   * @param pos the position for which the cell is needed
   * @return cell at the given position
   */
  ParticleCell &getContainingCell(const std::array<double, 3> &pos) const;

  /**
   * Get the lower and upper corner of the cell at the 1d index index1d
   * @param index1d the 1d index
   * @return std::pair of boxmin (lower corner) and boxmax (upper corner) of the box.
   */
  std::pair<std::array<double, 3>, std::array<double, 3>> getCellBoundingBox(index_t index1d) const;

  /**
   * Get the lower and upper corner of the cell at the 3d index index3d
   * @param index3d the 3d index
   * @return std::pair of boxmin (lower corner) and boxmax (upper corner) of the box.
   */
  std::pair<std::array<double, 3>, std::array<double, 3>> getCellBoundingBox(
      const std::array<index_t, 3> &index3d) const;

  /**
   * get the 3d index of the cellblock for a given position
   * @param pos the position
   * @return the 3d index
   */
  [[nodiscard]] std::array<index_t, 3> get3DIndexOfPosition(const std::array<double, 3> &pos) const;

  /**
   * get the 1d index of the cellblock for a given position
   * @param pos the position
   * @return the 1d index
   */
  [[nodiscard]] index_t get1DIndexOfPosition(const std::array<double, 3> &pos) const;

  /**
   * get the dimension of the cellblock including the haloboxes
   * @return the dimensions of the cellblock
   */
  [[nodiscard]] const std::array<index_t, 3> &getCellsPerDimensionWithHalo() const {
    return _cellsPerDimensionWithHalo;
  }

  /**
   * checks whether a given position is inside the halo region of the managed
   * cell block
   * @param position the given position
   * @return true if the position is inside the halo region
   */
  [[nodiscard]] bool checkInHalo(const std::array<double, 3> &position) const;

  /**
   * deletes all particles in the halo cells of the managed cell block
   */
  void clearHaloCells();

  /**
   * Get the nearby halo cells.
   * A list of halo cells is returned whose distance to position is at most
   * allowedDistance. If position is inside a halo cell that cell is also
   * returned. (1 norm is used, i.e. the distance is computed for each dimension
   * separately, Manhattan distance)
   * @param position cells close to this position are to be returned
   * @param allowedDistance the maximal distance to the position
   * @return a container of references to nearby halo cells
   */
  std::vector<ParticleCell *> getNearbyHaloCells(const std::array<double, 3> &position, double allowedDistance) const {
    auto index3d = get3DIndexOfPosition(position);

    std::vector<ParticleCell *> closeHaloCells;

    auto isHaloCell = [&](const auto &index) {
      for (size_t i = 0; i < 3; ++i) {
        if (index[i] < _cellsPerInteractionLength or
            index[i] >= _cellsPerDimensionWithHalo[i] - _cellsPerInteractionLength) {
          return true;
        }
      }
      return false;
    };

    // always add the cell the particle is currently in first, for that we test if it is a halo cell.
    if (isHaloCell(index3d)) {
      closeHaloCells.push_back(&getCell(index3d));
    }

    auto lowIndex3D = get3DIndexOfPosition(utils::ArrayMath::subScalar(position, allowedDistance));
    auto highIndex3D = get3DIndexOfPosition(utils::ArrayMath::addScalar(position, allowedDistance));
    // these indices are (already) at least 0 and at most _cellsPerDimensionWithHalo[i]-1

    auto currentIndex = lowIndex3D;

    for (currentIndex[0] = lowIndex3D[0]; currentIndex[0] <= highIndex3D[0]; ++currentIndex[0]) {
      for (currentIndex[1] = lowIndex3D[1]; currentIndex[1] <= highIndex3D[1]; ++currentIndex[1]) {
        for (currentIndex[2] = lowIndex3D[2]; currentIndex[2] <= highIndex3D[2]; ++currentIndex[2]) {
          // we have already added the cell which normally would belong to the particle, so skip here:
          if (currentIndex == index3d) {
            continue;
          }
          // we need to return the cell is it is a halo cell.
          if (isHaloCell(currentIndex)) {
            closeHaloCells.push_back(&getCell(currentIndex));
          }
        }
      }
    }

    return closeHaloCells;
  }

  /**
   * Get the lower corner of the halo region.
   * @return Coordinates of the lower corner.
   */
  [[nodiscard]] std::array<double, 3> getHaloBoxMin() const { return _haloBoxMin; }

  /**
   * Get the upper corner of the halo region.
   * @return Coordinates of the upper corner.
   */
  [[nodiscard]] std::array<double, 3> getHaloBoxMax() const { return _haloBoxMax; }

  /**
   * Get the number of cells per interaction length.
   * @return cells per interaction length.
   */
  [[nodiscard]] unsigned long getCellsPerInteractionLength() const { return _cellsPerInteractionLength; }

  /**
   * Get the cell lengths.
   * @return array of cell lengths.
   */
  [[nodiscard]] const std::array<double, 3> &getCellLength() const { return _cellLength; }

 private:
  [[nodiscard]] std::array<index_t, 3> index3D(index_t index1d) const;

  [[nodiscard]] index_t index1D(const std::array<index_t, 3> &index3d) const;

  std::array<index_t, 3> _cellsPerDimensionWithHalo;
  index_t _numCells;
  std::vector<ParticleCell> *_cells;

  std::array<double, 3> _boxMin, _boxMax;
  std::array<double, 3> _haloBoxMin, _haloBoxMax;

  double _interactionLength;

  unsigned long _cellsPerInteractionLength;

  std::array<double, 3> _cellLength;

  // 1 over above. Since this value is needed for sorting particles in cells, it
  // is computed quite often
  std::array<double, 3> _cellLengthReciprocal;
};

template <class ParticleCell>
inline typename CellBlock3D<ParticleCell>::index_t CellBlock3D<ParticleCell>::get1DIndexOfPosition(
    const std::array<double, 3> &pos) const {
  return index1D(get3DIndexOfPosition(pos));
}

template <class ParticleCell>
inline std::array<typename CellBlock3D<ParticleCell>::index_t, 3> CellBlock3D<ParticleCell>::get3DIndexOfPosition(
    const std::array<double, 3> &pos) const {
  std::array<typename CellBlock3D<ParticleCell>::index_t, 3> cellIndex{};

  for (size_t dim = 0; dim < 3; dim++) {
    const long int value = (static_cast<long int>(std::floor((pos[dim] - _boxMin[dim]) * _cellLengthReciprocal[dim]))) +
                           _cellsPerInteractionLength;
    const index_t nonNegativeValue = static_cast<index_t>(std::max(value, 0l));
    const index_t nonLargerValue = std::min(nonNegativeValue, _cellsPerDimensionWithHalo[dim] - 1);
    cellIndex[dim] = nonLargerValue;

    // this is a sanity check to prevent doubling of particles
    if (pos[dim] >= _boxMax[dim]) {
      // pos[dim] is located outside of the box. Make sure that cellIndex also references cell outside of the box.
      cellIndex[dim] = std::max(cellIndex[dim], _cellsPerDimensionWithHalo[dim] - _cellsPerInteractionLength);
    } else if (pos[dim] < _boxMin[dim] && cellIndex[dim] == _cellsPerInteractionLength) {
      // pos[dim] is located below the box (outside), but cellIndex references cell inside of the box. Correct cellIndex
      // by subtracting 1.
      --cellIndex[dim];
    } else if (pos[dim] < _boxMax[dim] &&
               cellIndex[dim] == _cellsPerDimensionWithHalo[dim] - _cellsPerInteractionLength) {
      // pos[dim] is located inside of the box, but cellIndex references cell outside of the box. Correct cellIndex to
      // last cell inside of the box.
      cellIndex[dim] = _cellsPerDimensionWithHalo[dim] - _cellsPerInteractionLength - 1;
    }
  }

  return cellIndex;
  // in very rare cases rounding is stupid, thus we need a check...
  /// @todo when the border and flag manager is there
}

template <class ParticleCell>
inline void CellBlock3D<ParticleCell>::rebuild(std::vector<ParticleCell> &vec, const std::array<double, 3> &bMin,
                                               const std::array<double, 3> &bMax, double interactionLength,
                                               double cellSizeFactor) {
  _cells = &vec;
  _boxMin = bMin;
  _boxMax = bMax;
  _interactionLength = interactionLength;

  if (cellSizeFactor >= 1.0) {
    _cellsPerInteractionLength = 1;
  } else {
    _cellsPerInteractionLength = ceil(1.0 / cellSizeFactor);
  }
  // compute cell length
  _numCells = 1;
  for (int d = 0; d < 3; ++d) {
    const double diff = _boxMax[d] - _boxMin[d];
    auto cellsPerDim = static_cast<index_t>(std::floor(diff / (_interactionLength * cellSizeFactor)));
    // at least one central cell
    cellsPerDim = std::max(cellsPerDim, 1ul);

    _cellsPerDimensionWithHalo[d] = cellsPerDim + 2 * _cellsPerInteractionLength;

    _cellLength[d] = diff / cellsPerDim;

    _cellLengthReciprocal[d] = cellsPerDim / diff;  // compute with least rounding possible

    _haloBoxMin[d] = _boxMin[d] - _cellsPerInteractionLength * _cellLength[d];
    _haloBoxMax[d] = _boxMax[d] + _cellsPerInteractionLength * _cellLength[d];

    _numCells *= _cellsPerDimensionWithHalo[d];
  }
  AutoPasLog(trace, "Box Length incl Halo : {}", autopas::utils::ArrayUtils::to_string(autopas::utils::ArrayMath::sub(_haloBoxMax, _haloBoxMin)));
  AutoPasLog(trace, "Cells/Dim  incl Halo : {}", autopas::utils::ArrayUtils::to_string(_cellsPerDimensionWithHalo));
  AutoPasLog(trace, "Cell Length          : {}", autopas::utils::ArrayUtils::to_string(_cellLength));
  AutoPasLog(trace, "Interaction Length   : {}", interactionLength);
  AutoPasLog(trace, "Cell Size Factor     : {}", cellSizeFactor);

  _cells->resize(_numCells);
  for (auto &cell : *_cells) {
    cell.setCellLength(_cellLength);
  }
}

template <class ParticleCell>
inline std::pair<std::array<double, 3>, std::array<double, 3>> CellBlock3D<ParticleCell>::getCellBoundingBox(
    const index_t index1d) const {
  return this->getCellBoundingBox(this->index3D(index1d));
}

template <class ParticleCell>
inline std::pair<std::array<double, 3>, std::array<double, 3>> CellBlock3D<ParticleCell>::getCellBoundingBox(
    const std::array<index_t, 3> &index3d) const {
  std::array<double, 3> boxmin{}, boxmax{};
  for (int d = 0; d < 3; d++) {
    // defaults
    boxmin[d] = index3d[d] * this->_cellLength[d] + _haloBoxMin[d];
    boxmax[d] = (index3d[d] + 1) * this->_cellLength[d] + _haloBoxMin[d];

    // stupid rounding errors (makes sure that the lower corner is set
    // correctly!
    if (index3d[d] == 0) {
      boxmin[d] = _haloBoxMin[d];
      boxmax[d] = this->_cellLength[d];
    } else if (index3d[d] == _cellsPerInteractionLength) {
      boxmin[d] = _boxMin[d];
    }
    // no else!, as this might ALSO be 1
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
  return _cells->at(index1d);
}

template <class ParticleCell>
inline ParticleCell &CellBlock3D<ParticleCell>::getCell(const std::array<index_t, 3> &index3d) const {
  return _cells->at(index1D(index3d));
}

template <class ParticleCell>
inline std::array<typename CellBlock3D<ParticleCell>::index_t, 3> CellBlock3D<ParticleCell>::index3D(
    index_t index1d) const {
  return utils::ThreeDimensionalMapping::oneToThreeD(index1d, _cellsPerDimensionWithHalo);
}

template <class ParticleCell>
inline typename CellBlock3D<ParticleCell>::index_t CellBlock3D<ParticleCell>::index1D(
    const std::array<index_t, 3> &index3d) const {
  return utils::ThreeDimensionalMapping::threeToOneD(index3d, _cellsPerDimensionWithHalo);
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
        index_t index = index1D({i, j, k});
        (*_cells)[index].clear();
      }
    }
  }
  // y: min and max of y
  std::iota(mid, haloSlices.end(), _cellsPerDimensionWithHalo[1] - _cellsPerInteractionLength);
  for (index_t i = 1; i < _cellsPerDimensionWithHalo[0] - 1; i++) {  // 0 and cells-1 already done in previous loop
    for (index_t j : haloSlices) {
      for (index_t k = 0; k < _cellsPerDimensionWithHalo[2]; k++) {
        index_t index = index1D({i, j, k});
        (*_cells)[index].clear();
      }
    }
  }
  // z: min and max of z
  std::iota(mid, haloSlices.end(), _cellsPerDimensionWithHalo[2] - _cellsPerInteractionLength);
  for (index_t i = 1; i < _cellsPerDimensionWithHalo[0] - 1; i++) {    // 0 and cells-1 already done in previous loop
    for (index_t j = 1; j < _cellsPerDimensionWithHalo[1] - 1; j++) {  // 0 and cells-1 already done in previous loop
      for (index_t k : haloSlices) {
        index_t index = index1D({i, j, k});
        (*_cells)[index].clear();
      }
    }
  }
}
}  // namespace autopas::internal
