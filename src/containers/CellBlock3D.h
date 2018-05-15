/*
 * CellBlock3D.h
 *
 *  Created on: 19 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_CELLBLOCK3D_H_
#define SRC_CONTAINERS_CELLBLOCK3D_H_

#include <array>
#include <cmath>
#include <vector>
#include "CellBoarderAndFlagManager.h"
#include "utils/ExceptionHandler.h"
#include "utils/ThreeDimensionalMapping.h"
#include "utils/inBox.h"

namespace autopas {
/**
 * Class that manages a block of ParticleCells.
 * It is used to resize the cellblock and to handle the conversion of 3d to 1d
 * indices
 * @tparam ParticleCell type of the handled ParticleCells
 */
template <class ParticleCell>
class CellBlock3D : public CellBoarderAndFlagManager {
 public:
  /**
   * the index type to access the particle cells
   */
  typedef std::size_t index_t;
  /**
   * Constructor of CellBlock3D
   * @param vec vector of ParticleCells that this class manages
   * @param bMin lower corner of the cellblock
   * @param bMax higher corner of the cellblock
   * @param interactionLength minimal size of each cell
   */
  CellBlock3D(std::vector<ParticleCell> &vec, const std::array<double, 3> bMin,
              const std::array<double, 3> bMax, double interactionLength) {
    rebuild(vec, bMin, bMax, interactionLength);
    for (int i = 0; i < 3; ++i) {
      if (bMax[i] < bMin[i] + interactionLength) {
        AutoPasLogger->error(
            "CellBlock3D: interaction length too large is {}, bmin {}, bmax {}",
            interactionLength, bMin[i], bMax[i]);
        utils::ExceptionHandler::exception(
            "Error in CellBlock3D: interaction Length too large!");
      }
    }
  }

  bool isHaloCell(index_t index1d) const override {
    auto index3d = index3D(index1d);
    bool isHaloCell = false;
    for (size_t i = 0; i < 3; i++) {
      if (index3d[i] == 0 or index3d[i] == _cellsPerDimensionWithHalo[i] - 1) {
        isHaloCell = true;
      }
    }
    return isHaloCell;
  }

  bool isOwningCell(index_t index1d) const override {
    return not isHaloCell(index1d);
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
   * @param interactionLength new minimal size of each cell
   */
  void rebuild(std::vector<ParticleCell> &vec,
               const std::array<double, 3> &bMin,
               const std::array<double, 3> &bMax, double interactionLength);

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
   * @param boxmin the lower corner (out)
   * @param boxmax the upper corner (out)
   */
  void getCellBoundingBox(index_t index1d, std::array<double, 3> &boxmin,
                          std::array<double, 3> &boxmax);

  /**
   * Get the lower and upper corner of the cell at the 3d index index3d
   * @param index3d the 3d index
   * @param boxmin the lower corner (out)
   * @param boxmax the upper corner (out)
   */
  void getCellBoundingBox(const std::array<index_t, 3> &index3d,
                          std::array<double, 3> &boxmin,
                          std::array<double, 3> &boxmax);

  /**
   * get the 3d index of the cellblock for a given position
   * @param pos the position
   * @return the 3d index
   */
  std::array<index_t, 3> get3DIndexOfPosition(
      const std::array<double, 3> &pos) const;

  /**
   * get the 1d index of the cellblock for a given position
   * @param pos the position
   * @return the 1d index
   */
  index_t get1DIndexOfPosition(const std::array<double, 3> &pos) const;

  /**
   * get the dimension of the cellblock including the haloboxes
   * @return the dimensions of the cellblock
   */
  const std::array<index_t, 3> &getCellsPerDimensionWithHalo() const {
    return _cellsPerDimensionWithHalo;
  }

  /**
   * checks whether a given position is inside the halo region of the managed
   * cell block
   * @param position the given position
   * @return true if the position is inside the halo region
   */
  bool checkInHalo(std::array<double, 3> position) const;

  /**
   * deletes all particles in the halo cells of the managed cell block
   */
  void clearHaloCells() {
    // x: min and max of x
    for (index_t i :
         {static_cast<index_t>(0),
          static_cast<index_t>(_cellsPerDimensionWithHalo[0] - 1)}) {
      for (index_t j = 0; j < _cellsPerDimensionWithHalo[1]; j++) {
        for (index_t k = 0; k < _cellsPerDimensionWithHalo[2]; k++) {
          index_t index = index1D({i, j, k});
          (*_vec1D)[index].clear();
        }
      }
    }
    // y: min and max of y
    for (index_t i = 1; i < _cellsPerDimensionWithHalo[0] - 1;
         i++) {  // 0 and cells-1 already done in previous loop
      for (index_t j :
           {static_cast<index_t>(0),
            static_cast<index_t>(_cellsPerDimensionWithHalo[1] - 1)}) {
        for (index_t k = 0; k < _cellsPerDimensionWithHalo[2]; k++) {
          index_t index = index1D({i, j, k});
          (*_vec1D)[index].clear();
        }
      }
    }
    // z: min and max of z
    for (index_t i = 1; i < _cellsPerDimensionWithHalo[0] - 1;
         i++) {  // 0 and cells-1 already done in previous loop
      for (index_t j = 1; j < _cellsPerDimensionWithHalo[1] - 1;
           j++) {  // 0 and cells-1 already done in previous loop
        for (index_t k :
             {static_cast<index_t>(0),
              static_cast<index_t>(_cellsPerDimensionWithHalo[2] - 1)}) {
          index_t index = index1D({i, j, k});
          (*_vec1D)[index].clear();
        }
      }
    }
  }

  /**
   * Get the nearby halo cells.
   * A list of halo cells is returned whose distance to position is at most
   * allowedDistance. If position is inside a halo cell that cell is also
   * returned. (1 norm is used, i.e. the distance is computed for each dimension
   * separately)
   * @param position cells close to this position are to be returned
   * @param allowedDistance the maximal distance to the position
   * @return a container of references to nearby halo cells
   */
  std::vector<ParticleCell *> getNearbyHaloCells(std::array<double, 3> position,
                                                 double allowedDistance) {
    auto index3d = get3DIndexOfPosition(position);
    std::array<int, 3> diff = {0, 0, 0};
    auto currentIndex = index3d;
    std::vector<ParticleCell *> closeHaloCells;
    for (diff[0] = -1; diff[0] < 2; diff[0]++) {
      currentIndex[0] = index3d[0] + diff[0];
      for (diff[1] = -1; diff[1] < 2; diff[1]++) {
        currentIndex[1] = index3d[1] + diff[1];
        for (diff[2] = -1; diff[2] < 2; diff[2]++) {
          currentIndex[2] = index3d[2] + diff[2];
          // check if there exists a cell with the specified coordinates
          bool isPossibleHaloCell = false;
          bool isValidCell = true;
          for (int i = 0; i < 3; i++) {
            isPossibleHaloCell |=
                currentIndex[i] == 0 ||
                currentIndex[i] == _cellsPerDimensionWithHalo[i] - 1;
            isValidCell &= currentIndex[i] < _cellsPerDimensionWithHalo[i] &&
                           currentIndex[i] >= 0;
          }
          if (isPossibleHaloCell && isValidCell) {
            std::array<std::array<double, 3>, 2> boxBound {};
            getCellBoundingBox(index3d, boxBound[0], boxBound[1]);
            bool close = true;
            for (int i = 0; i < 3; i++) {
              if (diff[i] < 0) {
                if (position[i] - boxBound[1][i] > allowedDistance) {
                  close = false;
                }
              } else if (diff[i] > 0) {
                if (boxBound[0][i] - position[i] > allowedDistance) {
                  close = false;
                }
              }
            }
            if (close) {
              closeHaloCells.push_back(&getCell(currentIndex));
            }
          }
        }
      }
    }

    return closeHaloCells;
  }

 private:
  std::array<index_t, 3> index3D(index_t index1d) const;
  index_t index1D(const std::array<index_t, 3> &index3d) const;

  std::array<index_t, 3> _cellsPerDimensionWithHalo;
  index_t _numCells;
  std::vector<ParticleCell> *_vec1D;

  std::array<double, 3> _boxMin, _boxMax;
  std::array<double, 3> _haloBoxMin, _haloBoxMax;

  double _interactionLength;
  //	int _cellsPerInteractionLength; = 1 hardcode to 1 for now, because flag
  // manager will also need to be adapted

  std::array<double, 3> _cellLength;

  // 1 over above. Since this value is needed for sorting particles in cells, it
  // is computed quite often
  std::array<double, 3> _cellLengthReciprocal;

  //	CellBorderAndFlagManager _cellBorderAndFlagManager;
};

template <class ParticleCell>
inline typename CellBlock3D<ParticleCell>::index_t
CellBlock3D<ParticleCell>::get1DIndexOfPosition(
    const std::array<double, 3> &pos) const {
  return index1D(get3DIndexOfPosition(pos));
}

template <class ParticleCell>
inline std::array<typename CellBlock3D<ParticleCell>::index_t, 3>
CellBlock3D<ParticleCell>::get3DIndexOfPosition(
    const std::array<double, 3> &pos) const {
  std::array<typename CellBlock3D<ParticleCell>::index_t, 3> cellIndex{};

  for (int dim = 0; dim < 3; dim++) {
    const long int value =
        (static_cast<long int>(
            floor((pos[dim] - _boxMin[dim]) * _cellLengthReciprocal[dim]))) +
        1l;
    const index_t nonnegativeValue = static_cast<index_t>(std::max(value, 0l));
    const index_t nonLargerValue =
        std::min(nonnegativeValue, _cellsPerDimensionWithHalo[dim] - 1);
    cellIndex[dim] = nonLargerValue;
    /// @todo this is a sanity check to prevent doubling of particles, but
    /// could be done better!
    if (pos[dim] >= _boxMax[dim]) {
      cellIndex[dim] = _cellsPerDimensionWithHalo[dim] - 1;
    } else if (pos[dim] < _boxMin[dim]) {
      cellIndex[dim] = 0;
    }
  }

  return cellIndex;
  // in very rare cases rounding is stupid, thus we need a check...
  /// @todo when the border and flag manager is there
}

template <class ParticleCell>
inline void CellBlock3D<ParticleCell>::rebuild(
    std::vector<ParticleCell> &vec, const std::array<double, 3> &bMin,
    const std::array<double, 3> &bMax, double interactionLength) {
  _vec1D = &vec;
  _boxMin = bMin;
  _boxMax = bMax;
  _interactionLength = interactionLength;

  // compute cell length
  _numCells = 1;
  for (int d = 0; d < 3; ++d) {
    double diff = _boxMax[d] - _boxMin[d];
    auto cellsPerDim =
        static_cast<index_t>(std::floor(diff / _interactionLength));
    // at least one central cell
    cellsPerDim = std::max(cellsPerDim, 1ul);

    _cellsPerDimensionWithHalo[d] = cellsPerDim + 2;

    _cellLength[d] = diff / cellsPerDim;
    _cellLengthReciprocal[d] =
        cellsPerDim / diff;  // compute with least rounding possible

    _haloBoxMin[d] = _boxMin[d] - _cellLength[d];
    _haloBoxMax[d] = _boxMax[d] + _cellLength[d];

    _numCells *= _cellsPerDimensionWithHalo[d];

    AutoPasLogger->debug("CellBlock3D: _cellsPerDimensionWithHalo[{}]={}", d,
                         _cellsPerDimensionWithHalo[d]);
  }

  _vec1D->resize(_numCells);

  //	ParticleCell::_cellBorderAndFlagManager.init(
  //		_cellsPerDimensionWithHalo,
  //		_haloBoxMin, _haloBoxMax,
  //		_boxMin, _boxMax,
  //		_cellLength);
}

template <class ParticleCell>
inline void CellBlock3D<ParticleCell>::getCellBoundingBox(
    const index_t index1d, std::array<double, 3> &boxmin,
    std::array<double, 3> &boxmax) {
  this->getCellBoundingBox(this->index3D(index1d), boxmin, boxmax);
}

template <class ParticleCell>
inline void CellBlock3D<ParticleCell>::getCellBoundingBox(
    const std::array<index_t, 3> &index3d, std::array<double, 3> &boxmin,
    std::array<double, 3> &boxmax) {
  for (int d = 0; d < 3; d++) {
    // defaults
    boxmin[d] = index3d[d] * this->_cellLength[d] + _haloBoxMin[d];
    boxmax[d] = (index3d[d] + 1) * this->_cellLength[d] + _haloBoxMin[d];

    // stupid rounding errors (makes sure that the lower corner is set
    // correctly!
    if (index3d[d] == 0) {
      boxmin[d] = _haloBoxMin[d];
      boxmax[d] = _boxMin[d];
    } else if (index3d[d] == 1) {
      boxmin[d] = _boxMin[d];
    }
    // no else!, as this might ALSO be 1
    if (index3d[d] == this->_cellsPerDimensionWithHalo[d] - 2) {
      boxmax[d] = _boxMax[d];
    } else if (index3d[d] == this->_cellsPerDimensionWithHalo[d] - 1) {
      boxmin[d] = _boxMax[d];
      boxmax[d] = _haloBoxMax[d];
    }
  }
}

template <class ParticleCell>
inline ParticleCell &CellBlock3D<ParticleCell>::getContainingCell(
    const std::array<double, 3> &pos) const {
  auto ind = get1DIndexOfPosition(pos);
  return getCell(ind);
}

template <class ParticleCell>
inline ParticleCell &CellBlock3D<ParticleCell>::getCell(index_t index1d) const {
  return _vec1D->at(index1d);
}

template <class ParticleCell>
inline ParticleCell &CellBlock3D<ParticleCell>::getCell(
    const std::array<index_t, 3> &index3d) const {
  return _vec1D->at(index1D(index3d));
}

template <class ParticleCell>
inline std::array<typename CellBlock3D<ParticleCell>::index_t, 3>
CellBlock3D<ParticleCell>::index3D(index_t index1d) const {
  return ThreeDimensionalMapping::oneToThreeD(index1d,
                                              _cellsPerDimensionWithHalo);
}

template <class ParticleCell>
inline typename CellBlock3D<ParticleCell>::index_t
CellBlock3D<ParticleCell>::index1D(
    const std::array<index_t, 3> &index3d) const {
  return ThreeDimensionalMapping::threeToOneD(index3d,
                                              _cellsPerDimensionWithHalo);
}

template <class ParticleCell>
bool CellBlock3D<ParticleCell>::checkInHalo(
    std::array<double, 3> position) const {
  return autopas::inBox(position, _haloBoxMin, _haloBoxMax) &&
         autopas::notInBox(position, _boxMin, _boxMax);
}
} /* namespace autopas */

#endif /* SRC_CONTAINERS_CELLBLOCK3D_H_ */
