/**
 * @file VerletClusterCellsParticleIterator.h
 * @author jspahl
 * @date 25.3.19
 */

#pragma once

#include <algorithm>
#include <vector>
#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"

namespace autopas {
namespace internal {
/**
 * ParticleIterator class to access particles inside a VerletClusterCells container efficiently.
 * The particles can be accessed using "iterator->" or "*iterator". The next
 * particle using the ++operator, e.g. "++iterator".
 * @tparam Particle Type of the particle that is accessed.
 * @tparam ParticleCell Type of the container
 */
template <class Particle, class ParticleCell, bool modifiable>
class VerletClusterCellsParticleIterator : public ParticleIteratorInterfaceImpl<Particle, modifiable> {
  using CellVecType = std::conditional_t<modifiable, std::vector<ParticleCell>, const std::vector<ParticleCell>>;
  using CellIteratorType = std::conditional_t<modifiable, typename std::vector<Particle>::iterator,
                                              typename std::vector<Particle>::const_iterator>;
  using DummyStartsType = std::conditional_t<modifiable, std::vector<size_t>, const std::vector<size_t>>;
  using ParticleType = std::conditional_t<modifiable, Particle, const Particle>;

 public:
  /**
   * Constructor
   * @param cont cells to iterate over
   * @param dummyStarts indices of the first dummy particle for each cell
   * @param offsetToDummy offset to add to create dummy particle
   * @param behavior iterator behavior
   */
  explicit VerletClusterCellsParticleIterator(CellVecType *cont, DummyStartsType &dummyStarts, double offsetToDummy,
                                              IteratorBehavior behavior = haloAndOwned)
      : _vectorOfCells(cont),
        _dummyStarts(dummyStarts),
        _cellId(0),
        _behavior(behavior),
        _offsetToDummy(offsetToDummy) {
    _cellIter = (*_vectorOfCells)[0]._particles.begin();
    _cellEnd = _cellIter + getDummyStartbyIndex(0);
    --_cellIter;

    VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable>::operator++();
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++()
   */
  inline VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable> &operator++() override {
    do {
      ++_cellIter;

      while (_cellIter == _cellEnd) {
        ++_cellId;

        if (not(_cellId < _vectorOfCells->size())) {
          return *this;
        }
        _cellIter = (*_vectorOfCells)[_cellId]._particles.begin();
        _cellEnd = _cellIter + getDummyStartbyIndex(_cellId);
      }
    } while (not fitsBehavior(*_cellIter));
    return *this;
  }

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  ParticleType &operator*() const override { return *_cellIter; }

  /**
   * Check whether the iterator currently points to a valid particle.
   * @return returns whether the iterator is valid
   */
  bool isValid() const override { return _cellId < _vectorOfCells->size(); }

  /**
   * @copydoc ParticleIteratorInterface::deleteCurrentParticle()
   */
  void deleteCurrentParticleImpl() override {
    if constexpr (modifiable) {
      auto pos = _cellIter->getR();
      pos[0] += _offsetToDummy;
      _cellIter->setR(pos);
      --(_dummyStarts)[_cellId];
      --_cellEnd;
      std::iter_swap(_cellIter, _cellEnd);
      --_cellIter;
    }
  }

  ParticleIteratorInterfaceImpl<Particle, modifiable> *clone() const override {
    return new VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable>(*this);
  }

 protected:
  /**
   * Checks if particle fits the selected behavior.
   * @param p particle
   * @return true if particle fits the behavior
   */
  bool fitsBehavior(const Particle &p) const {
    switch (_behavior) {
      case IteratorBehavior::haloAndOwned:
        return true;
      case IteratorBehavior::ownedOnly:
        return p.isOwned();
      case IteratorBehavior::haloOnly:
        return not p.isOwned();
    }
    return false;
  }

  /**
   * Returns dummy start for given cell.
   * @param index of cell
   * @return first index with dummy particle
   */
  size_t getDummyStartbyIndex(size_t index) { return _dummyStarts[index]; }
  /**
   * Pointer to the cell vector
   */
  CellVecType *_vectorOfCells;
  /**
   * Pointer to the vector with dummy indices
   */
  DummyStartsType &_dummyStarts;

  /**
   * Current cellId
   */
  size_t _cellId;

  /**
   * Current iterator within a cell
   */
  CellIteratorType _cellIter;
  /**
   * end for the current cell
   */
  CellIteratorType _cellEnd;

  /**
   * The behavior of the iterator.
   */
  IteratorBehavior _behavior;

  /**
   * Offset to add to set particle outside domain
   */
  double _offsetToDummy;
};

/**
 * VerletClusterCellsRegionParticleIterator to iterate over all particles within a specific region
 * @tparam Particle Particle type over which the iterator iterates
 * @tparam ParticleCell Cell type over which the iterator iterates
 */
template <class Particle, class ParticleCell, bool modifiable>
class VerletClusterCellsRegionParticleIterator
    : public VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable> {
  using CellVecType = std::conditional_t<modifiable, std::vector<ParticleCell>, const std::vector<ParticleCell>>;
  using DummyStartsType = std::conditional_t<modifiable, std::vector<size_t>, const std::vector<size_t>>;

 public:
  /**
   * Constructor of the VerletClusterCellsRegionParticleIterator.
   *
   * @param cont Container of particle cells.
   * @param dummyStarts indices of the first dummy particle for each cell
   * @param startRegion Lower corner of the region to iterate over.
   * @param endRegion Top corner of the region to iterate over.
   * @param indicesInRegion List of indices all threads will iterate over.
   * @param offsetToDummy offset to add to create dummy particle
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated through.
   * @param skin skin to which extend serch window
   */
  explicit VerletClusterCellsRegionParticleIterator(CellVecType *cont, DummyStartsType &dummyStarts,
                                                    std::array<double, 3> startRegion, std::array<double, 3> endRegion,
                                                    const std::vector<size_t> &indicesInRegion, double offsetToDummy,
                                                    IteratorBehavior behavior = haloAndOwned, double skin = 0.0)
      : VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable>(cont, dummyStarts, offsetToDummy,
                                                                               behavior),
        _startRegion(startRegion),
        _endRegion(endRegion),
        _indicesInRegion(indicesInRegion),
        _currentRegionIndex(0),
        _skin(skin) {
    _indicesInRegion.push_back(this->_vectorOfCells->size() + 1);
    this->_cellId = _indicesInRegion[_currentRegionIndex];
    this->_cellIter = (*this->_vectorOfCells)[this->_cellId]._particles.begin();
    this->_cellEnd = this->_cellIter + this->getDummyStartbyIndex(this->_cellId);
    --this->_cellIter;
    VerletClusterCellsRegionParticleIterator<Particle, ParticleCell, modifiable>::operator++();
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++()
   */
  inline VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable> &operator++() override {
    do {
      ++this->_cellIter;

      while (this->_cellIter == this->_cellEnd) {
        this->_cellId = _indicesInRegion[++_currentRegionIndex];
        if (not(this->_cellId < this->_vectorOfCells->size())) {
          return *this;
        }
        // Find lowest possible particle in the region.
        // Sorting in the array is at most off by skin/2 as the container is rebuilt if a particle moves more.
        // increasing the search area by skin guarantees particles in the region to be found.
        this->_cellIter = std::lower_bound(
            (*this->_vectorOfCells)[this->_cellId]._particles.begin(),
            (*this->_vectorOfCells)[this->_cellId]._particles.begin() + this->getDummyStartbyIndex(this->_cellId),
            _startRegion[2] - _skin, [](const Particle &a, const double b) { return a.getR()[2] < b; });
        this->_cellEnd = std::upper_bound(
            (*this->_vectorOfCells)[this->_cellId]._particles.begin(),
            (*this->_vectorOfCells)[this->_cellId]._particles.begin() + this->getDummyStartbyIndex(this->_cellId),
            _endRegion[2] + _skin, [](const double b, const Particle &a) { return b < a.getR()[2]; });
      }
    } while (
        (not VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable>::fitsBehavior(*this->_cellIter)) or
        utils::notInBox(this->_cellIter->getR(), _startRegion, _endRegion));
    return *this;
  }

  ParticleIteratorInterfaceImpl<Particle, modifiable> *clone() const override {
    return new VerletClusterCellsRegionParticleIterator<Particle, ParticleCell, modifiable>(*this);
  }

 protected:
  /**
   * lower corner of target region
   */
  const std::array<double, 3> _startRegion;
  
  /**
   * upper corner of target region
   */
  const std::array<double, 3> _endRegion;

  /**
   * indices of possible cells
   */
  std::vector<size_t> _indicesInRegion;

  /**
   * current index in _indicesInRegion
   */
  size_t _currentRegionIndex;

  /**
   * skin for search
   */
  const double _skin;
};
}  // namespace internal
}  // namespace autopas
