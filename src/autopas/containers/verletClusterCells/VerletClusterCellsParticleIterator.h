/**
 * @file VerletClusterCellsParticleIterator.h
 * @author jspahl
 * @date 25.3.19
 */

#pragma once

#include <algorithm>
#include <utility>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/ParticleDeletedObserver.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/iterators/RegionParticleIterator.h"

namespace autopas::internal {
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
   * @param particleDeletedObserver Observer that should be called, when a particle is deleted. Can be nullptr.
   */
  explicit VerletClusterCellsParticleIterator(CellVecType *cont, DummyStartsType &dummyStarts, double offsetToDummy,
                                              IteratorBehavior behavior = haloAndOwned,
                                              ParticleDeletedObserver *particleDeletedObserver = nullptr)
      : _vectorOfCells(cont),
        _dummyStarts(dummyStarts),
        _cellId(0),
        _behavior(behavior),
        _offsetToDummy(offsetToDummy),
        _particleDeletedObserver(particleDeletedObserver) {
    // 1. set _cellId to thread number.
    _cellId = autopas_get_thread_num();

    if (_cellId >= _vectorOfCells->size()) {
      // prevent segfaults if the _cellId is too large
      return;
    }

    // 2. set cell iterators to appropriate start
    _iteratorWithinOneCell = (*_vectorOfCells)[_cellId]._particles.begin();
    _cellEnd = _iteratorWithinOneCell + getDummyStartbyIndex(_cellId);

    // 3. do a -- for _iteratorWithinOneCell to be able to call operator++ and still end up at the front of everything.
    --_iteratorWithinOneCell;

    // 4. call operator++
    VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable>::operator++();
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++()
   */
  inline VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable> &operator++() override {
    do {
      // increase the cell-local iterator
      ++_iteratorWithinOneCell;

      while (_iteratorWithinOneCell == _cellEnd) {
        // increase _cellId by stride if we are at the end of a cell!
        auto stride = autopas_get_num_threads();
        _cellId += stride;

        if (_cellId >= _vectorOfCells->size()) {
          // if we are at or beyond the end of all cells, return (the iterator is now invalid)
          return *this;
        }
        _iteratorWithinOneCell = (*_vectorOfCells)[_cellId]._particles.begin();
        _cellEnd = _iteratorWithinOneCell + getDummyStartbyIndex(_cellId);
      }
    } while (not fitsBehavior(*_iteratorWithinOneCell));
    return *this;
  }

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  ParticleType &operator*() const override { return *_iteratorWithinOneCell; }

  /**
   * Check whether the iterator currently points to a valid particle.
   * @return returns whether the iterator is valid
   */
  bool isValid() const override { return _cellId < _vectorOfCells->size(); }

  /**
   * @copydoc ParticleIteratorInterface::deleteCurrentParticleImpl()
   */
  void deleteCurrentParticleImpl() override {
    if constexpr (modifiable) {
      auto pos = _iteratorWithinOneCell->getR();
      pos[0] += _offsetToDummy;
      _iteratorWithinOneCell->setR(pos);
      --(_dummyStarts)[_cellId];
      --_cellEnd;
      std::iter_swap(_iteratorWithinOneCell, _cellEnd);
      --_iteratorWithinOneCell;
      if (_particleDeletedObserver) {
        _particleDeletedObserver->notifyParticleDeleted();
      }
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
      case IteratorBehavior::haloOwnedAndDummy:
        return true;
      case IteratorBehavior::haloAndOwned:
        return not p.isDummy();
      case IteratorBehavior::ownedOnly:
        return p.isOwned();
      case IteratorBehavior::haloOnly:
        return p.isHalo();
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
  CellIteratorType _iteratorWithinOneCell;
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

  /**
   * Observer that should be called, when a particle is deleted.
   * Can be nullptr.
   */
  ParticleDeletedObserver *_particleDeletedObserver{nullptr};
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
   * @param particleDeletedObserver Observer that should be called, when a particle is deleted. Can be nullptr.
   */
  explicit VerletClusterCellsRegionParticleIterator(CellVecType *cont, DummyStartsType &dummyStarts,
                                                    std::array<double, 3> startRegion, std::array<double, 3> endRegion,
                                                    std::vector<size_t> indicesInRegion, double offsetToDummy,
                                                    IteratorBehavior behavior = haloAndOwned, double skin = 0.0,
                                                    ParticleDeletedObserver *particleDeletedObserver = nullptr)
      : VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable>(cont, dummyStarts, offsetToDummy,
                                                                               behavior, particleDeletedObserver),
        _startRegion(startRegion),
        _endRegion(endRegion),
        _indicesInRegion(std::move(indicesInRegion)),
        _currentRegionIndex(0),
        _skin(skin) {
    // 1. add a last entry within _indicesInRegion that points beyond the end of _vectorOfCells.
    _indicesInRegion.push_back(this->_vectorOfCells->size() + 1);

    // 2. set _currentRegionIndex to current thread id
    // also ensure that _currentIndex does not go beyond _indicesInRegion.
    _currentRegionIndex = std::min(static_cast<unsigned long>(autopas_get_thread_num()), _indicesInRegion.size() - 1);

    // 3. set the _cellID to the appropriate value
    this->_cellId = _indicesInRegion[_currentRegionIndex];

    if (this->_cellId >= this->_vectorOfCells->size()) {
      // prevent segfaults if the _cellId is too large.
      return;
    }

    // 4. set _iteratorWithinOneCell to begin of that cell
    this->_iteratorWithinOneCell = (*this->_vectorOfCells)[this->_cellId]._particles.begin();

    // 5. adapt _cellEnd to dummies.
    this->_cellEnd = this->_iteratorWithinOneCell + this->getDummyStartbyIndex(this->_cellId);

    // 6. do a -- to be able to do a ++. The ++ is needed to ensure a valid iterator (if possible), as the initial
    // iterator might be invalid.
    --this->_iteratorWithinOneCell;

    // 7. do the ++
    VerletClusterCellsRegionParticleIterator<Particle, ParticleCell, modifiable>::operator++();
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++()
   */
  inline VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable> &operator++() override {
    do {
      // increase the cell-local iterator
      ++this->_iteratorWithinOneCell;

      while (this->_iteratorWithinOneCell == this->_cellEnd) {
        // Increase _currentRegionIndex by stride if we are at the end of a cell!
        // Also ensure that it does not go beyond _indicesInRegion.size() - 1.
        // The last entry of _indicesInRegion is invalid (>= _vectorOfCells->size()), so this is fine!
        _currentRegionIndex = std::min(_currentRegionIndex + autopas_get_num_threads(), _indicesInRegion.size() - 1);
        this->_cellId = _indicesInRegion[_currentRegionIndex];
        if (this->_cellId >= this->_vectorOfCells->size()) {
          return *this;
        }
        if (not modifiable) {
          // We can't delete particles, so optimizations assuming sorted particles are in order!

          // Find lowest possible particle in the region.
          // Sorting in the array is at most off by skin/2 as the container is rebuilt if a particle moves more.
          // increasing the search area by skin guarantees particles in the region to be found.
          this->_iteratorWithinOneCell = std::lower_bound(
              (*this->_vectorOfCells)[this->_cellId]._particles.begin(),
              (*this->_vectorOfCells)[this->_cellId]._particles.begin() + this->getDummyStartbyIndex(this->_cellId),
              _startRegion[2] - _skin, [](const Particle &a, const double b) { return a.getR()[2] < b; });
          this->_cellEnd = std::upper_bound(
              (*this->_vectorOfCells)[this->_cellId]._particles.begin(),
              (*this->_vectorOfCells)[this->_cellId]._particles.begin() + this->getDummyStartbyIndex(this->_cellId),
              _endRegion[2] + _skin, [](const double b, const Particle &a) { return b < a.getR()[2]; });
        } else {
          // If the iterator is modifiable, we can delete particles, which breaks the sorted assumption and thus
          // some particles might be skipped.
          this->_iteratorWithinOneCell = (*this->_vectorOfCells)[this->_cellId]._particles.begin();
          this->_cellEnd = this->_iteratorWithinOneCell + this->getDummyStartbyIndex(this->_cellId);
        }
      }
    } while ((not VerletClusterCellsParticleIterator<Particle, ParticleCell, modifiable>::fitsBehavior(
                 *this->_iteratorWithinOneCell)) or
             utils::notInBox(this->_iteratorWithinOneCell->getR(), _startRegion, _endRegion));
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
}  // namespace autopas::internal
