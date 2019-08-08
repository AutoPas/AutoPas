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
template <class Particle, class ParticleCell>
class VerletClusterCellsParticleIterator : public ParticleIteratorInterfaceImpl<Particle> {
 public:
  /**
   * @param cont cells to iterate over
   * @param dummyStarts indices of the first dummy particle for each cell
   * @param offsetToDummy offset to add to create dummy particle
   * @param behavior iterator behavior
   */
  explicit VerletClusterCellsParticleIterator(std::vector<ParticleCell> *cont, std::vector<size_t> *dummyStarts,
                                              double offsetToDummy, IteratorBehavior behavior = haloAndOwned)
      : _vectorOfCells(cont), _dummyStarts(dummyStarts), cellId(0), _behavior(behavior), _offsetToDummy(offsetToDummy) {
    cellIter = (*_vectorOfCells)[0]._particles.begin();
    cellEnd = cellIter + (*_dummyStarts)[0];
    --cellIter;

    this->operator++();
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++()
   */
  inline VerletClusterCellsParticleIterator<Particle, ParticleCell> &operator++() override {
    do {
      ++cellIter;

      while (cellIter == cellEnd) {
        ++cellId;

        if (not isValid()) {
          return *this;
        }
        cellIter = (*_vectorOfCells)[cellId]._particles.begin();
        cellEnd = cellIter + (*_dummyStarts)[cellId];
      }
    } while (not fitsBehavior(*cellIter));
    return *this;
  }

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  Particle &operator*() const override { return *cellIter; }

  /**
   * Check whether the iterator currently points to a valid particle.
   * @return returns whether the iterator is valid
   */
  bool isValid() const override { return cellId < _vectorOfCells->size(); }

  /**
   * @copydoc ParticleIteratorInterface::deleteCurrentParticle()
   */
  void deleteCurrentParticle() override {
    auto pos = cellIter->getR();
    pos[0] += _offsetToDummy;
    cellIter->setR(pos);
    --(*_dummyStarts)[cellId];
    --cellEnd;
    std::iter_swap(cellIter, cellEnd);
    --cellIter;
  }

  ParticleIteratorInterfaceImpl<Particle> *clone() const override {
    return new VerletClusterCellsParticleIterator<Particle, ParticleCell>(*this);
  }

 protected:
  /**
   * Checks if particle fits the selected behavior
   * @param p particle
   * @return true if particle fits the behavior
   */
  bool fitsBehavior(Particle &p) {
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
   * Pointer to the cell vector
   */
  std::vector<ParticleCell> *_vectorOfCells;
  /**
   * Pointer to the vector with dummy indices
   */
  std::vector<size_t> *_dummyStarts;

  /**
   * Current cellId
   */
  size_t cellId;

  /**
   * Current iterator within a cell
   */
  typename std::vector<Particle>::iterator cellIter;
  /**
   * end for the current cell
   */
  typename std::vector<Particle>::iterator cellEnd;

  /**
   * The behavior of the iterator.
   */
  IteratorBehavior _behavior;

  /**
   * Offset to add to set particle aoutside domain
   */
  double _offsetToDummy;
};

/**
 * VerletClusterCellsRegionParticleIterator to iterate over all particles within a specific region
 * @tparam Particle Particle type over which the iterator iterates
 * @tparam ParticleCell Cell type over which the iterator iterates
 */
template <class Particle, class ParticleCell>
class VerletClusterCellsRegionParticleIterator : public VerletClusterCellsParticleIterator<Particle, ParticleCell> {
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
  explicit VerletClusterCellsRegionParticleIterator(std::vector<ParticleCell> *cont, std::vector<size_t> *dummyStarts,
                                                    std::array<double, 3> startRegion, std::array<double, 3> endRegion,
                                                    std::vector<size_t> &indicesInRegion, double offsetToDummy,
                                                    IteratorBehavior behavior = haloAndOwned, double skin = 0.0)
      : VerletClusterCellsParticleIterator<Particle, ParticleCell>(cont, dummyStarts, offsetToDummy, behavior),
        _startRegion(startRegion),
        _endRegion(endRegion),
        _indicesInRegion(indicesInRegion),
        _currentRegionIndex(0),
        _skin(skin) {
    _indicesInRegion.push_back(cont->size() + 1);
    this->cellId = _indicesInRegion[_currentRegionIndex];
    this->cellIter = (*this->_vectorOfCells)[this->cellId]._particles.begin();
    this->cellEnd = this->cellIter + (*this->_dummyStarts)[this->cellId];
    --this->cellIter;
    operator++();
  }

  /**
   * @copydoc ParticleIteratorInterface::operator++()
   */
  inline VerletClusterCellsParticleIterator<Particle, ParticleCell> &operator++() override {
    do {
      ++this->cellIter;

      while (this->cellIter == this->cellEnd) {
        this->cellId = _indicesInRegion[++_currentRegionIndex];
        if (not this->isValid()) {
          return *this;
        }
        this->cellIter = std::lower_bound(
            (*this->_vectorOfCells)[this->cellId]._particles.begin(),
            (*this->_vectorOfCells)[this->cellId]._particles.begin() + (*this->_dummyStarts)[this->cellId],
            _startRegion[2] - _skin, [](const Particle &a, const double b) { return a.getR()[2] < b; });
        this->cellEnd = std::upper_bound(
            (*this->_vectorOfCells)[this->cellId]._particles.begin(),
            (*this->_vectorOfCells)[this->cellId]._particles.begin() + (*this->_dummyStarts)[this->cellId],
            _endRegion[2] + _skin, [](const double b, const Particle &a) { return b < a.getR()[2]; });
      }
    } while ((not VerletClusterCellsParticleIterator<Particle, ParticleCell>::fitsBehavior(*this->cellIter)) or
             utils::notInBox(this->cellIter->getR(), _startRegion, _endRegion));
    return *this;
  }

  ParticleIteratorInterfaceImpl<Particle> *clone() const override {
    return new VerletClusterCellsRegionParticleIterator<Particle, ParticleCell>(*this);
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
   * indecies of possible cells
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
