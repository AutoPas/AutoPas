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
 * ParticleIterator class to access particles inside a VerletClusterCells container efficently.
 * The particles can be accessed using "iterator->" or "*iterator". The next
 * particle using the ++operator, e.g. "++iterator".
 * @tparam Particle Type of the particle that is accessed.
 * @tparam ParticleCell Type of the container
 */
template <class Particle, class ParticleCell>
class VerletClusterCellsParticleIterator : public ParticleIteratorInterfaceImpl<Particle> {
 public:
  explicit VerletClusterCellsParticleIterator(std::vector<ParticleCell> *cont, std::vector<size_t> *dummyStarts,
                                              IteratorBehavior behavior = haloAndOwned)
      : _vectorOfCells(cont), _dummyStarts(dummyStarts), cellId(0), _behavior(behavior) {
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
    auto iterCopy = cellIter;
    --cellIter;
    (*_vectorOfCells)[cellId]._particles.erase(iterCopy);
  }

  ParticleIteratorInterfaceImpl<Particle> *clone() const override {
    return new VerletClusterCellsParticleIterator<Particle, ParticleCell>(*this);
  }

 protected:
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

  std::vector<ParticleCell> *_vectorOfCells;
  std::vector<size_t> *_dummyStarts;

  size_t cellId;
  typename std::vector<Particle>::iterator cellIter;
  typename std::vector<Particle>::iterator cellEnd;

  /**
   * The behavior of the iterator.
   */
  IteratorBehavior _behavior;
};

/**
 * RegionParticleIterator to iterate over all particles within a specific region
 * @todo optimize the region particle iterater. Currently we iterate over all
 * particles
 * @tparam Particle Particle type over which the iterator iterates
 * @tparam ParticleCell Cell type over which the iterator iterates
 */
template <class Particle, class ParticleCell>
class VerletClusterCellsRegionParticleIterator : public VerletClusterCellsParticleIterator<Particle, ParticleCell> {
 public:
  /**
   * Constructor of the RegionParticleIterator.
   *
   * @param cont Container of particle cells.
   * @param startRegion Lower corner of the region to iterate over.
   * @param endRegion Top corner of the region to iterate over.
   * @param indicesInRegion List of indices all threads will iterate over.
   * @param flagManager The CellBorderAndFlagManager that shall be used to query the cell types.
   * Can be nullptr if the behavior is haloAndOwned.
   * @param behavior The IteratorBehavior that specifies which type of cells shall be iterated through.
   */
  explicit VerletClusterCellsRegionParticleIterator(std::vector<ParticleCell> *cont, std::vector<size_t> *dummyStarts,
                                                    std::array<double, 3> startRegion, std::array<double, 3> endRegion,
                                                    std::vector<size_t> &indicesInRegion,
                                                    IteratorBehavior behavior = haloAndOwned)
      : VerletClusterCellsParticleIterator<Particle, ParticleCell>(cont, dummyStarts, behavior),
        _startRegion(startRegion),
        _endRegion(endRegion),
        _indicesInRegion(indicesInRegion),
        _currentRegionIndex(0) {
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
            _startRegion[2], [](const Particle &a, const double b) { return a.getR()[2] <= b; });
        this->cellEnd = std::upper_bound(
            (*this->_vectorOfCells)[this->cellId]._particles.begin(),
            (*this->_vectorOfCells)[this->cellId]._particles.begin() + (*this->_dummyStarts)[this->cellId],
            _endRegion[2], [](const double b, const Particle &a) { return b <= a.getR()[2]; });
      }
    } while ((not VerletClusterCellsParticleIterator<Particle, ParticleCell>::fitsBehavior(*this->cellIter)) or
             utils::notInBox(this->cellIter->getR(), _startRegion, _endRegion));
    return *this;
  }

  ParticleIteratorInterfaceImpl<Particle> *clone() const override {
    return new VerletClusterCellsRegionParticleIterator<Particle, ParticleCell>(*this);
  }

 protected:
  const std::array<double, 3> _startRegion;
  const std::array<double, 3> _endRegion;
  std::vector<size_t> _indicesInRegion;
  size_t _currentRegionIndex;
};
}  // namespace internal
}  // namespace autopas
