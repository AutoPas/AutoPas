/*
 * ParticleIterator.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLEITERATOR_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLEITERATOR_H_

#include <vector>
#include "SingleCellIterator.h"

namespace autopas {

/**
 * ParticleIterator class to access particles inside of a container.
 * The particles can be accessed using "iterator->" or "*iterator". The next
 * particle using the ++operator, e.g. "++iterator"
 * @tparam Particle type of the particle that is accessed
 * @tparam ParticleCell type of the cell of the underlying data structure of the
 * container
 */
template <class Particle, class ParticleCell>
class ParticleIterator {
 public:
  /**
   * Constructor of the ParticleIterator class.
   *
   * @param cont linear data vector of ParticleCells
   */
  explicit ParticleIterator(std::vector<ParticleCell>* cont)
      : _vectorOfCells(cont),
        _iteratorAcrossCells(cont->begin()),
        _iteratorWithinOneCell() {
    if (_iteratorAcrossCells < cont->end()) {
      _iteratorWithinOneCell =
          (typename ParticleCell::iterator)(&(*_iteratorAcrossCells));
      if (not _iteratorWithinOneCell.isValid()) {
        next_non_empty_cell();
      }
    }
  }

  /**
   * Increment operator.
   * Used to jump to the next particle
   * @return next particle, usually ignored
   */
  inline ParticleIterator<Particle, ParticleCell>& operator++() {
    if (_iteratorWithinOneCell.isValid()) {
      ++_iteratorWithinOneCell;
    }

    // don't merge into if-else, _cell_iterator may becoeme invalid after ++

    if (not _iteratorWithinOneCell.isValid()) {
      next_non_empty_cell();
    }
    return *this;
  }

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  Particle& operator*() { return _iteratorWithinOneCell.operator*(); }

  /**
   * access particle using iterator->
   *
   * this is the member of pointer operator
   * @return current particle
   */
  Particle* operator->() { return &(this->operator*()); }  //

  /**
   * Deletes the current particle
   * @todo implement deletion
   */
  void deleteCurrentParticle() {
    //		cout << "deleteCurrentParticle is still ToDo" << endl;
  }

  /**
   * Check whether the iterator is valid
   * @return returns whether the iterator is valid
   */
  bool isValid() {
    return _vectorOfCells != nullptr and
           _iteratorAcrossCells < _vectorOfCells->end() and
           _iteratorWithinOneCell.isValid();
  }

 protected:
  /**
   * the next non-empty cell is selected
   */
  void next_non_empty_cell() {
    // find the next non-empty cell
    const int stride = 1;  // num threads
    for (_iteratorAcrossCells += stride;
         _iteratorAcrossCells < _vectorOfCells->end();
         _iteratorAcrossCells += stride) {
      const ParticleCell& c = *_iteratorAcrossCells;

      if (c.isNotEmpty()) {
        _iteratorWithinOneCell = (typename ParticleCell::iterator)(
            &(*_iteratorAcrossCells));
        break;
      }
    }
  }

 private:
  std::vector<ParticleCell>* _vectorOfCells;
  typename std::vector<ParticleCell>::iterator _iteratorAcrossCells;
  typename ParticleCell::iterator _iteratorWithinOneCell;
  //SingleCellIterator<Particle, ParticleCell> _iteratorWithinOneCell;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLEITERATOR_H_ */
