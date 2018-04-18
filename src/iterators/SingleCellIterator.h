/*
 * SingleCellIterator.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SINGLECELLITERATOR_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SINGLECELLITERATOR_H_

#include "cells/RMMParticleCell2T.h"

namespace autopas {

/**
 * SingleCellIterator class to loop over particles of a single cell.
 *
 * @tparam Particle type of the Particles
 * @tparam ParticleCell type of the ParticleCell
 */
template <class Particle, class ParticleCell>
class SingleCellIterator {
 public:
  /**
   * default constructor of SingleCellIterator
   */
  SingleCellIterator() : _cell(nullptr), _index(0), _deleted(false) {}

  /**
   * constructor of SingleCellIterator
   * @param cell_arg pointer to the cell of particles
   * @param ind index of the first particle
   */
  explicit SingleCellIterator(ParticleCell *cell_arg, int ind = 0)
      : _cell(cell_arg), _index(ind), _deleted(false) {}

  /**
   * destructor of SingleCellIterator
   */
  virtual ~SingleCellIterator() = default;

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  Particle &operator*() const {
    Particle *ptr = nullptr;
    //_cell->particleAt(_index, ptr);
    ptr = &(_cell->_particles.at(_index));
    return *ptr;
  }

  /**
   * access particle using "iterator->"
   *
   * this is the member of pointer operator
   * @return current particle
   */
  Particle *operator->() const { return &(this->operator*()); }

  /**
   * increment operator to get the next particle
   * @return the next particle, usually ignored
   */
  SingleCellIterator &operator++() {
    if (not _deleted) ++_index;
    _deleted = false;
  }

  /**
   * equality operator.
   * if both iterators are invalid or if they point to the same particle, this
   * returns true
   * @param rhs
   * @return
   */
  bool operator==(const SingleCellIterator &rhs) const {
    return (not rhs.isValid() and not this->isValid()) or
           (_cell == rhs._cell && _index == rhs._index);
  }

  /**
   * inequality operator
   * descrition see operator==
   * @param rhs
   * @return
   */
  bool operator!=(const SingleCellIterator &rhs) const {
    return !(rhs == *this);
  }
  /**
   * Check whether the iterator is valid
   * @return returns whether the iterator is valid
   */
  bool isValid() const {
    return _cell != nullptr and _index < _cell->numParticles();
  }

  /**
   * Get the index of the particle in the cell
   * @return index of the current particle
   */
  int getIndex() const { return _index; }

  /**
   * Deletes the current particle
   */
  void deleteCurrentParticle() {
    _cell->deleteByIndex(_index);
    _deleted = true;
  }

 private:
  ParticleCell *_cell;
  int _index;
  bool _deleted;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SINGLECELLITERATOR_H_ */
