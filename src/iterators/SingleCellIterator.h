/*
 * SingleCellIterator.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SINGLECELLITERATOR_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SINGLECELLITERATOR_H_

#include "cells/RMMParticleCell.h"

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
  SingleCellIterator() : _cell(nullptr), _index(0) {}

  /**
   * constructor of SingleCellIterator
   * @param cell_arg pointer to the cell of particles
   * @param ind index of the first particle
   */
  explicit SingleCellIterator(ParticleCell *cell_arg, int ind = 0)
      : _cell(cell_arg), _index(ind) {}
  Particle &operator*() const {
    Particle *ptr = nullptr;
    _cell->moleculesAt(_index, ptr);
    return *ptr;
  }

  /**
   * destructor of SingleCellIterator
   */
  virtual ~SingleCellIterator() = default;

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
  SingleCellIterator &operator++() { ++_index; }

  /**
   * Check whether the iterator is valid
   * @return returns whether the iterator is valid
   */
  bool isValid() { return _cell != nullptr and _index < _cell->numParticles(); }

  /**
   * Get the index of the particle in the cell
   * @return index of the current particle
   */
  int getIndex() const { return _index; }

 private:
  ParticleCell *_cell;
  int _index;
};

/**
 * specialization of the SingleCellIterator class for the RMMParticleCell
 * @tparam Particle
 */
template <class Particle>
class SingleCellIterator<Particle, autopas::RMMParticleCell<Particle>> {
 public:
  SingleCellIterator() : _cell(nullptr), _index(0) {}
  explicit SingleCellIterator(RMMParticleCell<Particle> *cell_arg, int ind = 0)
      : _cell(cell_arg), _index(ind) {}
  SingleCellIterator(const SingleCellIterator &cellIterator) {
    _cell = cellIterator._cell;
    _AoSReservoir = cellIterator._AoSReservoir;
    _index = cellIterator._index;
  }
  Particle &operator*() {
    // Particle * ptr = nullptr;
    // ptr = const_cast<Particle *>(& _AoSReservoir);
    Particle *ptr = &_AoSReservoir;
    _cell->moleculesAt(_index, ptr);
    return *ptr;
  }
  Particle *operator->() { return &(this->operator*()); }
  SingleCellIterator &operator++() {
    ++_index;
    return *this;
  }
  bool isValid() { return _cell != nullptr and _index < _cell->numParticles(); }

  int getIndex() const { return _index; }

 private:
  RMMParticleCell<Particle> *_cell;
  Particle _AoSReservoir;
  int _index;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_SINGLECELLITERATOR_H_ */
