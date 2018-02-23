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

template <class Particle, class ParticleCell>
class SingleCellIterator {
 public:
  SingleCellIterator() : _cell(nullptr), _index(0) {}
  explicit SingleCellIterator(ParticleCell *cell_arg, int ind = 0)
      : _cell(cell_arg), _index(ind) {}
  Particle &operator*() const {
    Particle *ptr = nullptr;
    _cell->moleculesAt(_index, ptr);
    return *ptr;
  }

  virtual ~SingleCellIterator() = default;

  Particle *operator->() const { return &(this->operator*()); }
  SingleCellIterator &operator++() { ++_index; }
  bool isValid() { return _cell != nullptr and _index < _cell->numParticles(); }

  int getIndex() const { return _index; }

 private:
  ParticleCell *_cell;
  int _index;
};

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
