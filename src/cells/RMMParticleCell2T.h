/*
 * RMMParticleCell2T.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef AUTOPAS_RMMPARTICLECELL_H_
#define AUTOPAS_RMMPARTICLECELL_H_

#include <utils/SoA.h>
#include "ParticleCell.h"

namespace autopas {

/**
 * Reduced Memory Mode ParticleCell
 * This cell type does not store particles explicitly. Instead, the particles
 * are stored directly in a structure of array.
 * @todo this currently does not store all information of the particles. Only
 * position and forces and thus does not work
 * @tparam Particle type of particle to be stored
 */
template <class Particle, class Iterator>
class RMMParticleCell2T
    : public ParticleCell<Particle, Iterator,
                          RMMParticleCell2T<Particle, Iterator>> {
 public:
  /**
   * Constructor of RMMParticleCell
   */
  RMMParticleCell2T() {
    _molsSoABuffer.initArrays(
        {Particle::AttributeNames::posX, Particle::AttributeNames::posY,
         Particle::AttributeNames::posZ, Particle::AttributeNames::forceX,
         Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ});
  }

  void addParticle(Particle &m) override {
    _molsSoABuffer.push(Particle::AttributeNames::posX, m.getR()[0]);
    _molsSoABuffer.push(Particle::AttributeNames::posY, m.getR()[1]);
    _molsSoABuffer.push(Particle::AttributeNames::posZ, m.getR()[2]);
    _molsSoABuffer.push(Particle::AttributeNames::forceX, m.getF()[0]);
    _molsSoABuffer.push(Particle::AttributeNames::forceY, m.getF()[1]);
    _molsSoABuffer.push(Particle::AttributeNames::forceZ, m.getF()[2]);
  }

  unsigned long numParticles() const override {
    return _molsSoABuffer.getNumParticles();
  }
  bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override { _molsSoABuffer.clear(); }

  void deleteByIndex(int index) override {
    assert(index >= 0 and index < numParticles());
    assert(numParticles() > 0);
    if (index < numParticles() - 1) {
      _molsSoABuffer.swap(index, numParticles() - 1);
    }
    _molsSoABuffer.pop_back();
  }

  /**
   * the soa buffer of the particle, all information is stored here.
   */
  SoA _molsSoABuffer;

  /**
   * iterator to iterate through ParticleCell
   * If you need to explicitly store this iterator use
   * typename RMMParticleCell<ParticleType>::iterator iter;
   */
  typedef Iterator iterator;

 private:
  void buildMoleculeFromSoA(unsigned int i, Particle *&rmm_or_not_pointer) {
    rmm_or_not_pointer->setR(_molsSoABuffer.read<3>(
        {Particle::AttributeNames::posX, Particle::AttributeNames::posY,
         Particle::AttributeNames::posZ},
        i));
    rmm_or_not_pointer->setF(_molsSoABuffer.read<3>(
        {Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
         Particle::AttributeNames::forceZ},
        i));
  }

  template <class ParticleType>
  friend class RMMParticleCellIterator;
};

/**
 * specialization of the  class for the RMMParticleCell
 * @tparam Particle
 */
template <class Particle>
class RMMParticleCellIterator {
 public:
  /**
   * default constructor of SingleCellIterator
   */
  RMMParticleCellIterator() : _cell(nullptr), _index(0), _deleted(false) {}

  /**
   * constructor of SingleCellIterator
   * @param cell_arg pointer to the cell of particles
   * @param ind index of the first particle
   */
  RMMParticleCellIterator(
      RMMParticleCell2T<Particle, RMMParticleCellIterator<Particle>> *cell_arg,
      int ind = 0)
      : _cell(cell_arg), _index(ind), _deleted(false) {}

  //  SingleCellIterator(const SingleCellIterator &cellIterator) {
  //    _cell = cellIterator._cell;
  //    _AoSReservoir = cellIterator._AoSReservoir;
  //    _index = cellIterator._index;
  //  }

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  Particle &operator*() {
    // Particle * ptr = nullptr;
    // ptr = const_cast<Particle *>(& _AoSReservoir);
    Particle *ptr = &_AoSReservoir;
    _cell->buildMoleculeFromSoA(_index, ptr);
    //_cell->particleAt(_index, ptr);
    return *ptr;
  }

  /**
   * access particle using "iterator->"
   *
   * this is the member of pointer operator
   * @return current particle
   */
  Particle *operator->() { return &(this->operator*()); }

  /**
   * increment operator to get the next particle
   * @return the next particle, usually ignored
   */
  RMMParticleCellIterator &operator++() {
    if (not _deleted) ++_index;
    _deleted = false;
    return *this;
  }

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

  /**
   * Deletes the current particle
   */
  void deleteCurrentParticle() {
    _cell->deleteByIndex(_index);
    _deleted = true;
  }

 private:
  RMMParticleCell2T<Particle, RMMParticleCellIterator<Particle>> *_cell;
  Particle _AoSReservoir;
  int _index;
  bool _deleted;
};

// provide a simpler template for RMMParticleCell, i.e.
// RMMParticleCell<Particle>
template <class Particle>
using RMMParticleCell =
    RMMParticleCell2T<Particle, RMMParticleCellIterator<Particle>>;

} /* namespace autopas */

#endif /* AUTOPAS_RMMPARTICLECELL_H_ */
