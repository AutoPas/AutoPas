/**
 * @file RMMParticleCell2T.h
 * @date 17.01.2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include "autopas/cells/ParticleCell.h"
#include "autopas/iterators/ParticleIteratorInterface.h"
#include "autopas/utils/SoA.h"

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
class RMMParticleCell2T : public ParticleCell<Particle> {
 public:
  /**
   * Constructor of RMMParticleCell
   */
  RMMParticleCell2T() = default;

  /**
   * @copydoc ParticleCell::addParticle(Particle&)
   */
  void addParticle(Particle &p) override {
    _particleSoABuffer.template push<Particle::AttributeNames::id>(p.getID());
    _particleSoABuffer.template push<Particle::AttributeNames::posX>(p.getR()[0]);
    _particleSoABuffer.template push<Particle::AttributeNames::posY>(p.getR()[1]);
    _particleSoABuffer.template push<Particle::AttributeNames::posZ>(p.getR()[2]);
    _particleSoABuffer.template push<Particle::AttributeNames::forceX>(p.getF()[0]);
    _particleSoABuffer.template push<Particle::AttributeNames::forceY>(p.getF()[1]);
    _particleSoABuffer.template push<Particle::AttributeNames::forceZ>(p.getF()[2]);
  }

  SingleCellIteratorWrapper<Particle> begin() override {
    return SingleCellIteratorWrapper<Particle>(new Iterator(this));
  }

  unsigned long numParticles() const override { return _particleSoABuffer.getNumParticles(); }
  bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override { _particleSoABuffer.clear(); }

  void deleteByIndex(size_t index) override {
    assert(index >= 0 and index < numParticles());
    assert(numParticles() > 0);
    if (index < numParticles() - 1) {
      _particleSoABuffer.swap(index, numParticles() - 1);
    }
    _particleSoABuffer.pop_back();
  }

  /**
   * the soa buffer of the particle, all information is stored here.
   */
  SoA<typename Particle::SoAArraysType> _particleSoABuffer;

 private:
  void buildParticleFromSoA(size_t i, Particle *&rmm_or_not_pointer) {
    rmm_or_not_pointer->setR(
        _particleSoABuffer.template readMultiple<Particle::AttributeNames::posX, Particle::AttributeNames::posY,
                                                 Particle::AttributeNames::posZ>(i));
    rmm_or_not_pointer->setF(
        _particleSoABuffer.template readMultiple<Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
                                                 Particle::AttributeNames::forceZ>(i));
  }

  void writeParticleToSoA(size_t index, Particle &particle) {
    _particleSoABuffer.template writeMultiple<Particle::AttributeNames::posX, Particle::AttributeNames::posY,
                                              Particle::AttributeNames::posZ>(index, particle.getR());
    _particleSoABuffer.template writeMultiple<Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
                                              Particle::AttributeNames::forceZ>(index, particle.getF());
  }

  /**
   * iterator friend class
   * @tparam ParticleType
   */
  template <class ParticleType>
  friend class RMMParticleCellIterator;

  /**
   * type of the internal iterator
   */
  typedef Iterator iterator_t;
};

/**
 * SingleCellIterator for the RMMParticleCell
 * @tparam Particle
 */
template <class Particle>
class RMMParticleCellIterator : public internal::SingleCellIteratorInterfaceImpl<Particle> {
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
  explicit RMMParticleCellIterator(RMMParticleCell2T<Particle, RMMParticleCellIterator<Particle>> *cell_arg,
                                   size_t ind = 0)
      : _cell(cell_arg), _index(ind), _deleted(false) {}

  //  SingleCellIterator(const SingleCellIterator &cellIterator) {
  //    _cell = cellIterator._cell;
  //    _AoSReservoir = cellIterator._AoSReservoir;
  //    _index = cellIterator._index;
  //  }

  /**
   * @copydoc ParticleIteratorInterface::operator*()
   */
  inline Particle &operator*() const override {
    // Particle * ptr = nullptr;
    // ptr = const_cast<Particle *>(& _AoSReservoir);
    Particle *ptr = &_AoSReservoir;
    _cell->buildParticleFromSoA(_index, ptr);
    //_cell->particleAt(_index, ptr);
    return *ptr;
  }

  /**
   * equality operator.
   * if both iterators are invalid or if they point to the same particle, this
   * returns true
   * @param rhs
   * @return
   */
  bool operator==(const SingleCellIteratorInterface<Particle> &rhs) const override {
    if (auto other = dynamic_cast<const RMMParticleCellIterator *>(&rhs)) {
      return (not this->isValid() and not rhs.isValid()) or (_cell == other->_cell && _index == other->_index);
    } else {
      return false;
    }
  }

  /**
   * inequality operator
   * descrition see operator==
   * @param rhs
   * @return
   */
  bool operator!=(const SingleCellIteratorInterface<Particle> &rhs) const override { return !(rhs == *this); }

  /**
   * increment operator to get the next particle
   * @return the next particle, usually ignored
   */
  inline RMMParticleCellIterator &operator++() override {
    if (not _deleted) {
      _cell->writeParticleToSoA(_index, _AoSReservoir);
      ++_index;
    }
    _deleted = false;
    return *this;
  }

  /**
   * Check whether the iterator is valid
   * @return returns whether the iterator is valid
   */
  bool isValid() const override { return _cell != nullptr and _index < _cell->numParticles(); }

  /**
   * Get the index of the particle in the cell
   * @return index of the current particle
   */
  size_t getIndex() const override { return _index; }

  /**
   * Deletes the current particle
   */
  void deleteCurrentParticle() override {
    _cell->deleteByIndex(_index);
    _deleted = true;
  }

  internal::SingleCellIteratorInterfaceImpl<Particle> *clone() const override {
    return new RMMParticleCellIterator<Particle>(*this);
  }

 private:
  RMMParticleCell2T<Particle, RMMParticleCellIterator<Particle>> *_cell;
  mutable Particle _AoSReservoir;
  size_t _index;
  bool _deleted;
};

// provide a simpler template for RMMParticleCell, i.e.
// RMMParticleCell<Particle>
/**
 * typedef for simpler access to RMMParticleCell
 */
template <class Particle>
using RMMParticleCell = RMMParticleCell2T<Particle, RMMParticleCellIterator<Particle>>;

}  // namespace autopas
