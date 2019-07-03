/**
 * @file SingleCellIterator.h
 *
 * @date 31 May 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/iterators/SingleCellIteratorInterface.h"

namespace autopas::internal {
/**
 * SingleCellIterator class to loop over particles of a single cell.
 *
 * @tparam Particle type of the Particles
 * @tparam ParticleCell type of the ParticleCell
 */
template <class Particle, class ParticleCell>
class SingleCellIterator : public SingleCellIteratorInterfaceImpl<Particle> {
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
  explicit SingleCellIterator(ParticleCell *cell_arg, size_t ind = 0) : _cell(cell_arg), _index(ind), _deleted(false) {}

  /**
   * destructor of SingleCellIterator
   */
  virtual ~SingleCellIterator() = default;

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  inline Particle &operator*() const override {
    Particle *ptr = nullptr;
    //_cell->particleAt(_index, ptr);
    ptr = &(_cell->_particles.at(_index));
    return *ptr;
  }

  /**
   * increment operator to get the next particle
   * @return the next particle, usually ignored
   */
  inline SingleCellIterator &operator++() override {
    if (not _deleted) ++_index;
    _deleted = false;
    return *this;
  }

  /**
   * equality operator.
   * if both iterators are invalid or if they point to the same particle, this
   * returns true
   * @param rhs
   * @return true if the iterators point to the same particle (in the same
   * cell), false otherwise
   */
  inline bool operator==(const SingleCellIteratorInterface<Particle> &rhs) const override {
    if (auto other = dynamic_cast<const SingleCellIterator<Particle, ParticleCell> *>(&rhs)) {
      return (not rhs.isValid() and not this->isValid()) or (_cell == other->_cell && _index == other->_index);
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
  inline bool operator!=(const SingleCellIteratorInterface<Particle> &rhs) const override { return !(rhs == *this); }
  /**
   * Check whether the iterator is valid
   * @return returns whether the iterator is valid
   */
  inline bool isValid() const override { return _cell != nullptr and _index < _cell->numParticles(); }

  /**
   * Get the index of the particle in the cell
   * @return index of the current particle
   */
  inline size_t getIndex() const override { return _index; }

  /**
   * Deletes the current particle
   */
  inline void deleteCurrentParticle() override {
    _cell->deleteByIndex(_index);
    _deleted = true;
  }

  inline SingleCellIteratorInterfaceImpl<Particle> *clone() const override {
    return new SingleCellIterator<Particle, ParticleCell>(*this);
  }

 private:
  ParticleCell *_cell;
  size_t _index;
  bool _deleted;
};

}  // namespace autopas::internal