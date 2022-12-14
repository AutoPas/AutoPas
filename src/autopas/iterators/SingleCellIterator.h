/**
 * @file SingleCellIterator.h
 *
 * @date 31 May 2018
 * @author tchipevn
 */

#pragma once

#include "autopas/iterators/SingleCellIteratorInterface.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas::internal {

/**
 * SingleCellIterator class to loop over particles of a single cell.
 *
 * @tparam Particle type of the Particles
 * @tparam ParticleCell type of the ParticleCell.
 * @tparam modifiable Defines whether the ParticleIterator is modifiable or not. If it is false, it points to a const
 * Particle.
 */
template <class Particle, class ParticleCell, bool modifiable>
class SingleCellIterator final : public SingleCellIteratorInterfaceImpl<Particle, modifiable> {
  using CellType = std::conditional_t<modifiable, ParticleCell, const ParticleCell>;
  using ParticleType = std::conditional_t<modifiable, Particle, const Particle>;

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
  explicit SingleCellIterator(CellType *cell_arg, size_t ind = 0) : _cell(cell_arg), _index(ind), _deleted(false) {}

  /**
   * destructor of SingleCellIterator
   */
  virtual ~SingleCellIterator() = default;

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  inline ParticleType &operator*() const override { return (*_cell)[_index]; }

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
  inline bool operator==(const SingleCellIteratorInterface<Particle, modifiable> &rhs) const override {
    if (auto other = dynamic_cast<const SingleCellIterator<Particle, ParticleCell, modifiable> *>(&rhs)) {
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
  inline bool operator!=(const SingleCellIteratorInterface<Particle, modifiable> &rhs) const override {
    return !(rhs == *this);
  }
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

  inline SingleCellIteratorInterfaceImpl<Particle, modifiable> *clone() const override {
    return new SingleCellIterator<Particle, ParticleCell, modifiable>(*this);
  }

 protected:
  inline void deleteCurrentParticleImpl() override {
    if constexpr (modifiable) {
      _cell->deleteByIndex(_index);
      _deleted = true;
    } else {
      utils::ExceptionHandler::exception("Error: Trying to delete a particle through a const iterator.");
    }
  }

 private:
  CellType *_cell;
  size_t _index;
  bool _deleted;
};

}  // namespace autopas::internal
