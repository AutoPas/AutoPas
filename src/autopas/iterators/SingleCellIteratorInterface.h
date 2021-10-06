/**
 * @file SingleCellIteratorInterface.h
 *
 * @date 31 May 2018
 * @author seckler
 */

#pragma once

#include <cstddef>

#include "autopas/iterators/ParticleIteratorInterface.h"

namespace autopas {

/**
 * SingleCellIteratorInterface class to loop over particles of a single cell.
 *
 * @tparam Particle type of the Particles
 * @tparam modifiable Defines whether the ParticleIterator is modifiable or not. If it is false, it points to a const
 * Particle.
 */
template <class Particle, bool modifiable>
class SingleCellIteratorInterface : public ParticleIteratorInterface<Particle, modifiable> {
 public:
  /**
   * default constructor of SingleCellIteratorInterface
   */
  SingleCellIteratorInterface() = default;

  /**
   * destructor of SingleCellIteratorInterface
   */
  virtual ~SingleCellIteratorInterface() = default;

  /**
   * equality operator.
   * if both iterators are invalid or if they point to the same particle, this
   * returns true
   * @param rhs
   * @return
   */
  virtual bool operator==(const SingleCellIteratorInterface<Particle, modifiable> &rhs) const = 0;

  /**
   * inequality operator
   * descrition see operator==
   * @param rhs
   * @return
   */
  virtual bool operator!=(const SingleCellIteratorInterface<Particle, modifiable> &rhs) const = 0;

  /**
   * Get the index of the particle in the cell
   * @return index of the current particle
   */
  virtual size_t getIndex() const = 0;
};

namespace internal {
/**
 * All implementations of the interface should inherit from this class. It extends the interface just by the clone
 * method, which is needed by the Wrapper.
 * @tparam Particle
 * @tparam modifiable
 */
template <class Particle, bool modifiable>
class SingleCellIteratorInterfaceImpl : public SingleCellIteratorInterface<Particle, modifiable> {
 public:
  /**
   * Clones the current object, should allocate new object and return it.
   * @return the clone
   */
  virtual SingleCellIteratorInterfaceImpl *clone() const = 0;
};
}  // namespace internal
}  // namespace autopas

/**
 * gets a static cell iterator from an iteratorwrapper of a cell.
 * @tparam cell_t
 * @param cell the cell for which the iterator should be get
 * @return static cell iterator
 */
template <typename cell_t>
typename cell_t::iterator_t getStaticCellIter(cell_t &cell) {
  auto _wrapper = cell.begin();
  auto _ptr = _wrapper.get();

  return static_cast<typename cell_t::iterator_t &>(*_ptr);
}
