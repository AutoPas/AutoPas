/**
 * @file SingleCellIteratorInterface.h
 *
 * @date 31 May 2018
 * @author seckler
 */

#pragma once
#include "ParticleIteratorInterface.h"

namespace autopas {

/**
 * SingleCellIteratorInterface class to loop over particles of a single cell.
 *
 * @tparam Particle type of the Particles
 */
template <class Particle>
class SingleCellIteratorInterface : public ParticleIteratorInterface<Particle> {
 public:
  /**
   * default constructor of SingleCellIteratorInterface
   */
  SingleCellIteratorInterface() {}

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
  virtual bool operator==(const SingleCellIteratorInterface<Particle> &rhs) const = 0;

  /**
   * inequality operator
   * descrition see operator==
   * @param rhs
   * @return
   */
  virtual bool operator!=(const SingleCellIteratorInterface<Particle> &rhs) const = 0;

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
 */
template <class Particle>
class SingleCellIteratorInterfaceImpl : public SingleCellIteratorInterface<Particle> {
 public:
  /**
   * Clones the current object, should allocate new object and return it.
   * @return the clone
   */
  virtual SingleCellIteratorInterfaceImpl<Particle> *clone() const = 0;
};
}  // namespace internal
}  // namespace autopas
