/**
 * ParticleIterator.h
 *
 *  Created on: 17 Jan 2018
 *      Author: seckler
 */
#pragma once

namespace autopas {

/**
 * ParticleIteratorInterface class to iterate over particles.
 * This class provides a basic interface for all iterators within AutoPas.
 * The particles can be accessed using "iterator->" or "*iterator". The next
 * particle using the ++operator, e.g. "++iterator"
 * @tparam Particle type of the particle that is accessed
 */
template <class Particle>
class ParticleIteratorInterface {
 public:
  /**
   * Increment operator.
   * Used to jump to the next particle
   * @return next particle, usually ignored
   */
  virtual inline ParticleIteratorInterface<Particle>& operator++() = 0;

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  virtual Particle& operator*() = 0;

  /**
   * access particle using iterator->
   *
   * this is the member of pointer operator
   * @return current particle
   */
  virtual inline Particle* operator->() final { return &(this->operator*()); }

  /**
   * Deletes the current particle
   */
  virtual void deleteCurrentParticle() = 0;

  /**
   * Check whether the iterator is valid
   * @return returns whether the iterator is valid
   */
  virtual bool isValid() const = 0;
};

} /* namespace autopas */