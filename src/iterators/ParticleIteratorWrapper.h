/**
 * ParticleIteratorWrapper.h
 *
 *  Created on: 29 May 2018
 *      Author: seckler
 */
#pragma once

#include "ParticleIteratorInterface.h"

namespace autopas {

/**
 * ParticleIteratorWrapper class is the main class visible to the user to iterate over particles of the AutoPas class.
 * The particles can be accessed using "iterator->" or "*iterator". The next particle using the ++operator, e.g.
 * "++iteratorWrapper".
 * @note The wrapper class provides an easy way to access the function of the ParticleIterators, e.g. allowing
 * "iteratorWrapper->getR()" or "++iteratorWrapper".
 * Without the wrapper class and using the ParticleIteratorInterface the calls would look like:
 * (*iteratorInterface)->getR() or "++(*iteratorWrapper)"
 * @tparam Particle type of the particle that is accessed
 */
template <class Particle>
class ParticleIteratorWrapper : public ParticleIteratorInterface<Particle> {
 public:
  /**
   * Constructor of the ParticleIteratorWrapper
   * @param particleIteratorInterface
   * @tparam InterfacePtrType type of the interface ptr
   */
  template <class InterfacePtrType>
  explicit ParticleIteratorWrapper(InterfacePtrType* particleIteratorInterface)
      : _particleIterator(static_cast<internal::ParticleIteratorInterfaceImpl<Particle>*>(particleIteratorInterface)) {}

  /**
   * copy operator
   * @param otherParticleIteratorWrapper the other ParticleIteratorWrapper
   */
  ParticleIteratorWrapper(const ParticleIteratorWrapper& otherParticleIteratorWrapper)
      : _particleIterator(otherParticleIteratorWrapper._particleIterator->clone()) {}

  /**
   * copy assignment operator to assign contents of otherParticleIteratorWrapper to this ParticleIteratorWrapper
   * @param otherParticleIteratorWrapper the other IteratorWrapper
   * @return *this
   */
  inline ParticleIteratorWrapper& operator=(const ParticleIteratorWrapper& otherParticleIteratorWrapper) {
    _particleIterator = std::unique_ptr<internal::ParticleIteratorInterfaceImpl<Particle>>(
        otherParticleIteratorWrapper._particleIterator->clone());
    return *this;
  }

  inline ParticleIteratorWrapper<Particle>& operator++() final {
    _particleIterator->operator++();
    return *this;
  }

  inline Particle& operator*() const final { return _particleIterator->operator*(); }

  inline void deleteCurrentParticle() final { _particleIterator->deleteCurrentParticle(); }

  inline bool isValid() const final { return _particleIterator->isValid(); }

 private:
  std::unique_ptr<autopas::internal::ParticleIteratorInterfaceImpl<Particle>> _particleIterator;
};

} /* namespace autopas */