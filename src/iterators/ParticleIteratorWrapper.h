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
 * ParticleIteratorInterface class to iterate over particles.
 * This class provides a basic interface for all iterators within AutoPas.
 * The particles can be accessed using "iterator->" or "*iterator". The next
 * particle using the ++operator, e.g. "++iterator"
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
  ParticleIteratorWrapper(InterfacePtrType* particleIteratorInterface)
      : _particleIterator(static_cast<internal::ParticleIteratorInterfaceImpl<Particle>*>(particleIteratorInterface)) {}

  /**
   * copy operator
   * @param otherParticleIteratorWrapper the other ParticleIteratorWrapper
   */
  ParticleIteratorWrapper(const ParticleIteratorWrapper& otherParticleIteratorWrapper)
      : _particleIterator(otherParticleIteratorWrapper._particleIterator->clone()) {}

  ParticleIteratorWrapper& operator=(const ParticleIteratorWrapper& otherParticleIteratorWrapper) {
    _particleIterator = std::unique_ptr<internal::ParticleIteratorInterfaceImpl<Particle>>(
        otherParticleIteratorWrapper._particleIterator->clone());
    return *this;
  }

  inline ParticleIteratorWrapper<Particle>& operator++() override {
    _particleIterator->operator++();
    return *this;
  }

  Particle& operator*() const override { return _particleIterator->operator*(); }

  void deleteCurrentParticle() override { _particleIterator->deleteCurrentParticle(); }

  bool isValid() const override { return _particleIterator->isValid(); }

 private:
  std::unique_ptr<internal::ParticleIteratorInterfaceImpl<Particle>> _particleIterator;
};

} /* namespace autopas */