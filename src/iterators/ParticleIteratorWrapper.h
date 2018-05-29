/**
 * ParticleIterator.h
 *
 *  Created on: 17 Jan 2018
 *      Author: seckler
 */
#pragma once

#include <utils/ExceptionHandler.h>
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
      : _particleIterator(static_cast<ParticleIteratorInterface<Particle>*>(
            particleIteratorInterface)) {}

  /**
   * copy operator
   * @param otherParticleIteratorWrapper the other ParticleIteratorWrapper
   */
  ParticleIteratorWrapper(
      const ParticleIteratorWrapper& otherParticleIteratorWrapper) {
    _particleIterator = std::unique_ptr<ParticleIteratorInterface<Particle>>(
        otherParticleIteratorWrapper._particleIterator->clone());
  }

  inline ParticleIteratorWrapper<Particle>& operator++() override {
    _particleIterator->operator++();
    return *this;
  }

  Particle& operator*() override { return _particleIterator->operator*(); }

  void deleteCurrentParticle() override {
    _particleIterator->deleteCurrentParticle();
  }

  bool isValid() const override { return _particleIterator->isValid(); }

  ParticleIteratorInterface<Particle>* clone() const override {
    autopas::utils::ExceptionHandler::exception(
        "ParticleIteratorWrapper::clone() should never be called");
    return nullptr;
  }

 private:
  std::unique_ptr<ParticleIteratorInterface<Particle>> _particleIterator;
};

} /* namespace autopas */