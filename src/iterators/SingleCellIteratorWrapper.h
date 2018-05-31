/**
 * SingleCellIteratorWrapper.h
 *
 *  Created on: 31 May 2018
 *      Author: seckler
 */
#pragma once

#include <utils/ExceptionHandler.h>
#include "SingleCellIteratorInterface.h"

namespace autopas {

/**
 * ParticleIteratorInterface class to iterate over particles.
 * This class provides a basic interface for all iterators within AutoPas.
 * The particles can be accessed using "iterator->" or "*iterator". The next
 * particle using the ++operator, e.g. "++iterator"
 * @tparam Particle type of the particle that is accessed
 */
template <class Particle>
class SingleCellIteratorWrapper : public SingleCellIteratorInterface<Particle> {
 public:
  /**
   * Constructor of the ParticleIteratorWrapper
   * @param particleIteratorInterface
   * @tparam InterfacePtrType type of the interface ptr
   */
  template <class InterfacePtrType>
  SingleCellIteratorWrapper(InterfacePtrType* particleIteratorInterface)
      : _particleIterator(
            static_cast<internal::SingleCellIteratorInterfaceImpl<Particle>*>(
                particleIteratorInterface)) {}

  /**
   * copy constructor
   * @param otherParticleIteratorWrapper the other ParticleIteratorWrapper
   */
  SingleCellIteratorWrapper(
      const SingleCellIteratorWrapper& otherParticleIteratorWrapper) {
    _particleIterator =
        std::unique_ptr<internal::SingleCellIteratorInterfaceImpl<Particle>>(
            otherParticleIteratorWrapper._particleIterator->clone());
  }

  /**
   * copy constructor
   * @param otherParticleIteratorWrapper the other ParticleIteratorWrapper
   */
  SingleCellIteratorWrapper& operator=(
      const SingleCellIteratorWrapper& otherParticleIteratorWrapper) {
    _particleIterator =
        std::unique_ptr<internal::SingleCellIteratorInterfaceImpl<Particle>>(
            otherParticleIteratorWrapper._particleIterator->clone());
    return *this;
  }

  inline SingleCellIteratorWrapper<Particle>& operator++() override {
    _particleIterator->operator++();
    return *this;
  }

  Particle& operator*() const override {
    return _particleIterator->operator*();
  }

  void deleteCurrentParticle() override {
    _particleIterator->deleteCurrentParticle();
  }

  bool isValid() const override { return _particleIterator->isValid(); }

  bool operator==(
      const SingleCellIteratorInterface<Particle>& rhs) const override {
    return _particleIterator->operator==(rhs);
  }

  bool operator!=(
      const SingleCellIteratorInterface<Particle>& rhs) const override {
    return _particleIterator->operator!=(rhs);
  }

  int getIndex() const override { return _particleIterator->getIndex(); }

 private:
  std::unique_ptr<internal::SingleCellIteratorInterfaceImpl<Particle>>
      _particleIterator;
};

} /* namespace autopas */