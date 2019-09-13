/**
 * @file ParticleIteratorWrapper.h
 * Contains ParticleIteratorWrapper class.
 * @author seckler
 * @date 29.05.2018
 */
#pragma once

#include <memory>
#include "autopas/iterators/ParticleIteratorInterface.h"

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
  ParticleIteratorWrapper() : _particleIterator(nullptr) {}

  /**
   * Constructor of the ParticleIteratorWrapper
   * @param particleIteratorInterface
   * @tparam InterfacePtrType type of the interface ptr
   */
  template <class InterfacePtrType>
  explicit ParticleIteratorWrapper(InterfacePtrType *particleIteratorInterface)
      : _particleIterator(static_cast<internal::ParticleIteratorInterfaceImpl<Particle> *>(particleIteratorInterface)) {
  }

  /**
   * copy operator
   * @param otherParticleIteratorWrapper the other ParticleIteratorWrapper
   */
  ParticleIteratorWrapper(const ParticleIteratorWrapper &otherParticleIteratorWrapper)
      : _particleIterator(otherParticleIteratorWrapper._particleIterator->clone()) {}

  /**
   * copy assignment operator to assign contents of otherParticleIteratorWrapper to this ParticleIteratorWrapper
   * @param otherParticleIteratorWrapper the other IteratorWrapper
   * @return *this
   */
  inline ParticleIteratorWrapper &operator=(const ParticleIteratorWrapper &otherParticleIteratorWrapper) {
    _particleIterator = std::unique_ptr<internal::ParticleIteratorInterfaceImpl<Particle>>(
        otherParticleIteratorWrapper._particleIterator->clone());
    return *this;
  }

  inline ParticleIteratorWrapper<Particle> &operator++() override final {
    _particleIterator->operator++();
    return *this;
  }

  inline Particle &operator*() const override final { return _particleIterator->operator*(); }

  inline void deleteCurrentParticle() override final { _particleIterator->deleteCurrentParticle(); }

  inline bool isValid() const override final {
    if (_particleIterator) {
      return _particleIterator->isValid();
    } else {
      return false;
    }
  }

  /**
   * Equality operator that compares with a bool.
   * Needed to be able to compare with AutoPas::end().
   * @param input normally: AutoPas::end()
   * @return true if isValid == input, false otherwise.
   */
  bool operator==(const bool &input) const { return isValid() == input; }

  /**
   * Inequality operator that compares with a bool.
   * Needed to be able to compare with end().
   * @param input normally: autoPas.end()
   * @return true if isValid != input, false otherwise.
   */
  bool operator!=(const bool &input) const { return not(*this == input); }

 private:
  std::unique_ptr<autopas::internal::ParticleIteratorInterfaceImpl<Particle>> _particleIterator;
};

}  // namespace autopas