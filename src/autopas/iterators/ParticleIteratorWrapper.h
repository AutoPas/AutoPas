/**
 * @file ParticleIteratorWrapper.h
 * Contains ParticleIteratorWrapper class.
 * @author seckler
 * @date 29.05.2018
 */
#pragma once

#include <memory>

#include "autopas/iterators/ParticleIteratorInterface.h"
#include "autopas/utils/ExceptionHandler.h"

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
 * @tparam modifiable Defines whether the ParticleIterator is modifiable or not. If it is false, it points to a const
 * Particle.
 */
template <class Particle, bool modifiable>
class ParticleIteratorWrapper : public ParticleIteratorInterface<Particle, modifiable> {
  using ParticleType = std::conditional_t<modifiable, Particle, const Particle>;

 public:
  ParticleIteratorWrapper() : _particleIterator(nullptr) {}

  /**
   * Constructor of the ParticleIteratorWrapper
   * @param particleIteratorInterface
   * @tparam InterfacePtrType type of the interface ptr
   */
  template <class InterfacePtrType>
  explicit ParticleIteratorWrapper(InterfacePtrType *particleIteratorInterface)
      : _particleIterator(
            static_cast<internal::ParticleIteratorInterfaceImpl<Particle, modifiable> *>(particleIteratorInterface)) {}

  /**
   * copy operator
   * @param otherParticleIteratorWrapper the other ParticleIteratorWrapper
   */
  ParticleIteratorWrapper(const ParticleIteratorWrapper &otherParticleIteratorWrapper) : _particleIterator(nullptr) {
    if (otherParticleIteratorWrapper._particleIterator) {
      _particleIterator.reset(otherParticleIteratorWrapper._particleIterator->clone());
    }
  }

  /**
   * copy assignment operator to assign contents of otherParticleIteratorWrapper to this ParticleIteratorWrapper
   * @param otherParticleIteratorWrapper the other IteratorWrapper
   * @return *this
   */
  inline ParticleIteratorWrapper &operator=(const ParticleIteratorWrapper &otherParticleIteratorWrapper) {
    _particleIterator = std::unique_ptr<internal::ParticleIteratorInterfaceImpl<Particle, modifiable>>(
        otherParticleIteratorWrapper._particleIterator->clone());
    return *this;
  }

  inline ParticleIteratorWrapper<Particle, modifiable> &operator++() override final {
    _particleIterator->operator++();
    return *this;
  }

  inline ParticleType &operator*() const override final { return _particleIterator->operator*(); }

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

  /**
   * Adds an additional particle vector to the iterator.
   * @param additionalVector The additional particle vector that should also be iterated over.
   */
  void addAdditionalVector(
      std::conditional_t<modifiable, std::vector<Particle> &, const std::vector<Particle> &> additionalVector) {
    _particleIterator->addAdditionalVector(additionalVector);
  }

 protected:
  inline void deleteCurrentParticleImpl() override final {
    if constexpr (modifiable) {
      internal::deleteParticle(*_particleIterator);
    } else {
      utils::ExceptionHandler::exception("Error: Trying to delete a particle through a const iterator.");
    }
  }

 private:
  std::unique_ptr<autopas::internal::ParticleIteratorInterfaceImpl<Particle, modifiable>> _particleIterator;
};

/**
 * Provides type aliases for iterator types, to provide a fixed interface that is not influenced by changes in the
 * iterator classes.
 * @tparam Particle The type of the particle.
 */
template <typename Particle>
class IteratorTraits {
 public:
  /**
   * Alias for normal (non-const) iterator.
   */
  using iterator_t = autopas::ParticleIteratorWrapper<Particle, true>;
  /**
   * Alias for const iterator.
   */
  using const_iterator_t = autopas::ParticleIteratorWrapper<Particle, false>;
};

}  // namespace autopas