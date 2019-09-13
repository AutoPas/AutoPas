/**
 * @file SingleCellIteratorWrapper.h
 *
 * @date 31 May 2018
 * @author seckler
 */
#pragma once

#include "autopas/iterators/SingleCellIteratorInterface.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * SingleCellIteratorWrapper class is the main class visible to the user to iterate over particles in one cell.
 * The particles can be accessed using "iterator->" or "*iterator". The next particle using the ++operator, e.g.
 * "++iteratorWrapper".
 * @note The wrapper class provides an easy way to access the functions of the SingleCellIteratorInterface, e.g.
 * allowing "iteratorWrapper->getR()" or "++iteratorWrapper". Without the wrapper class and using the
 * ParticleIteratorInterface the calls would require de-referencing like:
 * (*iteratorInterface)->getR() or "++(*iteratorWrapper)"
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
  explicit SingleCellIteratorWrapper(InterfacePtrType *particleIteratorInterface)
      : _particleIterator(
            static_cast<internal::SingleCellIteratorInterfaceImpl<Particle> *>(particleIteratorInterface)) {}

  /**
   * copy constructor
   * @param otherParticleIteratorWrapper the other ParticleIteratorWrapper
   */
  SingleCellIteratorWrapper(const SingleCellIteratorWrapper &otherParticleIteratorWrapper) {
    _particleIterator = std::unique_ptr<internal::SingleCellIteratorInterfaceImpl<Particle>>(
        otherParticleIteratorWrapper._particleIterator->clone());
  }

  /**
   * copy assignment constructor
   * @param otherParticleIteratorWrapper the other ParticleIteratorWrapper
   * @return the modified SingleCellIteratorWrapper
   */
  SingleCellIteratorWrapper &operator=(const SingleCellIteratorWrapper &otherParticleIteratorWrapper) {
    _particleIterator = std::unique_ptr<internal::SingleCellIteratorInterfaceImpl<Particle>>(
        otherParticleIteratorWrapper._particleIterator->clone());
    return *this;
  }

  inline SingleCellIteratorWrapper<Particle> &operator++() override final {
    _particleIterator->operator++();
    return *this;
  }

  inline Particle &operator*() const override final { return _particleIterator->operator*(); }

  inline void deleteCurrentParticle() override final { _particleIterator->deleteCurrentParticle(); }

  inline bool isValid() const override final { return _particleIterator->isValid(); }

  inline bool operator==(const SingleCellIteratorInterface<Particle> &rhs) const override final {
    return _particleIterator->operator==(rhs);
  }

  inline bool operator!=(const SingleCellIteratorInterface<Particle> &rhs) const override final {
    return _particleIterator->operator!=(rhs);
  }

  inline size_t getIndex() const override final { return _particleIterator->getIndex(); }

  /**
   * Returns the stored single cell iterator.
   * @return
   */
  inline SingleCellIteratorInterface<Particle> *get() const { return _particleIterator.get(); }

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
  std::unique_ptr<internal::SingleCellIteratorInterfaceImpl<Particle>> _particleIterator;
};

}  // namespace autopas