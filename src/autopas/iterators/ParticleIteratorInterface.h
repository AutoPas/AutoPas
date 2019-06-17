/**
 * @file ParticleIteratorInterface.h
 *
 * @date 29 May 2018
 * @author seckler
 */
#pragma once

namespace autopas {

/**
 * Enum to specify the behavior of an iterator.
 */
enum IteratorBehavior {
  haloOnly,     /// iterate only over halo
  ownedOnly,    /// iterate only over inner cells
  haloAndOwned  /// iterate over both halo and inner cells
};

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
  virtual ParticleIteratorInterface<Particle> &operator++() = 0;

  virtual ~ParticleIteratorInterface(){};

  /**
   * access the particle using *iterator
   * this is the indirection operator
   * @return current particle
   */
  virtual Particle &operator*() const = 0;

  /**
   * access particle using iterator->
   *
   * this is the member of pointer operator
   * @return current particle
   */
  virtual inline Particle *operator->() const final { return &(this->operator*()); }

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

inline namespace internal {
/**
 * All implementations of the interface should inherit from this class. It extends the interface just by the clone
 * method, which is needed by the Wrapper.
 * @tparam Particle
 */
template <class Particle>
class ParticleIteratorInterfaceImpl : public ParticleIteratorInterface<Particle> {
 public:
  /**
   * Clones the current object, should allocate new object and return it.
   * @return the clone
   */
  virtual ParticleIteratorInterfaceImpl<Particle> *clone() const = 0;
};

}  // namespace internal

}  // namespace autopas