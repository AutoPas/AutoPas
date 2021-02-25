/**
 * @file ParticleIteratorInterface.h
 *
 * @date 29 May 2018
 * @author seckler
 */
#pragma once

namespace autopas {

namespace internal {

/**
 * Function to access hidden iterator.deleteCurrentParticle() to mark it as internal.
 * @tparam ParticleIterator
 * @param iterator
 */
template <class ParticleIterator>
void deleteParticle(ParticleIterator &iterator) {
  iterator.deleteCurrentParticle();
}
}  // namespace internal
/**
 * ParticleIteratorInterface class to iterate over particles.
 * This class provides a basic interface for all iterators within AutoPas.
 * The particles can be accessed using "iterator->" or "*iterator". The next
 * particle using the ++operator, e.g. "++iterator"
 * @tparam Particle type of the particle that is accessed
 * @tparam modifiable Defines whether the ParticleIterator is modifiable or not. If it is false, it points to a const
 * Particle.
 */
template <class Particle, bool modifiable>
class ParticleIteratorInterface {
  using ParticleType = std::conditional_t<modifiable, Particle, const Particle>;

  /**
   * Function to access hidden iterator.deleteCurrentParticle() to mark it as internal.
   * @tparam ParticleIterator
   */
  template <class T>
  friend void internal::deleteParticle(T &);

 public:
  /**
   * Increment operator.
   * Used to jump to the next particle
   * @return next particle, usually ignored
   */
  virtual ParticleIteratorInterface<Particle, modifiable> &operator++() = 0;

  virtual ~ParticleIteratorInterface() = default;

  /**
   * Access the particle using *iterator
   * This is the indirection operator
   * @return current particle
   */
  virtual ParticleType &operator*() const = 0;

  /**
   * Access particle using iterator->
   *
   * This is the member of pointer operator
   * @return current particle
   */
  virtual inline ParticleType *operator->() const final { return &(this->operator*()); }

  /**
   * Check whether the iterator is valid
   * @return returns whether the iterator is valid
   */
  virtual bool isValid() const = 0;

 protected:
  /**
   * Deletes the current particle.
   * @note Do NOT use this function from outside the autopas namespace. Use AutoPas::deleteParticle instead.
   *
   * @note From inside of autopas use internal::deleteParticle(iterator).
   *
   * @note This function is disabled for const iterators.
   * @return void
   */
  template <typename Dummy = void>
  inline std::enable_if_t<modifiable, Dummy> deleteCurrentParticle() {
    deleteCurrentParticleImpl();
  }

  /**
   * Implementation of the deletion. The split is needed to disable it.
   */
  virtual void deleteCurrentParticleImpl() = 0;
};

namespace internal {
/**
 * All implementations of the interface should inherit from this class. It extends the interface just by the clone
 * method, which is needed by the Wrapper.
 * @tparam Particle
 * @tparam modifiable
 */
template <class Particle, bool modifiable>
class ParticleIteratorInterfaceImpl : public ParticleIteratorInterface<Particle, modifiable> {
 public:
  /**
   * Clones the current object, should allocate new object and return it.
   * @return the clone
   */
  virtual ParticleIteratorInterfaceImpl<Particle, modifiable> *clone() const = 0;
};

}  // namespace internal

}  // namespace autopas
