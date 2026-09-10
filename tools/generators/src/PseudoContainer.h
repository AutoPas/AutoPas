/**
 * @file PseudoContainer.h
 * @author F. Gratl
 * @date 01.08.24
 */

#pragma once

namespace autopasTools {
/**
 * Wrapper class for things like std::vector to be used in templated methods that take an AutoPas or ParticleContainer.
 *
 * This wrapper for example maps functions like std::vector::push_back to AutoPas::addParticle().
 *
 * It does not take ownership of the wrapped object but only holds a pointer to it!
 *
 * @tparam T Full type of the wrapped container
 */
template <class T>
class PseudoContainer {
 public:
  /**
   * Constructor
   * @param actualContainer Reference to container to be wrapped.
   */
  explicit PseudoContainer(T &actualContainer) : actualContainer(&actualContainer) {}

  /**
   * Actual type of the container without qualifiers or pointers.
   */
  using ActualT = std::remove_pointer_t<std::decay_t<T>>;
  /**
   * Particle / value type.
   */
  using ParticleType = typename ActualT::value_type;

  /**
   * Wrapper for push_back
   * @param p
   */
  void addParticle(const ParticleType &p) { actualContainer->push_back(p); }

 private:
  /**
   * Pointer to the wrapped container
   */
  ActualT *actualContainer;
};
}  // namespace autopasTools