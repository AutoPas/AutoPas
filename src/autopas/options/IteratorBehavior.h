/**
 * @file IteratorBehavior.h
 * @author F. Gratl
 * @date 25.02.21
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"
#include "autopas/particles/OwnershipState.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the choices for behaviors of iterators.
 * Choices are a bit vector and can thus be combined via a logical OR.
 */
class IteratorBehavior : public Option<IteratorBehavior> {
 public:
  /**
   * Type used for the internal enum.
   */
  using Value_t = unsigned int;

  /**
   * Different possibilities for iterator behaviors.
   */
  enum Value : Value_t {
    /**
     * Iterate only over owned particles.
     */
    owned = 0b00001,
    /**
     * Iterate only over halo particles.
     */
    halo = 0b00010,
    /**
     * Iterate over both halo and owned particles. Defined fore ease of access.
     */
    ownedOrHalo = 0b00011,
    /**
     * Iterate only over dummy particles.
     */
    dummy = 0b00100,
    /**
     * Iterate over both halo and owned particles and also dummy particles. Defined fore ease of access.
     */
    ownedOrHaloOrDummy = 0b00111,
    /**
     * Force the iterator to behave like a sequential iterator even when created in a parallel region.
     */
    forceSequential = 0b01000,
    /**
     * Force the iterator to iterate over Container only
     */
    containerOnly = 0b10000,
  };

  // sanity checks
  // combined values should be equivalent to the disjunction of their elements
  static_assert((owned | halo) == ownedOrHalo, "Iterator behaviors are defined with non matching values!");
  static_assert((owned | halo | dummy) == ownedOrHaloOrDummy,
                "Iterator behaviors are defined with non matching values!");
  // forceSequential must not overlap with anything else
  static_assert((ownedOrHaloOrDummy & forceSequential) == 0,
                "Iterator behaviors are defined with non matching values!");
  // containerOnly must not overlap with anything else
  static_assert(((ownedOrHaloOrDummy | forceSequential) & containerOnly) == 0,
                "Iterator behaviors are defined with non matching values!");

  /**
   * Constructor.
   */
  IteratorBehavior() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr IteratorBehavior(Value option) : _value(option) {}

  /**
   * Constructor from number value.
   * This is useful when combining values and directly using the result as argument of type IteratorBehavior.
   * This is necessary since e.g. owned & halo results in an object of Value_t instead of Value.
   * @param option
   */
  constexpr IteratorBehavior(Value_t option) : _value(static_cast<IteratorBehavior::Value>(option)) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * Technically these options are not discouraged but when fetching a set of desired behaviors these are probably
   * not useful and would be excluded anyway.
   * @return
   */
  static std::set<IteratorBehavior> getDiscouragedOptions() {
    return {IteratorBehavior::dummy, IteratorBehavior::ownedOrHaloOrDummy, IteratorBehavior::forceSequential,
            IteratorBehavior::containerOnly};
  }

  /**
   * Provides a way to iterate over the possible choices of AcquisitionFunction.
   * @return map option -> string representation
   */
  static std::map<IteratorBehavior, std::string> getOptionNames() {
    return {
        {IteratorBehavior::owned, "owned"},
        {IteratorBehavior::halo, "halo"},
        {IteratorBehavior::ownedOrHalo, "ownedOrHalo"},
        {IteratorBehavior::containerOnly, "containerOnly"},
        {IteratorBehavior::dummy, "dummy"},
        {IteratorBehavior::ownedOrHaloOrDummy, "ownedOrHaloOrDummy"},
        {IteratorBehavior::forceSequential, "forceSequential"},
    };
  }

  /**
   * Check whether this iterator behavior covers the given particle
   * @tparam ParticleType
   * @param particle particle to be checked
   * @return true if this iterator behavior covers the given particle, false otherwise
   */
  template <typename ParticleType>
  bool contains(ParticleType &particle) {
    switch (this->_value) {
      case options::IteratorBehavior::ownedOrHaloOrDummy:
        return true;
      case options::IteratorBehavior::ownedOrHalo:
        return not particle.isDummy();
      case options::IteratorBehavior::halo:
        return particle.isHalo();
      case options::IteratorBehavior::owned:
        return particle.isOwned();
      case options::IteratorBehavior::dummy:
        return particle.isDummy();
      default:
        utils::ExceptionHandler::exception("IteratorBehavior::contains() Unknown IteratorBehavior: {}.",
                                           getOptionNames()[this->_value]);
        return false;
    }
  }

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
