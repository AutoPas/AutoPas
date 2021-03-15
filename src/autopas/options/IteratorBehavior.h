/**
 * @file IteratorBehavior.h
 * @author F. Gratl
 * @date 25.02.21
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

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
    owned = 0b0001,
    /**
     * Iterate only over halo particles.
     */
    halo = 0b0010,
    /**
     * Iterate over both halo and owned particles. Defined fore ease of access.
     */
    ownedOrHalo = 0b0011,
    /**
     * Iterate only over dummy particles.
     */
    dummy = 0b0100,
    /**
     * Iterate over both halo and owned particles and also dummy particles. Defined fore ease of access.
     */
    ownedOrHaloOrDummy = 0b0111,
    /**
     * Force the iterator to behave like a sequential iterator even when created in a parallel region.
     */
    forceSequential = 0b1000,
  };

  // sanity checks
  // combined values should be equivalent to the disjunction of their elements
  static_assert((owned | halo) == ownedOrHalo, "Iterator behaviors are defined with non matching values!");
  static_assert((owned | halo | dummy) == ownedOrHaloOrDummy,
                "Iterator behaviors are defined with non matching values!");
  // forceSequential must does not overlap with anything else
  static_assert((ownedOrHaloOrDummy & forceSequential) == 0,
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
   * Technically these options are nor discouraged but when fetching a set of desired behaviors these are probably
   * not useful and would be excluded anyway.
   * @return
   */
  static std::set<IteratorBehavior> getDiscouragedOptions() {
    return {IteratorBehavior::dummy, IteratorBehavior::ownedOrHaloOrDummy, IteratorBehavior::forceSequential};
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
        {IteratorBehavior::dummy, "dummy"},
        {IteratorBehavior::ownedOrHaloOrDummy, "ownedOrHaloOrDummy"},
        {IteratorBehavior::forceSequential, "forceSequential"},
    };
  }

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
