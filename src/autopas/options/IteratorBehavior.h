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
 * Class representing the acquisition function choices for the Bayesian search.
 */
class IteratorBehavior : public Option<IteratorBehavior> {
 public:
  /**
   * Different possibilities for iterator behaviors
   */
  enum Value {
    /**
     * Iterate only over halo particles.
     */
    haloOnly,
    /**
     * Iterate only over owned particles.
     */
    ownedOnly,
    /**
     * Iterate over both halo and owned particles.
     */
    haloAndOwned,
    /**
     * Iterate over both halo and owned particles and also dummy particles.
     */
    haloOwnedAndDummy,
  };

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
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<IteratorBehavior> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of AcquisitionFunction.
   * @return map option -> string representation
   */
  static std::map<IteratorBehavior, std::string> getOptionNames() {
    return {
        {IteratorBehavior::ownedOnly, "ownedOnly"},
        {IteratorBehavior::haloOnly, "haloOnly"},
        {IteratorBehavior::haloAndOwned, "haloAndOwned"},
        {IteratorBehavior::haloOwnedAndDummy, "haloOwnedAndDummy"},
    };
  }

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
