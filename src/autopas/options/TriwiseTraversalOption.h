/**
 * @file TraversalOption.h
 * @author F. Gratl
 * @date 1/18/19
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the traversal choices for 3-body interactions.
 */
class TriwiseTraversalOption : public Option<TriwiseTraversalOption> {
 public:
  /**
   * Possible choices for 3-body traversals. Try to maintain lexicographic ordering.
   */
  enum Value {
    // DirectSum Traversals:
    /**
     * DSSequentialTripletTraversal : Sequential triple loop over all particles.
     */
    ds_sequential_triplet
  };

  /**
   * Constructor.
   */
  TriwiseTraversalOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr TriwiseTraversalOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<TriwiseTraversalOption> getDiscouragedOptions() {
    return {};
  }

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   * @return map option -> string representation
   */
  static std::map<TriwiseTraversalOption, std::string> getOptionNames() {
    return {
        // DirectSum Traversals:
        {TriwiseTraversalOption::ds_sequential_triplet, "ds_sequential_triplet"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
