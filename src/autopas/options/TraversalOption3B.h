/**
 * @file TriwiseTraversalOption.h
 * @author S. Newcome
 * @date 31/05/2023
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the traversal choices for 3-body interactions.
 */
class TraversalOption3B : public Option<TraversalOption3B> {
 public:
  /**
   * Possible choices for 3-body traversals. Try to maintain lexicographic ordering.
   */
  enum Value {
    // DirectSum Traversals:
    /**
     * DSSequentialTripletTraversal : Sequential triple loop over all particles.
     */
    ds_sequential_3b
  };

  /**
   * Constructor.
   */
  TraversalOption3B() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr TraversalOption3B(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<TraversalOption3B> getDiscouragedOptions() {
    return {};
  }

  /**
   * Provides a way to iterate over the possible choices of TraversalOption.
   * @return map option -> string representation
   */
  static std::map<TraversalOption3B, std::string> getOptionNames() {
    return {
        // DirectSum Traversals:
        {TraversalOption3B::ds_sequential_3b, "ds_sequential_3b"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
