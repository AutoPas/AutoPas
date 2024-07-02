/**
 * @file InteractionTypeOption.h
 * @author M. Muehlhaeusser
 * @date 14.08.2023
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {

/**
 * Class representing the interaction type choices.
 */
class InteractionTypeOption : public Option<InteractionTypeOption> {
 public:
  /**
   * Possible choices for the interaction type.
   */
  enum Value {
    /**
     * Pairwise interactions.
     */
    pairwise,
    /**
     * 3-Body interactions.
     */
    triwise,
  };

  /**
   * Constructor.
   */
  InteractionTypeOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr InteractionTypeOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @return
   */
  static std::set<InteractionTypeOption> getDiscouragedOptions() { return {}; }

  /**
   * Provides a way to iterate over the possible choices of InteractionTypeOption.
   * @return map option -> string representation
   */
  static std::map<InteractionTypeOption, std::string> getOptionNames() {
    return {
        {InteractionTypeOption::pairwise, "pairwise"},
        {InteractionTypeOption::triwise, "triwise"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace autopas
