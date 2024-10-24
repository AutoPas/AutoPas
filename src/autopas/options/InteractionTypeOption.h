/**
 * @file InteractionTypeOption.h
 * @author M. Muehlhaeusser
 * @date 14.08.2023
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Class representing the interaction type choices.
 */
class InteractionTypeOption : public Option<InteractionTypeOption> {
 public:
  /**
   * Type used for the internal enum.
   */
  using Value_t = unsigned int;

  /**
   * Possible choices for the interaction type.
   */
  enum Value : Value_t {
    /**
     * Pairwise interactions.
     */
    pairwise = 0b0001,
    /**
     * Triwise interactions.
     */
    triwise = 0b0010,
    /**
     * All interaction types. Used e.g. for setting autopas options for all interactions.
     */
    all = 0b0011,
  };

  // sanity check
  // combined values should be equivalent to InteractionTypeOption::all.
  static_assert((pairwise | triwise) == all, "InteractionTypeOptions are defined with non matching values!");

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
   * Constructor from number value.
   * This is useful when combining values and directly using the result as argument of type InteractionTypeOption.
   * This is necessary since e.g. pairwise & triwise results in an object of Value_t instead of Value.
   * @param option
   */
  constexpr InteractionTypeOption(Value_t option) : _value(static_cast<InteractionTypeOption::Value>(option)) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Set of options that are very unlikely to be interesting.
   * @note Here, this specifically means InteractionTypeOption::all
   * @return
   */
  static std::set<InteractionTypeOption> getDiscouragedOptions() { return {Value::all}; }

  /**
   * Provides a way to iterate over the possible choices of InteractionTypeOption.
   * @return map option -> string representation
   */
  static std::map<InteractionTypeOption, std::string> getOptionNames() {
    return {
        {InteractionTypeOption::pairwise, "pairwise"},
        {InteractionTypeOption::triwise, "triwise"},
        {InteractionTypeOption::all, "all"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
