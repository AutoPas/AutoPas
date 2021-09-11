/**
 * @file LoadBalancerOption.h
 * @author J. Körner
 * @date 08.09.2021
 */
#pragma once

#include "autopas/options/Option.h"

class LoadBalancerOption : public autopas::Option<LoadBalancerOption> {
 public:
  /**
   * Different load balancer options
   */
  enum Value {
    /**
     * Use the inverted pressure load balancer.
     */
    invertedPressure,

    /**
     * Use the ALL load balance.
     */
    all,

    /** * Don't use any load balancer.
     */
    none
  };

  /**
   * Constructor.
   */
  LoadBalancerOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr LoadBalancerOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provides a way to iterate over the possible choices of LoadBalancerOptions.
   * @return map option -> string representation
   */
  static std::map<LoadBalancerOption, std::string> getOptionNames() {
    return {{LoadBalancerOption::invertedPressure, "Inverted Pressure"},
            {LoadBalancerOption::all, "ALL"},
            {LoadBalancerOption::none, "None"}};
  }

 private:
  /**
   * Stores the value of the option.
   */
  Value _value{Value(-1)};
};
