/**
 * @file EnergySensorOption.h
 * @author Maximilian Praus
 * @date 06.07.2024
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {

/**
 * Class representing the different energy sensor options
 */

class EnergySensorOption : public Option<EnergySensorOption> {
 public:
  /**
   * Possible choices for energy sensor
   */
  enum Value {
    /**
     * Use LIKWID
     */
    likwid,
    /**
     * Use RAPL
     */
    rapl,
  };

  /**
   * Cosntructor
   */
  EnergySensorOption() = default;

  /**
   * Cosntructor with selected option
   * @param option
   */
  constexpr EnergySensorOption(Value option) : _value(option) {}

  /**
   * Cast to value
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provide a way to iterate over the options of EnergySensorOption
   * @return map option -> string representation
   */
  static std::map<EnergySensorOption, std::string> getOptionNames() {
    return {
        {EnergySensorOption::likwid, "likwid"},
        {EnergySensorOption::rapl, "rapl"},
    };
  };

 private:
  Value _value{Value(-1)};
};

}  // namespace autopas