/**
 * @file ComputationLoadOption.h
 * @author Ashutosh Solanki
 * @date 02.07.2025
 */
#pragma once

#include "autopas/options/Option.h"

/**
 * Class representing the choices for computation load tracking during the MD-Flexible simulation.
 * Users can choose between different computation load metrics to track.
 */
class ComputationLoadOption : public autopas::Option<ComputationLoadOption> {
 public:
  /**
   * Different computation load options
   */
  enum Value {
    /**
     * Track complete cycle computation load.
     */
    completeCycle,

    /**
     * Track force update computation load.
     */
    forceUpdate,

    /**
     * Track MPI communication computation load.
     */
    MPICommunication,

    /**
     * Track particle count computation load.
     */
    particleCount
  };

  /**
   * Constructor.
   */
  ComputationLoadOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr ComputationLoadOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provides a way to iterate over the possible choices of ComputationLoadOptions.
   * @return map option -> string representation
   */
  static std::map<ComputationLoadOption, std::string> getOptionNames() {
    return {{ComputationLoadOption::completeCycle, "CompleteCycle"},
            {ComputationLoadOption::forceUpdate, "ForceUpdate"},
            {ComputationLoadOption::MPICommunication, "MPICommunication"},
            {ComputationLoadOption::particleCount, "ParticleCount"}};
  }

 private:
  /**
   * Stores the value of the option.
   */
  Value _value{Value(-1)};
}; 