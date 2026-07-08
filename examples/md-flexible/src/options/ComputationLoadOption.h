/**
 * @file ComputationLoadOption.h
 * @author Ashutosh Solanki
 * @date 02.07.2025
 */
#pragma once

#include "autopas/options/Option.h"

/**
 * Class representing the choices for computational load tracking during the MD-Flexible simulation.
 * Users can choose between different computational load metrics to track.
 */
class ComputationLoadOption : public autopas::Option<ComputationLoadOption> {
 public:
  /**
   * Different computation load options
   */
  enum Value {
    /**
     * Use the complete cycle as the computational load.
     */
    completeCycle,

    /**
     * Use the non-boundary calculations as the computational load.
     */
    nonBoundaryCalculations,

    /**
     * Use the force update as the computational load.
     */
    forceUpdate,

    /**
     * Use the MPI communication as the computational load.
     */
    MPICommunication,

    /**
     * Use the particle count as the computational load.
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
    return {{completeCycle, "CompleteCycle"},
            {nonBoundaryCalculations, "NonBoundaryCalculations"},
            {forceUpdate, "ForceUpdate"},
            {MPICommunication, "MPICommunication"},
            {particleCount, "ParticleCount"}};
  }

 private:
  /**
   * Stores the value of the option.
   */
  Value _value{Value(-1)};
};