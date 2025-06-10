/**
 * @file ComputationLoadOption.h
 * @author Ashutosh Solanki
 * @date 2025-06-10
 */

#pragma once

#include <map>
#include <set>

#include "autopas/options/Option.h"

/**
 * Class representing the computation load options for load balancing.
 */
class ComputationLoadOption : public autopas::Option<ComputationLoadOption> {
 public:
  /**
   * Possible choices for computation load metrics.
   */
  enum Value {
    /**
     * Use complete 1 cycle time as computation load.
     */
    completeCycle,
    /**
     * Use position update time as computation load.
     */
    positionUpdate,
    /**
     * Use container update time as computation load.
     */
    updateContainer,
    /**
     * Use halo particle exchange time as computation load.
     */
    haloParticleExchange,
    /**
     * Use migrating particle exchange time as computation load.
     */
    migratingParticleExchange,
    /**
     * Use reflection at boundaries time as computation load.
     */
    reflectParticlesAtBoundaries,
    /**
     * Use total force update time as computation load.
     */
    forceUpdateTotal,
    /**
     * Use velocity update time as computation load.
     */
    velocityUpdate
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
   * Provides a way to iterate over the possible choices of ComputationLoadOption.
   * @return map option -> string representation
   */
  static std::map<ComputationLoadOption, std::string> getOptionNames() {
    return {{ComputationLoadOption::completeCycle, "CompleteCycle"},
            {ComputationLoadOption::positionUpdate, "PositionUpdate"},
            {ComputationLoadOption::updateContainer, "UpdateContainer"},
            {ComputationLoadOption::haloParticleExchange, "HaloParticleExchange"},
            {ComputationLoadOption::migratingParticleExchange, "MigratingParticleExchange"},
            {ComputationLoadOption::reflectParticlesAtBoundaries, "ReflectParticlesAtBoundaries"},
            {ComputationLoadOption::forceUpdateTotal, "ForceUpdateTotal"},
            {ComputationLoadOption::velocityUpdate, "VelocityUpdate"}};
  }

 private:
  Value _value{Value::completeCycle};
}; 