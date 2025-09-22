/**
 * @file TuningTriggerInterface.h
 * @author Niklas Ladurner
 * @date 06.06.2025
 */

#pragma once

#include "autopas/options/TuningTriggerOption.h"

namespace autopas {

/**
 * Abstract class representing a trigger for dynamic tuning intervals.
 */
class TuningTriggerInterface {
 public:
  virtual ~TuningTriggerInterface() = default;

  /**
   * Checks whether a new tuning phase should be triggered in the current iteration.
   *
   * @param currentIteration The current iteration.
   * @param tuningInterval The tuningInterval in the currentIteration.
   * @return Boolean value signalling whether a new tuning phase should be triggered in the current iteration.
   */
  virtual bool shouldStartTuningPhase([[maybe_unused]] size_t currentIteration,
                                      [[maybe_unused]] size_t tuningInterval) = 0;

  /**
   * Gives the return value of shouldStartTuningPhase() for the next iteration.
   *
   * This is used to look ahead if a new tuning phase might start in the next iteration, such that containers may be
   * adapted.
   *
   * @param currentIteration The current iteration.
   * @param tuningInterval The tuningInterval in the currentIteration.
   * @return Boolean value signalling whether shouldStartTuningPhase() returns true in the next iteration.
   */
  inline bool shouldStartTuningPhaseNextIteration([[maybe_unused]] size_t currentIteration,
                                                  [[maybe_unused]] size_t tuningInterval) const {
    return (currentIteration = _nextTriggeringIteration - 1);
  };

  /**
   * Checks wether a new tuning phase will be triggered in less than 10 iterations from now.
   *
   * This is used to facilitate liveinfo collection before a new tuning phase starts.
   *
   * @param currentIteration The current iteration.
   * @param tuningInterval The tuningInterval in the currentIteration.
   * @return Boolean value signalling whether shouldStartTuningPhase() returns true in less than 10 iterations from now.
   */
  inline bool shouldStartTuningPhaseSoon([[maybe_unused]] size_t currentIteration,
                                         [[maybe_unused]] size_t tuningInterval) const {
    return (currentIteration > _nextTriggeringIteration - 10);
  };

  /**
   * Passes a singular iteration runtime to the trigger.
   *
   * @param sample The runtime of a singular iteration.
   */
  virtual void passRuntimeSample(unsigned long sample) = 0;

  /**
   * Get this object's associated TuningTriggerOption type.
   * @return TuningTriggerOption
   */
  virtual TuningTriggerOption getOptionType() const = 0;

 protected:
  unsigned long _nextTriggeringIteration = 10;
};
}  // namespace autopas
