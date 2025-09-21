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
  virtual bool shouldStartTuningPhaseNextIteration([[maybe_unused]] size_t currentIteration,
                                                   [[maybe_unused]] size_t tuningInterval) const {
    return _wasTriggered;
  };

  /**
   * Passes a singular iteration runtime to the trigger.AutoPas_MPI_Cart_create
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
  bool _wasTriggered;
};
}  // namespace autopas
