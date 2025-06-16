/**
 * @file StaticSimpleTrigger.h
 * @author Niklas Ladurner
 * @date 06.06.2025
 */

#pragma once

#include "autopas/tuning/triggers/TuningTriggerInterface.h"

namespace autopas {

/**
 * Represents a trigger with the same behavior as static tuning, i.e. it triggers a tuning phase on each
 * multiple of tuningInterval.
 */
class StaticSimpleTrigger : public TuningTriggerInterface {
 public:
  /**
   * Constructor
   */
  StaticSimpleTrigger() = default;

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) override {
    bool oldTriggerState = _wasTriggered;
    _wasTriggered = (currentIteration % tuningInterval == 0);
    return oldTriggerState;
  }

  inline bool needsRuntimeSample() const override { return false; }
  void passRuntimeSample(unsigned long sample) override {}

  TuningTriggerOption getOptionType() const override { return TuningTriggerOption::staticSimple; }
};
}  // namespace autopas
