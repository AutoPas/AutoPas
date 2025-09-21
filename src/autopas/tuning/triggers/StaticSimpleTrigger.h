/**
 * @file StaticSimpleTrigger.h
 * @author Niklas Ladurner
 * @date 06.06.2025
 */

#pragma once

#include "autopas/tuning/triggers/TuningTriggerInterface.h"

namespace autopas {

/**
 * Represents a trigger with the same behavior as static tuning intervals, i.e. it triggers a tuning phase on each
 * multiple of tuningInterval.
 */
class StaticSimpleTrigger : public TuningTriggerInterface {
 public:
  /**
   * Constructor
   */
  StaticSimpleTrigger() = default;

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) override {
    return (currentIteration % tuningInterval) == 0;
  }

  inline bool shouldStartTuningPhaseNextIteration(size_t currentIteration, size_t tuningInterval) const override {
    return ((currentIteration + 1) % tuningInterval) == 0;
  };

  void passRuntimeSample([[maybe_unused]] unsigned long sample) override {}

  TuningTriggerOption getOptionType() const override { return TuningTriggerOption::staticSimple; }
};
}  // namespace autopas
