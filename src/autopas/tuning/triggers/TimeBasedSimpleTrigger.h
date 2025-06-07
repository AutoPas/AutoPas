/**
 * @file TimeBasedSimpleTrigger.h
 * @author Niklas Ladurner
 * @date 06.06.2025
 */

#pragma once

#include "autopas/tuning/triggers/TuningTriggerInterface.h"

namespace autopas {

/**
 * Represents a naive time-based trigger. A new tuning phase is triggered, if the runtime of the current iteration
 * is greater or equal to the last iteration's runtime times triggerFactor.
 */
class TimeBasedSimpleTrigger : public TuningTriggerInterface {
 public:
  /**
   * Constructor
   *
   * @param triggerFactor The factor by which the current iteration runtime must exceed the last iteration runtime to
   * trigger a new tuning phase.
   */
  TimeBasedSimpleTrigger(float triggerFactor) : _triggerFactor(triggerFactor) {};

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) const override {
    return _currentIterationRuntime >= (_triggerFactor * _lastIterationRuntime);
  }

  inline bool needsRuntimeSample() const override { return true; }

  void passRuntimeSample(unsigned long sample) override {
    _lastIterationRuntime = _currentIterationRuntime;
    _currentIterationRuntime = sample;
  }

  TuningTriggerOption getOptionType() const override { return TuningTriggerOption::timeBasedSimple; }

 private:
  float _triggerFactor;
  unsigned long _currentIterationRuntime;
  unsigned long _lastIterationRuntime;
};
}  // namespace autopas
