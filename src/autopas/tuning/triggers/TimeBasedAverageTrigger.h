/**
 * @file TimeBasedAverageTrigger.h
 * @author Niklas Ladurner
 * @date 07.06.2025
 */

#pragma once

#include <vector>

#include "autopas/tuning/triggers/TuningTriggerInterface.h"

namespace autopas {

/**
 * Represents a time-based trigger. A new tuning phase is triggered, if the runtime of the current iteration
 * is greater or equal to the running average of the last n iteration's runtime times triggerFactor.
 */
class TimeBasedAverageTrigger : public TuningTriggerInterface {
 public:
  /**
   * Constructor
   *
   * @param triggerFactor The factor by which the current iteration runtime must exceed the last iteration runtime to
   * trigger a new tuning phase.
   * @param nSamples The number of past runtime samples to take the average of.
   */
  TimeBasedAverageTrigger(float triggerFactor, unsigned nSamples) : _triggerFactor(triggerFactor), _nSamples(nSamples) {
    if (triggerFactor < 0) {
      AutoPasLog(WARN, "triggerFactor for TimeBasedAverageTrigger is {}, but has to be >= 0. Defaulted to 1.",
                 triggerFactor);
      _triggerFactor = 1.0;
    }
    if (_nSamples < 1) {
      AutoPasLog(WARN, "nSamples for TimeBasedAverageTrigger is {}, but has to be > 0. Defaulted to 10.", _nSamples);
      _nSamples = 10;
    }
  };

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) override {
    bool oldTriggerState = _wasTriggered;
    _wasTriggered = false;

    if (_lastIterationRuntimes.size() < _nSamples) return oldTriggerState;

    unsigned long average = std::reduce(_lastIterationRuntimes.begin(), _lastIterationRuntimes.end()) / _nSamples;
    _wasTriggered = (_currentIterationRuntime >= (_triggerFactor * average));
    return oldTriggerState;
  }

  inline bool needsRuntimeSample() const override { return true; }

  void passRuntimeSample(unsigned long sample) override {
    if (_lastIterationRuntimes.size() == _nSamples) _lastIterationRuntimes.erase(_lastIterationRuntimes.begin());
    _lastIterationRuntimes.push_back(_currentIterationRuntime);
    _currentIterationRuntime = sample;
  }

  TuningTriggerOption getOptionType() const override { return TuningTriggerOption::timeBasedAverage; }

 private:
  float _triggerFactor;
  unsigned _nSamples;
  unsigned long _currentIterationRuntime;
  std::vector<unsigned long> _lastIterationRuntimes;
};
}  // namespace autopas
