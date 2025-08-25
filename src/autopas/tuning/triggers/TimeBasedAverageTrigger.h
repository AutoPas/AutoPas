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
 * is greater or equal to the running average of the last n iterations runtime times triggerFactor.
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
      AutoPasLog(WARN, "triggerFactor for TimeBasedAverageTrigger is {}, but has to be >= 0. Defaulted to 1.25.",
                 triggerFactor);
      _triggerFactor = 1.25;
    }
    if (_nSamples < 1) {
      AutoPasLog(WARN, "nSamples for TimeBasedAverageTrigger is {}, but has to be > 0. Defaulted to 500.", _nSamples);
      _nSamples = 500;
    }

    _runtimeSamples.resize(_nSamples);
    _headElement = 0;
    _sampleSum = 0;
    _bufferFull = false;
  };

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) override {
    bool oldTriggerState = _wasTriggered;
    _wasTriggered = false;

    if (!_bufferFull) [[unlikely]]
      return oldTriggerState;

    unsigned long average = _sampleSum / _nSamples;
    _wasTriggered = (_currentIterationRuntime >= (_triggerFactor * average));
    return oldTriggerState;
  }

  inline bool needsRuntimeSample() const override { return true; }

  void passRuntimeSample(unsigned long sample) override {
    // do not compare to stale samples from before tuning phase
    if (_wasTriggered) [[unlikely]] {
      _runtimeSamples.clear();
      _runtimeSamples.resize(_nSamples);
      _headElement = 0;
      _bufferFull = false;
      return;
    }

    // we use a simple circular buffer to store the last _nSamples samples
    if (_bufferFull) {
      _sampleSum -= _runtimeSamples[_headElement];
    }
    _runtimeSamples[_headElement] = _currentIterationRuntime;
    _sampleSum += _currentIterationRuntime;

    _headElement = (_headElement + 1) % _nSamples;
    _bufferFull = _bufferFull || (_headElement == 0);
    _currentIterationRuntime = sample;
  }

  TuningTriggerOption getOptionType() const override { return TuningTriggerOption::timeBasedAverage; }

 private:
  float _triggerFactor;
  unsigned _nSamples;
  unsigned long _currentIterationRuntime;
  std::vector<unsigned long> _runtimeSamples;
  size_t _headElement;
  unsigned long _sampleSum;
  bool _bufferFull;
};
}  // namespace autopas
