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
      AutoPasLog(WARN, "triggerFactor for TimeBasedAverageTrigger is {}, but has to be >= 0. Defaulted to 1.5.\n",
                 triggerFactor);
      _triggerFactor = 1.5;
    }
    if (_nSamples < 1) {
      AutoPasLog(WARN, "nSamples for TimeBasedAverageTrigger is {}, but has to be > 0. Defaulted to 1000.\n",
                 _nSamples);
      _nSamples = 1000;
    }

    _runtimeSamples.resize(_nSamples);
    _headElement = 0;
    _sampleSum = 0;
    _bufferFull = false;
  };

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) override {
    bool _wasTriggered =
        (currentIteration > _nextTriggeringIteration - 10) && (currentIteration <= _nextTriggeringIteration);

    if (!_wasTriggered) {
      if (!_bufferFull) [[unlikely]]
        return false;

      unsigned long average = _sampleSum / _nSamples;
      if (_currentIterationRuntime >= (_triggerFactor * average)) _nextTriggeringIteration = currentIteration + 10;
      return false;
    }

    if (currentIteration == _nextTriggeringIteration) {
      // Do not compare to stale samples from before tuning phase.
      _headElement = 0;
      _sampleSum = 0;
      _bufferFull = false;
      return true;
    }
    return false;
  }

  void passRuntimeSample(unsigned long sample) override {
    // Use a simple circular buffer to store the last _nSamples samples.
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
  unsigned long _headElement;
  unsigned long _sampleSum;
  bool _bufferFull;
};
}  // namespace autopas
