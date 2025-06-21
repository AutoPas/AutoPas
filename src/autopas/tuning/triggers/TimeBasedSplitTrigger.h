/**
 * @file TimeBasedSplitTrigger.h
 * @author Niklas Ladurner
 * @date 17.06.2025
 */

#pragma once

#include <vector>

#include "autopas/tuning/triggers/TuningTriggerInterface.h"

namespace autopas {

/**
 * Represents a time-based trigger. This trigger considers the last n runtime samples and divides them up into two
 * intervals, A and B. A new tuning phase is triggered, if avg(B) is greater or equal to the avg(A) times triggerFactor.
 *
 * The detailed construction of A, B is as follows:  A := [t_{i−n}, t_{i−j}], B := [t_{i−j+1}, t_i], where j =
 * floor(n/2).
 */
class TimeBasedSplitTrigger : public TuningTriggerInterface {
 public:
  /**
   * Constructor
   *
   * @param triggerFactor The factor by which the current iteration runtime must exceed the last iteration runtime to
   * trigger a new tuning phase.
   * @param nSamples The number of runtime samples to consider.
   */
  TimeBasedSplitTrigger(float triggerFactor, unsigned nSamples) : _triggerFactor(triggerFactor) {
    if (triggerFactor < 0) {
      AutoPasLog(WARN, "triggerFactor for TimeBasedSplitTrigger is {}, but has to be >= 0. Defaulted to 1.",
                 triggerFactor);
      _triggerFactor = 1.0;
    }
    if (nSamples < 1) {
      AutoPasLog(WARN, "nSamples for TimeBasedSplitTrigger is {}, but has to be > 0. Defaulted to 10.", nSamples);
      nSamples = 10;
    }

    _intervalLengthA = (nSamples + 1) / 2;
    _intervalLengthB = nSamples / 2;
    _runtimeSamples.reserve(_intervalLengthA + _intervalLengthB);
  };

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) override {
    bool oldTriggerState = _wasTriggered;
    _wasTriggered = false;

    if (_runtimeSamples.size() < _intervalLengthA + _intervalLengthB) return oldTriggerState;

    unsigned long averageA =
        std::reduce(_runtimeSamples.begin(), _runtimeSamples.begin() + _intervalLengthA) / _intervalLengthA;
    unsigned long averageB =
        std::reduce(_runtimeSamples.begin() + _intervalLengthA, _runtimeSamples.end()) / _intervalLengthB;

    _wasTriggered = (averageB >= (_triggerFactor * averageA));
    return oldTriggerState;
  }

  inline bool needsRuntimeSample() const override { return true; }

  void passRuntimeSample(unsigned long sample) override {
    if (_runtimeSamples.size() == _intervalLengthA + _intervalLengthB) _runtimeSamples.erase(_runtimeSamples.begin());
    _runtimeSamples.push_back(sample);

    // if we start a new tuning phase, clear runtime samples such that we do not compare runtimes between different
    // configurations
    if (_wasTriggered) _runtimeSamples.clear();
  }

  TuningTriggerOption getOptionType() const override { return TuningTriggerOption::timeBasedSplit; }

 private:
  float _triggerFactor;
  unsigned _intervalLengthA;
  unsigned _intervalLengthB;
  std::vector<unsigned long> _runtimeSamples;
};
}  // namespace autopas
