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
      AutoPasLog(WARN, "triggerFactor for TimeBasedSplitTrigger is {}, but has to be >= 0. Defaulted to 1.25.",
                 triggerFactor);
      _triggerFactor = 1.25;
    }
    if (nSamples < 1) {
      AutoPasLog(WARN, "nSamples for TimeBasedSplitTrigger is {}, but has to be > 0. Defaulted to 500.", nSamples);
      nSamples = 500;
    }

    _intervalLengthA = (nSamples + 1) / 2;
    _intervalLengthB = nSamples / 2;
    _nSamples = nSamples;
    _runtimeSamples.resize(_nSamples);
    _headElement = 0;
    _sampleSumA = 0;
    _sampleSumB = 0;
    _bufferFull = false;
  };

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) override {
    bool oldTriggerState = _wasTriggered;
    _wasTriggered = false;

    if (!_bufferFull) return oldTriggerState;

    unsigned long averageA = _sampleSumA / _intervalLengthA;
    unsigned long averageB = _sampleSumB / _intervalLengthB;

    _wasTriggered = (averageB >= (_triggerFactor * averageA));
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

    // use a simple circular buffer to store the last _nSamples samples
    if (!_bufferFull && (_headElement < _intervalLengthA)) {
      _runtimeSamples[_headElement] = sample;
      _sampleSumA += sample;
    } else {
      if (_bufferFull) {
        _sampleSumA -= _runtimeSamples[_headElement];
      }
      _runtimeSamples[_headElement] = sample;
      _sampleSumB += sample;

      // the sample passing from interval B to A
      unsigned long midpoint = _runtimeSamples[(_headElement + _intervalLengthA) % _nSamples];
      _sampleSumA += midpoint;
      _sampleSumB -= midpoint;
    }

    _headElement = (_headElement + 1) % _nSamples;
    _bufferFull = _bufferFull || (_headElement == 0);
  }

  TuningTriggerOption getOptionType() const override { return TuningTriggerOption::timeBasedSplit; }

 private:
  float _triggerFactor;
  unsigned _intervalLengthA;
  unsigned _intervalLengthB;
  unsigned long _nSamples;
  std::vector<unsigned long> _runtimeSamples;
  size_t _headElement;
  unsigned long _sampleSumA;
  unsigned long _sampleSumB;
  bool _bufferFull;
};
}  // namespace autopas
