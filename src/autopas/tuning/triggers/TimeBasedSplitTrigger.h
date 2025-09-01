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
 * Represents a time-based trigger. This trigger considers the last n runtime samples together with the current
 * iteration runtime and divides them up into twointervals, A and B. A new tuning phase is triggered, if avg(B) is
 * greater or equal to the avg(A) times triggerFactor.
 *
 * The detailed construction of A, B is as follows:  A := [t_{i−n}, t_{i−j}], B := [t_{i−j+1}, t_i], where j =
 * ceil(n/2).
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
    // Interval B contains the current iteration, therefore its length is increased by 1.
    _intervalLengthB = nSamples / 2 + 1;
    _runtimeSamplesA.resize(_intervalLengthA);
    _runtimeSamplesB.resize(_intervalLengthB);
    _headElementA = _headElementB = 0;
    _sampleSumA = _sampleSumB = 0;
    _bufferFullA = _bufferFullB = false;
  };

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) override {
    bool oldTriggerState = _wasTriggered;
    _wasTriggered = false;

    if (!_bufferFullB) return oldTriggerState;

    unsigned long averageA = _sampleSumA / _intervalLengthA;
    unsigned long averageB = _sampleSumB / _intervalLengthB;

    _wasTriggered = (averageB >= (_triggerFactor * averageA));
    return oldTriggerState;
  }

  inline bool needsRuntimeSample() const override { return true; }

  void passRuntimeSample(unsigned long sample) override {
    // Do not compare to stale samples from before tuning phase.
    if (_wasTriggered) [[unlikely]] {
      _headElementA = _headElementB = 0;
      _sampleSumA = _sampleSumB = 0;
      _bufferFullA = _bufferFullB = false;
    }

    // Use two circular buffers to store the last (_nSamples+1) samples.
    if (!_bufferFullA) {
      _runtimeSamplesA[_headElementA] = sample;
      _sampleSumA += sample;
      _headElementA = (_headElementA + 1) % _intervalLengthA;

      _bufferFullA = _headElementA == 0;
    } else {
      if (_bufferFullB) {
        _sampleSumA -= _runtimeSamplesA[_headElementA];

        // Handles the sample passing from interval B to A.
        unsigned long passing_sample = _runtimeSamplesB[_headElementB];
        _runtimeSamplesA[_headElementA] = passing_sample;
        _sampleSumA += passing_sample;
        _headElementA = (_headElementA + 1) % _intervalLengthA;

        _sampleSumB -= passing_sample;
      }

      _runtimeSamplesB[_headElementB] = sample;
      _sampleSumB += sample;
      _headElementB = (_headElementB + 1) % _intervalLengthB;

      _bufferFullB = _bufferFullB || (_headElementB == 0);
    }
  }

  TuningTriggerOption getOptionType() const override { return TuningTriggerOption::timeBasedSplit; }

 private:
  float _triggerFactor;
  unsigned _intervalLengthA;
  unsigned _intervalLengthB;
  std::vector<unsigned long> _runtimeSamplesA;
  std::vector<unsigned long> _runtimeSamplesB;
  unsigned long _headElementA;
  unsigned long _headElementB;
  unsigned long _sampleSumA;
  unsigned long _sampleSumB;
  bool _bufferFullA;
  bool _bufferFullB;
};
}  // namespace autopas
