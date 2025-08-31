/**
 * @file TimeBasedRegressionTrigger.h
 * @author Niklas Ladurner
 * @date 13.07.2025
 */

#pragma once

#include <vector>

#include "autopas/tuning/triggers/TuningTriggerInterface.h"

namespace autopas {

/**
 * Represents a time-based trigger. This trigger considers the last n runtime samples and performs a simple linear
 * regression on them to get an estimate of the normalized slope of runtime increase. A new tuning phase is triggered,
 * if this slope is greater or equal to triggerFactor. This means that a slope > 1.0 can be interpreted as a projected
 * increase in iteration runtime, and a slope < 1.0 as projected decrease in iteration runtime.
 *
 * The detailed method is as follows:
 * In the following, $n$ is the number of samples, $t_k$ the runtime at iteration $k$, $i$ the current iteration and
 * $\bar t$ the average runtime.
 * As we only care for the slope factor estimate \hat\beta_1, we can transform the standard simple linear regression to
 * \hat\beta_1' = \frac{1}{C_2}\sum_{k=0}^{n}(k-C_1)(t_{i-n+k}-\bar t)
 * Where C_1=\frac{n}{2}, C_2=\sum_{k=0}^{n}(k-C_1)^2=\frac{n(n+1)(n+2)}{12} can be precomputed at initialization.
 *
 * Then, we normalize $\hat\beta_1'$ to make it independent of the scenario:
 * \hat\beta_{\text{norm}} = 1+\frac{(n+1)\hat\beta_1'}{2\bar t}
 */
class TimeBasedRegressionTrigger : public TuningTriggerInterface {
 public:
  /**
   * Constructor
   *
   * @param triggerFactor The positive threshold value at which, if exceeded by the normalized slope of the regression,
   * a new tuning phase is triggered.
   * @param nSamples The number of runtime samples to consider.
   */
  TimeBasedRegressionTrigger(float triggerFactor, unsigned nSamples)
      : _triggerFactor(triggerFactor), _nSamples(nSamples) {
    if (triggerFactor <= 0) {
      AutoPasLog(WARN, "triggerFactor for TimeBasedRegressionTrigger is {}, but has to be > 0. Defaulted to 1.5.",
                 triggerFactor);
      _triggerFactor = 1.5;
    }
    if (_nSamples < 1) {
      AutoPasLog(WARN, "nSamples for TimeBasedRegressionTrigger is {}, but has to be > 0. Defaulted to 1000.",
                 _nSamples);
      _nSamples = 1000;
    }

    _C1 = _nSamples / 2.0;
    _C2_reciprocal = 12.0 / (_nSamples * (_nSamples + 1) * (_nSamples + 2));

    // The vector also contains the current iteration, therefore its length is _nSamples +1.
    _runtimeSamples.resize(_nSamples + 1);
    _headElement = 0;
    _sampleSum = 0;
    _bufferFull = false;
  };

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) override {
    bool oldTriggerState = _wasTriggered;
    _wasTriggered = false;

    if (!_bufferFull) [[unlikely]]
      return oldTriggerState;

    double t_avg = _sampleSum / (_nSamples + 1.0);

    double beta = 0.0;
    for (size_t k = 0; k < (_nSamples + 1); k++) {
      beta += (k - _C1) * (_runtimeSamples.at((_headElement + k) % (_nSamples + 1)) - t_avg);
    }

    beta *= _C2_reciprocal;
    double beta_normalized = 1 + (0.5 * (_nSamples + 1) * beta) / t_avg;

    _wasTriggered = beta_normalized >= _triggerFactor;
    return oldTriggerState;
  }

  inline bool needsRuntimeSample() const override { return true; }

  void passRuntimeSample(unsigned long sample) override {
    // Do not compare to stale samples from before tuning phase.
    if (_wasTriggered) [[unlikely]] {
      _headElement = 0;
      _sampleSum = 0;
      _bufferFull = false;
    }

    // Use a simple circular buffer to store the last (_nSamples+1) samples.
    if (_bufferFull) {
      _sampleSum -= _runtimeSamples[_headElement];
    }
    _runtimeSamples[_headElement] = sample;
    _sampleSum += sample;

    _headElement = (_headElement + 1) % (_nSamples + 1);
    _bufferFull = _bufferFull || (_headElement == 0);
  }

  TuningTriggerOption getOptionType() const override { return TuningTriggerOption::timeBasedRegression; }

 private:
  float _triggerFactor;
  unsigned _nSamples;
  double _C1;
  double _C2_reciprocal;
  std::vector<unsigned long> _runtimeSamples;
  size_t _headElement;
  unsigned long _sampleSum;
  bool _bufferFull;
};
}  // namespace autopas
