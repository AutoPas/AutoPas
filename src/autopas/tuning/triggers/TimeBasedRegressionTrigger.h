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
 * \hat\beta_1' = \frac{1}{B}\sum_{k=0}^{n-1}(k-A)(t_{i-n-1+k}-\bar t)
 * Where A = \frac{n-1}{2}, B=\sum_{k=0}^{n-1}(k-A)^2 can be precomputed at initialization.
 *
 * Then, we normalize $\hat\beta_1'$ to make it independent of the scenario by dividing it by the average iteration
 * runtime in the current interval.
 */
class TimeBasedRegressionTrigger : public TuningTriggerInterface {
 public:
  /**
   * Constructor
   *
   * @param triggerFactor The positive threshold value at which, if exceeded by the normalized slope of the regression, a new
   * tuning phase is triggered.
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
      AutoPasLog(WARN, "nSamples for TimeBasedRegressionTrigger is {}, but has to be > 0. Defaulted to 10.", _nSamples);
      _nSamples = 10;
    }

    _A = (_nSamples - 1.0) / 2.0;
    double B = .0;
    for (size_t k = 0; k < _nSamples; k++) B += (k - _A) * (k - _A);
    _B_reciprocal = 1.0 / B;

    _runtimeSamples.reserve(nSamples);
  };

  inline bool shouldStartTuningPhase(size_t currentIteration, size_t tuningInterval) override {
    bool oldTriggerState = _wasTriggered;
    _wasTriggered = false;

    if (_runtimeSamples.size() < _nSamples) return oldTriggerState;

    double t_avg = std::reduce(_runtimeSamples.begin(), _runtimeSamples.end()) / (_nSamples * 1.0);

    double beta = 0.0;
    for (size_t k = 0; k < _nSamples; k++) {
      beta += (k - _A) * (_runtimeSamples.at(k) - t_avg);
    }

    beta *= _B_reciprocal;
    double beta_normalized = 1 + _nSamples * beta / t_avg;

    _wasTriggered = beta_normalized >= _triggerFactor;

    return oldTriggerState;
  }

  inline bool needsRuntimeSample() const override { return true; }

  void passRuntimeSample(unsigned long sample) override {
    if (_runtimeSamples.size() == _nSamples) _runtimeSamples.erase(_runtimeSamples.begin());
    _runtimeSamples.push_back(sample);

    // if we start a new tuning phase, clear runtime samples such that we do not compare runtimes between different
    // configurations
    if (_wasTriggered) _runtimeSamples.clear();
  }

  TuningTriggerOption getOptionType() const override { return TuningTriggerOption::timeBasedRegression; }

 private:
  float _triggerFactor;
  unsigned _nSamples;
  double _A;
  double _B_reciprocal;
  std::vector<unsigned long> _runtimeSamples;
};
}  // namespace autopas
