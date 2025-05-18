/**
 * @file AutoTunerInfo.h
 * @author F. Gratl
 * @date 27.06.23
 */

#pragma once

#include "autopas/options/EnergySensorOption.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TuningMetricOption.h"
namespace autopas {

/**
 * Helper struct encapsulating most minor information for the auto tuner.
 */
struct AutoTunerInfo {
  /**
   * Strategy how to select the optimum from the collected samples.
   */
  SelectorStrategyOption selectorStrategy{SelectorStrategyOption::fastestAbs};
  /**
   * Metric used to rate configurations (time or energy).
   */
  TuningMetricOption tuningMetric{TuningMetricOption::time};
  /**
   * Number of time steps after which the auto-tuner shall reevaluate the optimum.
   */
  unsigned int tuningInterval{5000};
  /**
   * Number of samples that shall be collected per combination.
   */
  unsigned int maxSamples{3};

  /**
   * EarlyStoppingFactor for the auto-tuner.
   */
  double earlyStoppingFactor{std::numeric_limits<double>::infinity()};
  /**
   * Used energy sensor of energy metric selected
   */
  EnergySensorOption energySensor{EnergySensorOption::rapl};
  /**
   * Flag for whether LOESS Smoothening is used to smoothen the tuning results.
   */
  bool useLOESSSmoothening{false};
};
}  // namespace autopas