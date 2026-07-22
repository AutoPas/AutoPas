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

  /**
   * Default number of particles in two cells from which AoS sorting should be performed, used to seed this tuner's
   * SortingThresholdBenchmark before/without a benchmark run. Mirrors LogicHandlerInfo::aosSortingThreshold.
   * For details on the chosen default threshold see: https://github.com/AutoPas/AutoPas/pull/619
   */
  size_t aosSortingThreshold{8};

  /**
   * Default number of particles in two SoA buffers from which SoA sorting should be performed, used to seed this
   * tuner's SortingThresholdBenchmark before/without a benchmark run. Mirrors LogicHandlerInfo::soaSortingThreshold.
   * Default comes from the LJFunctorHWY Benchmarks.
   */
  size_t soaSortingThreshold{100};
};
}  // namespace autopas