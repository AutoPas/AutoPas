/**
 * @file AutoTunerInfo.h
 * @author F. Gratl
 * @date 27.06.23
 */

#pragma once

#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TuningMetricOption.h"
namespace autopas {

/**
 * Helper struct encapsulating most minor information for the auto tuner.
 */
struct AutoTunerInfo {
  /**
   * For MPI-tuning: Maximum of the relative difference in the comparison metric for two ranks which exchange their
   * tuning information.
   */
  double MPITuningMaxDifferenceForBucket{3.0};
  /**
   * For MPI-tuning: Weight for maxDensity in the calculation for bucket distribution.
   */
  double MPITuningWeightForMaxDensity{0.0};
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
   * The verlet rebuild frequency this AutoPas instance uses.
   */
  unsigned int rebuildFrequency{10};
  /**
   * Suffix for all output files produced by this class.
   */
  std::string outputSuffix{};
  /**
   * Whether to use the tuning strategy logger proxy to log tuning information.
   */
  bool useTuningStrategyLoggerProxy{false};
};
}  // namespace autopas