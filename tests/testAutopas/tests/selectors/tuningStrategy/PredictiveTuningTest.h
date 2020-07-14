/**
 * @file PredictiveTuningTest.h
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/selectors/tuningStrategy/PredictiveTuning.h"

class PredictiveTuningTest : public AutoPasTestBase {
 protected:
  /**
   * Tunes for a few iterations (length of times vector) and checks whether all possible configurations were tuned.
   * @param predictiveTuning The PredictiveTuning strategy.
   * @param evidences Array of times to add to the evidences.
   * @param returnConfigIndex Configuration with this index in the times vector is returned.
   * @param iteration The current number of iterations, will be increased accordingly.
   * @return The configuration with index returnConfigIndex.
   */
  autopas::Configuration tuneForSomeIterationsAndCheckAllTuned(autopas::PredictiveTuning &predictiveTuning,
                                                               const std::array<double, 3> &evidences,
                                                               size_t returnConfigIndex, size_t &iteration);

  void testGeneric(autopas::ExtrapolationMethodOption extrapolationMethodOption,
                   const std::vector<std::array<long, 3>> &evidences, size_t optimalPredictionIndex);

  const autopas::Configuration configurationLC_C01 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);
  const autopas::Configuration configurationLC_C08 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  const autopas::Configuration configurationLC_Sliced = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
      autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  std::vector<autopas::Configuration> allConfigs{configurationLC_C01, configurationLC_C08, configurationLC_Sliced};

  static constexpr double relativeOptimumRange{1.2};
  static constexpr unsigned int maxTuningIterationsWithoutTest{5};
  static constexpr unsigned int evidenceFirstPrediction{2};
};
