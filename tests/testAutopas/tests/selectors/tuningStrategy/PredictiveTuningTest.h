/**
 * @file PredictiveTuningTest.h
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#pragma once

#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/options/ExtrapolationMethodOption.h"
#include "autopas/selectors/tuningStrategy/PredictiveTuning.h"

class PredictiveTuningTest : public AutoPasTestBase {
 protected:
  /**
   * Tunes for a few iterations (length of evidenceList) and checks whether all possible configurations were tuned.
   * @param predictiveTuning The PredictiveTuning strategy.
   * @param evidenceList Array of times to add to the evidenceList.
   * @param iteration The current number of iterations, will be increased accordingly.
   * @return Vector of the tested configurations.
   */
  static auto tuneForSomeIterationsAndCheckAllTuned(
      autopas::PredictiveTuning &predictiveTuning, const std::vector<long> &evidenceList, size_t &iteration,
      const std::vector<autopas::Configuration> &allConfigurations = _allConfigs) {
    std::vector<autopas::Configuration> testedConfigs(evidenceList.size());
    autopas::Configuration optimalConfiguration;
    auto minTime = std::numeric_limits<size_t>::max();
    for (size_t index = 0ul; index < evidenceList.size(); ++index) {
      testedConfigs[index] = predictiveTuning.getCurrentConfiguration();
      if (evidenceList[index] < minTime) {
        optimalConfiguration = predictiveTuning.getCurrentConfiguration();
        minTime = evidenceList[index];
      }
      predictiveTuning.addEvidence(evidenceList[index], iteration);
      ++iteration;
      predictiveTuning.tune();
    }
    EXPECT_THAT(allConfigurations, testing::UnorderedElementsAreArray(testedConfigs));

    EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());

    return testedConfigs;
  }

  /**
   * Creates an instance of a Predictive Tuner with the given parameters.
   * The search space is fixed to LC, no N3, SoA, CSF = 1, LoadEstimator = none, and all given traversal options.
   * @param testsUntilFirstPrediction
   * @param extrapolationMethodOption
   * @param blacklistRange
   * @param allowedLCTraversalOptions This is basically the whole search space.
   * @return Preconfigured PredictiveTuning object.
   */
  static autopas::PredictiveTuning getPredictiveTuning(
      unsigned int testsUntilFirstPrediction, autopas::ExtrapolationMethodOption extrapolationMethodOption,
      double blacklistRange = 0,
      const std::set<autopas::TraversalOption> &allowedLCTraversalOptions = {
          autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced});

  /**
   * Simulates multiple tuning phases. The number of phases depends on the length of the provided vector of evidence.
   *
   * A predictive tuning object is generated that contains a fixed set of configurations (see _allConfigs).
   * The provided evidence is fed to the tuning strategy as if it was tuning for several phases with the assumption
   * that the Nth entry in the sub-vectors refer to the same configurations.
   * In the end it is checked whether the tuner decides on the configuration behind the optimalPredictionIndex.
   *
   * @param extrapolationMethodOption
   * @param evidenceVectors 2D vector containing all evidence that shall be added per tuning phase.
   * @param optimalPredictionIndex Index referring to the second dimension of the evidenceVectors that
   * should produce the optimal prediction.
   */
  void simulateTuningPhases(autopas::ExtrapolationMethodOption extrapolationMethodOption,
                            const std::vector<std::vector<long>> &evidenceVectors, size_t optimalPredictionIndex);

  static constexpr autopas::Configuration _configurationLC_C01 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);
  static constexpr autopas::Configuration _configurationLC_C08 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  static constexpr autopas::Configuration _configurationLC_Sliced = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
      autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  ///@todo c++20: this can be made constexpr, as std::vector will get to be constexpr.
  inline static std::vector<autopas::Configuration> _allConfigs{_configurationLC_C01, _configurationLC_C08,
                                                                _configurationLC_Sliced};

  static constexpr double _relativeOptimumRange{1.2};
  static constexpr unsigned int _maxTuningIterationsWithoutTest{5};
  static constexpr unsigned int _evidenceFirstPrediction{2};
};
