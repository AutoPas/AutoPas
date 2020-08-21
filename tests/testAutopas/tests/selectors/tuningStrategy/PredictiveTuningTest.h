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
   * Tunes for a few iterations (length of times vector) and checks whether all possible configurations were tuned.
   * @param predictiveTuning The PredictiveTuning strategy.
   * @param evidenceList Array of times to add to the evidenceList.
   * @param iteration The current number of iterations, will be increased accordingly.
   * @return Vector of the tested configurations.
   */
  static auto tuneForSomeIterationsAndCheckAllTuned(
      autopas::PredictiveTuning &predictiveTuning, const std::vector<long> &evidenceList, size_t &iteration,
      const std::vector<autopas::Configuration> &allConfigurations = allConfigs) {
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

  static autopas::PredictiveTuning getPredictiveTuning(
      unsigned int testsUntilFirstPrediction, autopas::ExtrapolationMethodOption extrapolationMethodOption,
      const std::set<autopas::TraversalOption> &allowedTraversalOptions = {
          autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced});

  void testGeneric(autopas::ExtrapolationMethodOption extrapolationMethodOption,
                   const std::vector<std::vector<long>> &evidenceVectors, size_t optimalPredictionIndex);

  static constexpr autopas::Configuration configurationLC_C01 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);
  static constexpr autopas::Configuration configurationLC_C08 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  static constexpr autopas::Configuration configurationLC_Sliced = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
      autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  ///@todo c++20: this can be made constexpr, as std::vector will get to be constexpr.
  inline static std::vector<autopas::Configuration> allConfigs{configurationLC_C01, configurationLC_C08,
                                                               configurationLC_Sliced};

  static constexpr double relativeOptimumRange{1.2};
  static constexpr unsigned int maxTuningIterationsWithoutTest{5};
  static constexpr unsigned int relativeRangeForBlacklist{0};
  static constexpr unsigned int evidenceFirstPrediction{2};
};
