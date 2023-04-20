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

class PredictiveTuningTest : public AutoPasTestBase,
                             public ::testing::WithParamInterface<autopas::ExtrapolationMethodOption> {
 public:
  struct PrintToStringParamName {
    template <class ParamType>
    std::string operator()(const testing::TestParamInfo<ParamType> &info) const {
      auto extrapolationOption = static_cast<ParamType>(info.param);
      auto str = extrapolationOption.to_string();
      std::replace(str.begin(), str.end(), '-', '_');
      std::replace(str.begin(), str.end(), '.', '_');
      return str;
    }
  };

 protected:
  /**
   * Tunes for a few iterations (length of evidence) and checks whether all possible configurations were tuned.
   * @param predictiveTuning The PredictiveTuning strategy.
   * @param evidence Map of configurations to evidence. These simulate the performance of every configuration.
   * @param expectedPredictions Map of configurations to expected predictions. Pass an empty map if this should be
   * skipped.
   * @param iteration The current number of iterations, will be increased accordingly.
   * @return Vector of the tested configurations.
   */
  static auto simulateTuningPhase(autopas::PredictiveTuning &predictiveTuning,
                                  const std::map<autopas::Configuration, long> &evidence, size_t &iteration) {
    // input sanity checks
    ASSERT_NE(evidence.size(), 0);

    // collects all configurations that are tested in this phase
    std::vector<autopas::Configuration> testedConfigs;
    testedConfigs.reserve(evidence.size());

    // identify optimum of the provided input
    auto [bestConfig, bestEvidence] =
        *(std::min_element(evidence.begin(), evidence.end(),
                           [](const auto &pairA, const auto &pairB) { return pairA.second < pairB.second; }));

    // strategies are reset at the start of every phase.
    predictiveTuning.reset(iteration);

    // simulate the iterations of the tuning phase
    bool stillTuning = true;
    while (stillTuning) {
      // the configuration for this iteration
      auto config = predictiveTuning.getCurrentConfiguration();
      testedConfigs.push_back(config);

      EXPECT_NO_THROW(evidence.at(config)) << "Did not expect to test configuration " << config;

      predictiveTuning.addEvidence(evidence.at(config), iteration);

      stillTuning = predictiveTuning.tune();
      ++iteration;
    }

    // gather vector of all configs that are expected to be tested
    std::vector<autopas::Configuration> allConfigurations;
    allConfigurations.reserve(evidence.size());
    for (const auto &[config, _] : evidence) {
      allConfigurations.push_back(config);
    }
    EXPECT_THAT(testedConfigs, testing::UnorderedElementsAreArray(allConfigurations))
        << "Test did not test all expected configurations in tuning phase.";

    // at the end of a tuning phase the strategy should point to the optimum
    EXPECT_EQ(predictiveTuning.getCurrentConfiguration(), bestConfig)
        << "Did not select the fastest configuration. At the end of the tuning phase.";
  }

  /**
   * Creates an instance of a Predictive Tuner with the given parameters.
   * The search space is fixed to LC, no N3, SoA, CSF = 1, LoadEstimator = none, and all given traversal options.
   * @param testsUntilFirstPrediction
   * @param extrapolationMethodOption
   * @param blacklistRange Set to zero to disable blacklisting.
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
   * The provided evidence is fed to the tuning strategy as if it was tuning for several phases.
   * In the end it is checked whether the tuner decides on the configuration behind the optimalPredictionIndex.
   *
   * @param extrapolationMethodOption
   * @param evidencePerPhase
   * @param tuningInterval Number of iterations between the last tuning iteration and evaluation of expectations.
   * @param expectedPredictions
   */
  void simulateTuningPhasesAndCheckPrediction(
      autopas::ExtrapolationMethodOption extrapolationMethodOption,
      const std::vector<std::map<autopas::Configuration, long>> &evidencePerPhase, unsigned int tuningInterval,
      const std::map<autopas::Configuration, long> &expectedPredictions);

  static constexpr autopas::Configuration _configurationLC_C01 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled, 15);
  static constexpr autopas::Configuration _configurationLC_C08 = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled, 15);

  static constexpr autopas::Configuration _configurationLC_Sliced = autopas::Configuration(
      autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
      autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled, 15);

  ///@todo c++20: this can be made constexpr, as std::vector will get to be constexpr.
  inline static std::vector<autopas::Configuration> _allConfigs{_configurationLC_C01, _configurationLC_C08,
                                                                _configurationLC_Sliced};

  static constexpr double _relativeOptimumRange{1.2};
  static constexpr unsigned int _maxTuningIterationsWithoutTest{5};
  static constexpr unsigned int _evidenceFirstPrediction{2};
};
