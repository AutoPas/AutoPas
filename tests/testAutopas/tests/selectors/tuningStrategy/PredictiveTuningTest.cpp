/**
 * @file PredictiveTuningTest.cpp
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#include "PredictiveTuningTest.h"

#include <gmock/gmock-more-matchers.h>

autopas::PredictiveTuning PredictiveTuningTest::getPredictiveTuning(
    unsigned int testsUntilFirstPrediction, autopas::ExtrapolationMethodOption extrapolationMethodOption,
    double blacklistRange, const std::set<autopas::TraversalOption> &allowedLCTraversalOptions) {
  return autopas::PredictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.}, allowedLCTraversalOptions, {autopas::LoadEstimatorOption::none},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, _relativeOptimumRange,
      _maxTuningIterationsWithoutTest, blacklistRange, testsUntilFirstPrediction, extrapolationMethodOption);
}

void PredictiveTuningTest::simulateTuningPhasesAndCheckPrediction(
    autopas::ExtrapolationMethodOption extrapolationMethodOption,
    const std::vector<std::map<autopas::Configuration, long>> &evidencePerPhase, unsigned int tuningInterval,
    const std::map<autopas::Configuration, long> &expectedPredictions) {
  // setup sanity check
  ASSERT_EQ(evidencePerPhase.back().size(), expectedPredictions.size())
      << "Provided number of evidence and expectation maps do not match! "
         "The test needs one prediction for each configuration that gets evidences until the last phase.";

  size_t iteration = 0;

  auto predictiveTuning = getPredictiveTuning(evidencePerPhase.size(), extrapolationMethodOption);

  // First reset tuning.
  predictiveTuning.reset(iteration);

  // simulate multiple tuning phases
  for (const auto &evidence : evidencePerPhase) {
    // the tuning phase
    simulateTuningPhase(predictiveTuning, evidence, iteration);

    // during these first phases it is expected that no predictions are made
    EXPECT_THAT(predictiveTuning.getConfigurationPredictions(), ::testing::IsEmpty());
  }

  iteration += tuningInterval;

  // fake start of another tuning phase which triggers calculation of predictions.
  predictiveTuning.reset(iteration);

  // check predictions
  const auto &actualPredictions = predictiveTuning.getConfigurationPredictions();
  EXPECT_THAT(actualPredictions, ::testing::UnorderedElementsAreArray(expectedPredictions))
      << "Predictions do not match expectations!";
}

/**
 * Tests the prediction and the selection of the right optimum.
 * Three different configurations:
 *      - first that starts optimal but becomes slowly more expensive.
 *      - second starts expensive but becomes quickly better. (expected optimum)
 *      - third is constant.
 * @note All expected values in this test were calculated by hand and are precise if not stated otherwise.
 */
TEST_P(PredictiveTuningTest, testPredictions) {
  auto extrapolationOption = GetParam();
  // we do the generation + switch case here to ensure a test is generated for every method.
  switch (extrapolationOption) {
    case autopas::ExtrapolationMethodOption::linePrediction: {
      simulateTuningPhasesAndCheckPrediction(
          extrapolationOption,
          {
              {{_configurationLC_C01, 94}, {_configurationLC_C08, 109}, {_configurationLC_Sliced, 101}},
              {{_configurationLC_C01, 97}, {_configurationLC_C08, 103}, {_configurationLC_Sliced, 101}},
          },
          // all predictions are evaluated for the seventh iteration (iteration==6)
          0, {{_configurationLC_C01, 100}, {_configurationLC_C08, 99}, {_configurationLC_Sliced, 101}});
      break;
    }
    case autopas::ExtrapolationMethodOption::lagrange: {
      simulateTuningPhasesAndCheckPrediction(
          extrapolationOption,
          {
              {{_configurationLC_C01, 79}, {_configurationLC_C08, 122}, {_configurationLC_Sliced, 101}},
              {{_configurationLC_C01, 90}, {_configurationLC_C08, 111}, {_configurationLC_Sliced, 101}},
              {{_configurationLC_C01, 97}, {_configurationLC_C08, 103}, {_configurationLC_Sliced, 101}},
          },
          // all predictions are evaluated for the tenth iteration (iteration==9)
          // for _configurationLC_C08 we actually expect 99.2, however due to internal rounding errors we end up with 98
          0, {{_configurationLC_C01, 100}, {_configurationLC_C08, 98}, {_configurationLC_Sliced, 101}});
      break;
    }
    case autopas::ExtrapolationMethodOption::newton: {
      simulateTuningPhasesAndCheckPrediction(
          extrapolationOption,
          {
              {{_configurationLC_C01, 79}, {_configurationLC_C08, 115}, {_configurationLC_Sliced, 101}},
              {{_configurationLC_C01, 90}, {_configurationLC_C08, 109}, {_configurationLC_Sliced, 101}},
              {{_configurationLC_C01, 97}, {_configurationLC_C08, 103}, {_configurationLC_Sliced, 101}},
          },
          // all predictions are evaluated for the tenth iteration (iteration==9)
          // for _configurationLC_C08 we actually expect 99.33
          0, {{_configurationLC_C01, 100}, {_configurationLC_C08, 99}, {_configurationLC_Sliced, 101}});
      break;
    }
    case autopas::ExtrapolationMethodOption::linearRegression: {
      simulateTuningPhasesAndCheckPrediction(
          extrapolationOption,
          {
              // values along the desired line offset so that the regression should still return the same line
              {{_configurationLC_C01, 91 + 1}, {_configurationLC_C08, 115 - 1}, {_configurationLC_Sliced, 101 - 2}},
              {{_configurationLC_C01, 94 - 2}, {_configurationLC_C08, 109 + 2}, {_configurationLC_Sliced, 101 + 4}},
              {{_configurationLC_C01, 97 + 1}, {_configurationLC_C08, 103 - 1}, {_configurationLC_Sliced, 101 - 2}},
          },
          // all predictions are evaluated for the tenth iteration (iteration==9)
          0, {{_configurationLC_C01, 100}, {_configurationLC_C08, 99}, {_configurationLC_Sliced, 101}});
      break;
    }
    default:
      GTEST_FAIL() << "No test implemented for extrapolation method " << extrapolationOption.to_string();
  }
}

/**
 * Tests the prediction and the selection of the right optimum. Values are chosen so that the predictions result in
 * under or overflow.
 *
 * @note We expect underflows to be set to 1 and overflows to max size_t - 1.
 *
 * (max-value itself is used as a indicator that no prediction was computed).
 */
TEST_P(PredictiveTuningTest, testUnderAndOverflow) {
  auto extrapolationOption = GetParam();
  constexpr auto maxSizeT = std::numeric_limits<size_t>::max();
  simulateTuningPhasesAndCheckPrediction(
      extrapolationOption,
      {
          {{_configurationLC_C01, 100}, {_configurationLC_C08, 1}, {_configurationLC_Sliced, 101}},
          {{_configurationLC_C01, 10}, {_configurationLC_C08, maxSizeT / 4}, {_configurationLC_Sliced, 101}},
      },
      // all predictions are evaluated for the seventh iteration (iteration==6)
      100, {{_configurationLC_C01, 1}, {_configurationLC_C08, maxSizeT - 1}, {_configurationLC_Sliced, 101}});
}

INSTANTIATE_TEST_SUITE_P(Generated, PredictiveTuningTest,
                         ::testing::ValuesIn(autopas::ExtrapolationMethodOption::getAllOptions()),
                         PredictiveTuningTest::PrintToStringParamName());

/**
 * Tests the selection of the right optimumSearchSpace in the regard of the relativeOptimum and tuning.
 * Three different configurations:
 *      - first is constant near the optimum (11).
 *      - second is constant the optimum (10).
 *      - third is constant and not in the optimum range (20).
 * In the third iteration the first and second configuration should be in _optimalSearchSpace and after the tuning
 phase
 * the second should be the optimal configuration.
 */
TEST_F(PredictiveTuningTest, testLinearPredictionTuningThreeIterations) {
  size_t iteration = 0;
  auto predictiveTuning =
      getPredictiveTuning(_evidenceFirstPrediction, autopas::ExtrapolationMethodOption::linePrediction);

  predictiveTuning.reset(iteration);

  std::map<autopas::Configuration, long> evidence{
      {_configurationLC_C01, 11}, {_configurationLC_C08, 10}, {_configurationLC_Sliced, 20}};
  std::vector<std::pair<autopas::Configuration, long>> evidenceSorted(evidence.begin(), evidence.end());
  std::sort(evidenceSorted.begin(), evidenceSorted.end(),
            [](const auto &pairA, const auto &pairB) { return pairA.second < pairB.second; });
  auto optimalConfiguration = evidenceSorted[0].first;
  auto nearOptimalConfiguration = evidenceSorted[1].first;

  simulateTuningPhase(predictiveTuning, evidence, iteration);

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  simulateTuningPhase(predictiveTuning, evidence, iteration);

  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  // This test if a configuration near the optimum gets selected into _optimalSearchSpace.
  EXPECT_EQ(nearOptimalConfiguration, predictiveTuning.getCurrentConfiguration());

  predictiveTuning.addEvidence(11, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);

  // This tests that the right optimum configuration gets selected at the end of a tuning phase.
  predictiveTuning.tune();
  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
}

/**
 * Tests the selection of the right optimumSearchSpace in the regard of the number of iterations without a test and
 * tuning. Two different configurations:
 *      - first is constant out of the optimum range (20).
 *      - second is constant the optimum (10).
 * In iteration three to six only the second configuration should be in _optimalSearchSpace.
 * In the seventh iteration both configuration should be tested in the tuning phase.
 */
TEST_F(PredictiveTuningTest, testLinearPredictionTooLongNotTested) {
  size_t iteration = 0;
  std::vector<autopas::Configuration> configurationsToCompare{_configurationLC_C08, _configurationLC_Sliced};
  auto predictiveTuning =
      getPredictiveTuning(_evidenceFirstPrediction, autopas::ExtrapolationMethodOption::linePrediction, 0,
                          {autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_c08});

  predictiveTuning.reset(iteration);

  std::map<autopas::Configuration, long> evidence{{_configurationLC_C01, 20}, {_configurationLC_C08, 10}};

  std::vector<std::pair<autopas::Configuration, long>> evidenceSorted(evidence.begin(), evidence.end());
  std::sort(evidenceSorted.begin(), evidenceSorted.end(),
            [](const auto &pairA, const auto &pairB) { return pairA.second < pairB.second; });
  auto optimalConfiguration = evidenceSorted[0].first;
  auto badConfiguration = evidenceSorted[1].first;

  simulateTuningPhase(predictiveTuning, evidence, iteration);

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(badConfiguration, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;
  predictiveTuning.tune();

  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;
  predictiveTuning.tune();

  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());

  // Iterating through the maxNumberOfIterationsWithoutTest iterations where C08 should not be tested.
  for (int i = 0; i < _maxTuningIterationsWithoutTest; i++) {
    // End of the (i + 2)-th tuning phase.
    predictiveTuning.reset(iteration);

    EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(10, iteration);
    ++iteration;

    EXPECT_FALSE(predictiveTuning.tune());
    EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
  }

  // End of the (_maxTuningIterationsWithoutTest + 2)-th tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  // Tests that a configuration gets selected into _tooLongNotTested.
  predictiveTuning.tune();
  EXPECT_EQ(badConfiguration, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);

  predictiveTuning.tune();
  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
}

/**
 * Tests the recognition of an invalid _optimalSearchSpace in tuning and the reselection of _optimalSearchSpace
 * Three different configurations:
 *      - first is constant out of the optimum range (15) and invalid in the third iteration.
 *      - second is constant the optimum (10) and invalid in the third iteration.
 *      - third is constant out of the optimum range (20).
 * In the third iteration reselectOptimalSearchSpace should be called twice and the third configuration should be
 * selected after the tuning phase.
 */
TEST_F(PredictiveTuningTest, testInvalidOptimalSearchSpaceTwice) {
  size_t iteration = 0;
  auto predictiveTuning =
      getPredictiveTuning(_evidenceFirstPrediction, autopas::ExtrapolationMethodOption::linePrediction);

  predictiveTuning.reset(iteration);

  std::map<autopas::Configuration, long> evidence{
      {_configurationLC_C01, 15}, {_configurationLC_C08, 10}, {_configurationLC_Sliced, 20}};

  std::vector<std::pair<autopas::Configuration, long>> evidenceSorted(evidence.begin(), evidence.end());
  std::sort(evidenceSorted.begin(), evidenceSorted.end(),
            [](const auto &pairA, const auto &pairB) { return pairA.second < pairB.second; });

  simulateTuningPhase(predictiveTuning, evidence, iteration);

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  simulateTuningPhase(predictiveTuning, evidence, iteration);
  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(evidenceSorted[0].first, predictiveTuning.getCurrentConfiguration());

  // mark best as invalid
  predictiveTuning.tune(true);
  ++iteration;

  // Tests if secondBestConfiguration gets selected if the first one is invalid.
  EXPECT_EQ(evidenceSorted[1].first, predictiveTuning.getCurrentConfiguration());

  // mark second best as invalid
  predictiveTuning.tune(true);
  ++iteration;

  // Tests if thirdBestConfiguration gets selected if the second one is invalid.
  EXPECT_EQ(evidenceSorted[2].first, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);

  // Tests if thirdBestConfiguration still gets selected if the second one is invalid.
  predictiveTuning.tune();
  EXPECT_EQ(evidenceSorted[2].first, predictiveTuning.getCurrentConfiguration());
}

/**
 * Tests the correct use of the blacklist:
 *      - first is constant the optimum (1).
 *      - second is constant in the blacklist range (5).
 *      - third is  out of the optimum range (20) and should not be tested after the first tuning phase.
 * In the second iteration the third configuration should not be tested in the tuning phase.
 */
TEST_F(PredictiveTuningTest, testBlacklist) {
  size_t iteration = 0;
  auto predictiveTuning =
      getPredictiveTuning(_evidenceFirstPrediction, autopas::ExtrapolationMethodOption::linePrediction, 10);

  std::map<autopas::Configuration, long> evidence{
      {_configurationLC_C01, 1}, {_configurationLC_C08, 5}, {_configurationLC_Sliced, 20}};

  std::vector<std::pair<autopas::Configuration, long>> evidenceSorted(evidence.begin(), evidence.end());
  std::sort(evidenceSorted.begin(), evidenceSorted.end(),
            [](const auto &pairA, const auto &pairB) { return pairA.second < pairB.second; });
  auto bestConfiguration = evidenceSorted[0].first;
  auto secondBestConfiguration = evidenceSorted[1].first;

  simulateTuningPhase(predictiveTuning, evidence, iteration);

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  // remove blacklisted evidence
  evidence.erase(_configurationLC_Sliced);

  simulateTuningPhase(predictiveTuning, evidence, iteration);
  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(bestConfiguration, predictiveTuning.getCurrentConfiguration());
}
