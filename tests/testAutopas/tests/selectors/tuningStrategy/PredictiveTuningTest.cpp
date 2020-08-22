/**
 * @file PredictiveTuningTest.cpp
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#include "PredictiveTuningTest.h"

autopas::PredictiveTuning PredictiveTuningTest::getPredictiveTuning(
    unsigned int testsUntilFirstPrediction, autopas::ExtrapolationMethodOption extrapolationMethodOption,
    int blacklistRange, const std::set<autopas::TraversalOption> &allowedTraversalOptions) {
  return autopas::PredictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.}, allowedTraversalOptions, {autopas::LoadEstimatorOption::none},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, relativeOptimumRange,
      maxTuningIterationsWithoutTest, blacklistRange, testsUntilFirstPrediction, extrapolationMethodOption);
}

void PredictiveTuningTest::testGeneric(autopas::ExtrapolationMethodOption extrapolationMethodOption,
                                       const std::vector<std::vector<long>> &evidenceVectors,
                                       size_t optimalPredictionIndex) {
  size_t iteration = 0;

  auto predictiveTuning = getPredictiveTuning(evidenceVectors.size(), extrapolationMethodOption);

  // First reset tuning.
  predictiveTuning.reset(iteration);

  autopas::Configuration optimalPrediction;

  for (size_t index = 0; index < evidenceVectors.size(); ++index) {
    if (index == 0) {
      optimalPrediction = tuneForSomeIterationsAndCheckAllTuned(predictiveTuning, evidenceVectors[index],
                                                                iteration)[optimalPredictionIndex];
    } else {
      tuneForSomeIterationsAndCheckAllTuned(predictiveTuning, evidenceVectors[index], iteration);
    }
    // End of the first tuning phase.
    predictiveTuning.reset(iteration);
  }

  // This tests the actual prediction.
  EXPECT_EQ(optimalPrediction, predictiveTuning.getCurrentConfiguration());
}

/*
 * Tests the prediction and the selection of the right optimum.
 * Three different configurations:
 *      - first that is not optimal but it getting less expensive.
 *      - second is optimal in the beginning but is getting more expensive.
 *      - third is constant and not in the optimum range.
 * In the third iteration the first configuration should be predicted to be the optimum.
 */
TEST_F(PredictiveTuningTest, testLinePrediction) {
  testGeneric(autopas::ExtrapolationMethodOption::linePrediction, {{4, 1, 20}, {3, 2, 20}}, 0);
}

/*
 * Tests the prediction and the selection of the right optimum with the linear regression method.
 * Three different configurations:
 *      - first that is not optimal but it getting less expensive.
 *      - second is optimal in the beginning but is getting more expensive.
 *      - third is constant and not in the optimum range.
 * In the third iteration the first configuration should be predicted to be the optimum.
 */
TEST_F(PredictiveTuningTest, testLinearRegression) {
  testGeneric(autopas::ExtrapolationMethodOption::linearRegression, {{375, 300, 2000}, {350, 325, 2000}}, 0);
}

/*
 * Tests the prediction and the selection of the right optimum with the lagrange interpolation method.
 * Three different configurations:
 *      - first that is not optimal but it getting less expensive.
 *      - second is optimal in the beginning but is getting more expensive.
 *      - third is constant and not in the optimum range.
 * In the third iteration the first configuration should be predicted to be the optimum.
 */
TEST_F(PredictiveTuningTest, testLagrange) {
  testGeneric(autopas::ExtrapolationMethodOption::lagrange, {{6, 1, 20}, {5, 2, 20}, {4, 3, 20}}, 0);
}

/*
 * Tests the prediction and the selection of the right optimum with the newton interpolation method.
 * Three different configurations:
 *      - first that is not optimal but it getting less expensive.
 *      - second is optimal in the beginning but is getting more expensive.
 *      - third is constant and not in the optimum range.
 * In the third iteration the first configuration should be predicted to be the optimum.
 */
TEST_F(PredictiveTuningTest, testNewton) {
  testGeneric(autopas::ExtrapolationMethodOption::newton, {{60, 10, 200}, {50, 20, 200}, {40, 30, 200}}, 0);
}

/*
 * Tests the selection of the right optimumSearchSpace in the regard of the relativeOptimum and tuning.
 * Three different configurations:
 *      - first is constant near the optimum (11).
 *      - second is constant the optimum (10).
 *      - third is constant and not in the optimum range (20).
 * In the third iteration the first and second configuration should be in _optimalSearchSpace and after the tuning phase
 * the second should be the optimal configuration.
 */
TEST_F(PredictiveTuningTest, testLinearPredictionTuningThreeIterations) {
  size_t iteration = 0;
  auto predictiveTuning =
      getPredictiveTuning(evidenceFirstPrediction, autopas::ExtrapolationMethodOption::linePrediction);

  predictiveTuning.reset(iteration);

  auto testedConfigs = tuneForSomeIterationsAndCheckAllTuned(predictiveTuning, {11, 10, 20}, iteration);
  auto nearOptimalConfiguration = testedConfigs[0];
  auto optimalConfiguration = testedConfigs[1];

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  tuneForSomeIterationsAndCheckAllTuned(predictiveTuning, {11, 10, 20}, iteration);

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

/*
 * Tests the selection of the right optimumSearchSpace in the regard of the number of iterations without a test and
 * tuning. Two different configurations:
 *      - first is constant out of the optimum range (20).
 *      - second is constant the optimum (10).
 * In iteration three to six only the second configuration should be in _optimalSearchSpace.
 * In the seventh iteration both configuration should be tested in the tuning phase.
 */
TEST_F(PredictiveTuningTest, testLinearPredictionTooLongNotTested) {
  size_t iteration = 0;
  std::vector<autopas::Configuration> configurationsToCompare{configurationLC_C08, configurationLC_Sliced};
  auto predictiveTuning =
      getPredictiveTuning(evidenceFirstPrediction, autopas::ExtrapolationMethodOption::linePrediction, 0,
                          {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_sliced});

  predictiveTuning.reset(iteration);

  auto testedConfigs =
      tuneForSomeIterationsAndCheckAllTuned(predictiveTuning, {20, 10}, iteration, configurationsToCompare);
  auto badConfiguration = testedConfigs[0];
  auto optimalConfiguration = testedConfigs[1];

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
  for (int i = 0; i < maxTuningIterationsWithoutTest; i++) {
    // End of the (i + 2)-th tuning phase.
    predictiveTuning.reset(iteration);

    EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(10, iteration);
    ++iteration;

    EXPECT_FALSE(predictiveTuning.tune());
    EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
  }

  // End of the (maxTuningIterationsWithoutTest + 2)-th tuning phase.
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

/*
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
      getPredictiveTuning(evidenceFirstPrediction, autopas::ExtrapolationMethodOption::linePrediction);

  predictiveTuning.reset(iteration);

  auto testedConfigs = tuneForSomeIterationsAndCheckAllTuned(predictiveTuning, {15, 10, 20}, iteration);

  auto secondBestConfiguration = testedConfigs[0];
  auto bestConfiguration = testedConfigs[1];
  auto thirdBestConfiguration = testedConfigs[2];

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  tuneForSomeIterationsAndCheckAllTuned(predictiveTuning, {15, 10, 20}, iteration);
  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(bestConfiguration, predictiveTuning.getCurrentConfiguration());

  // mark best as invalid
  predictiveTuning.tune(true);
  ++iteration;

  // Tests if secondBestConfiguration gets selected if the first one is invalid.
  EXPECT_EQ(secondBestConfiguration, predictiveTuning.getCurrentConfiguration());

  // mark second best as invalid
  predictiveTuning.tune(true);
  ++iteration;

  // Tests if thirdBestConfiguration gets selected if the second one is invalid.
  EXPECT_EQ(thirdBestConfiguration, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);

  // Tests if thirdBestConfiguration still gets selected if the second one is invalid.
  predictiveTuning.tune();
  EXPECT_EQ(thirdBestConfiguration, predictiveTuning.getCurrentConfiguration());
}

/*
 * Tests the correct use of the blacklist:
 *      - first is constant the optimum (1).
 *      - second is constant in the blacklist range (5).
 *      - third is  out of the optimum range (20) and should not be tested after the first tuning phase.
 * In the second iteration the third configuration should not be tested in the tuning phase.
 */
TEST_F(PredictiveTuningTest, testBlacklist) {
  size_t iteration = 0;
  auto predictiveTuning =
      getPredictiveTuning(evidenceFirstPrediction, autopas::ExtrapolationMethodOption::linePrediction, 10);

  auto testedConfigs = tuneForSomeIterationsAndCheckAllTuned(predictiveTuning, {1, 5, 20}, iteration);

  auto bestConfiguration = testedConfigs[0];
  auto secondBestConfiguration = testedConfigs[1];

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  tuneForSomeIterationsAndCheckAllTuned(predictiveTuning, {1, 5}, iteration,
                                        {bestConfiguration, secondBestConfiguration});
  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(bestConfiguration, predictiveTuning.getCurrentConfiguration());
}
