/**
 * @file PredictiveTuningTest.cpp
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#include "PredictiveTuningTest.h"

#include <gmock/gmock-matchers.h>

autopas::PredictiveTuning PredictiveTuningTest::getPredictiveTuning(
    unsigned int testsUntilFirstPrediction, autopas::ExtrapolationMethodOption extrapolationMethodOption) {
  return autopas::PredictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, testsUntilFirstPrediction, extrapolationMethodOption);
}

void PredictiveTuningTest::testGeneric(autopas::ExtrapolationMethodOption extrapolationMethodOption,
                                       const std::vector<std::array<long, 3>> &evidences,
                                       size_t optimalPredictionIndex) {
  size_t iteration = 0;

  auto predictiveTuning = getPredictiveTuning(evidences.size(), extrapolationMethodOption);

  // First reset tuning.
  predictiveTuning.reset(iteration);

  autopas::Configuration optimalPrediction;

  for (size_t index = 0; index < evidences.size(); ++index) {
    if (index == 0) {
      optimalPrediction = tuneForSomeIterationsAndCheckAllTuned<1>(predictiveTuning, evidences[index],
                                                                   {optimalPredictionIndex}, iteration)[0];
    } else {
      tuneForSomeIterationsAndCheckAllTuned<0>(predictiveTuning, evidences[index], {}, iteration)[0];
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
 * In the third iteration lc_c08 should be predicted to be the optimum.
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
 * In the third iteration lc_c08 should be predicted to be the optimum.
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
 * In the third iteration lc_c08 should be predicted to be the optimum.
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
 * In the third iteration lc_c08 should be predicted to be the optimum.
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
 * In the third iteration lc_c08 and sliced should be in _optimalSearchSpace and after the tuning phase sliced should be
 * the optimal configuration.
 */
TEST_F(PredictiveTuningTest, testLinearPredictionTuningThreeIterations) {
  size_t iteration = 0;
  auto predictiveTuning =
      getPredictiveTuning(evidenceFirstPrediction, autopas::ExtrapolationMethodOption::linePrediction);

  predictiveTuning.reset(iteration);

  auto [nearOptimalConfiguration, optimalConfiguration] =
      tuneForSomeIterationsAndCheckAllTuned<2,3>(predictiveTuning, {11, 10, 20}, {0, 1}, iteration);
  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  tuneForSomeIterationsAndCheckAllTuned<0,3>(predictiveTuning, {11, 10, 20}, {}, iteration);

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
 * In iteration three to six only sliced should be in _optimalSearchSpace.
 * In the seventh iteration lc_c08 and sliced should be in _optimalSearchSpace.
 */
TEST_F(PredictiveTuningTest, testLinearPredictionTooLongNotTested) {
  unsigned int iteration = 0;
  std::vector<autopas::Configuration> configurationsToCompare{configurationLC_C08, configurationLC_Sliced};
  std::vector<autopas::Configuration> testedConfigs;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_sliced}, {autopas::LoadEstimatorOption::none},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, relativeOptimumRange,
      maxTuningIterationsWithoutTest, evidenceFirstPrediction, autopas::ExtrapolationMethodOption::linePrediction);

  predictiveTuning.reset(iteration);

  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  auto badConfiguration = predictiveTuning.getCurrentConfiguration();
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  auto optimalConfiguration = predictiveTuning.getCurrentConfiguration();
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  EXPECT_THAT(configurationsToCompare, testing::UnorderedElementsAreArray(testedConfigs));

  predictiveTuning.tune();
  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());

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

    predictiveTuning.tune();
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
 *      - first is constant out of the optimum range (15).
 *      - second is constant the optimum (10) and invalid in the third iteration.
 *      - third is constant out of the optimum range (20).
 * In the third iteration reselectOptimalSearchSpace should be called and C08 should be selected after the tuning phase.
 */
TEST_F(PredictiveTuningTest, testInvalidOptimalSearchSpaceOnce) {
  unsigned int iteration = 0;
  std::vector<autopas::Configuration> testedConfigs;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, evidenceFirstPrediction,
      autopas::ExtrapolationMethodOption::linePrediction);

  predictiveTuning.reset(iteration);

  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  auto newOptimalConfiguration = predictiveTuning.getCurrentConfiguration();
  predictiveTuning.addEvidence(15, iteration);
  ++iteration;

  predictiveTuning.tune();
  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  auto optimalConfiguration = predictiveTuning.getCurrentConfiguration();
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));

  predictiveTuning.tune();
  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);
  testedConfigs.clear();

  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(15, iteration);
  ++iteration;

  predictiveTuning.tune();
  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));

  predictiveTuning.tune();
  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());

  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
  ++iteration;

  predictiveTuning.tune(true);

  // Tests if a new _optimalSearchSpace gets selected if the first one is invalid.
  EXPECT_EQ(newOptimalConfiguration, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(15, iteration);

  predictiveTuning.tune();
  EXPECT_EQ(newOptimalConfiguration, predictiveTuning.getCurrentConfiguration());
}

/*
 * Tests the recognition of an invalid _optimalSearchSpace in tuning and the reselection of _optimalSearchSpace
 * Three different configurations:
 *      - first is constant out of the optimum range (15) and invalid in the third iteration.
 *      - second is constant the optimum (10) and invalid in the third iteration.
 *      - third is constant out of the optimum range (20).
 * In the third iteration reselectOptimalSearchSpace should be called twice and C01 should be selected after the
 * tuning phase.
 */
TEST_F(PredictiveTuningTest, testInvalidOptimalSearchSpaceTwice) {
  unsigned int iteration = 0;
  std::vector<autopas::Configuration> testedConfigs;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, evidenceFirstPrediction,
      autopas::ExtrapolationMethodOption::linePrediction);

  predictiveTuning.reset(iteration);

  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(15, iteration);
  ++iteration;

  predictiveTuning.tune();
  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  auto optimalConfiguration = predictiveTuning.getCurrentConfiguration();
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  auto newOptimalConfiguration = predictiveTuning.getCurrentConfiguration();
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));

  predictiveTuning.tune();
  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);
  testedConfigs.clear();

  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(15, iteration);
  ++iteration;

  predictiveTuning.tune();
  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  testedConfigs.emplace_back(predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));

  predictiveTuning.tune();
  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());

  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(optimalConfiguration, predictiveTuning.getCurrentConfiguration());
  ++iteration;

  predictiveTuning.tune(true);

  // Tests if a new _optimalSearchSpace gets selected if the first one is invalid.
  // EXPECT_EQ(, predictiveTuning.getCurrentConfiguration());
  ++iteration;

  predictiveTuning.tune(true);

  // Tests if a new _optimalSearchSpace gets selected if the second one is invalid.
  EXPECT_EQ(newOptimalConfiguration, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);

  predictiveTuning.tune();
  EXPECT_EQ(newOptimalConfiguration, predictiveTuning.getCurrentConfiguration());
}
