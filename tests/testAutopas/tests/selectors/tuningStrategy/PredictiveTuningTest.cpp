/**
 * @file PredictiveTuningTest.cpp
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#include "PredictiveTuningTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

/*
 * Tests the prediction and the selection of the right optimum.
 * Three different configurations:
 *      - C08 that is not optimal but it getting less expensive.
 *      - Sliced is optimal in the beginning but is getting more expensive.
 *      - C01 is constant an not in the optimum range.
 * In the third iteration lc_c08 should be predicted to be the optimum.
 */
TEST_F(PredictiveTuningTest, testSelectPossibleConfigurations) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, evidenceFirstPrediction, linePrediction);

  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(4, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(1, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(3, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(2, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  // This tests the actual prediction.
  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
}

/*
 * Tests the prediction and the selection of the right optimum with the linear regression method.
 * Three different configurations:
 *      - C08 that is not optimal but it getting less expensive.
 *      - Sliced is optimal in the beginning but is getting more expensive.
 *      - C01 is constant an not in the optimum range.
 * In the third iteration lc_c08 should be predicted to be the optimum.
 */
TEST_F(PredictiveTuningTest, testLinearRegression) {
  unsigned int iteration = 1;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, 2, autopas::ExtrapolationMethodOption::linearRegression);

  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(375, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(300, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(2000, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(350, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(325, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(2000, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  // This tests the actual prediction.
  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
}

/*
 * Tests the prediction and the selection of the right optimum with the lagrange interpolation method.
 * Three different configurations:
 *      - C08 that is not optimal but it getting less expensive.
 *      - Sliced is optimal in the beginning but is getting more expensive.
 *      - C01 is constant an not in the optimum range.
 * In the third iteration lc_c08 should be predicted to be the optimum.
 */
TEST_F(PredictiveTuningTest, testLagrange) {
  unsigned int iteration = 1;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, 3, autopas::ExtrapolationMethodOption::lagrange);

  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(6, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(1, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(5, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(2, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(4, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(3, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the third tuning phase.
  predictiveTuning.reset(iteration);

  // This tests the actual prediction.
  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
}

/*
 * Tests the prediction and the selection of the right optimum with the newton interpolation method.
 * Three different configurations:
 *      - C08 that is not optimal but it getting less expensive.
 *      - Sliced is optimal in the beginning but is getting more expensive.
 *      - C01 is constant an not in the optimum range.
 * In the third iteration lc_c08 should be predicted to be the optimum.
 */
TEST_F(PredictiveTuningTest, testNewton) {
  unsigned int iteration = 1;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, 3, autopas::ExtrapolationMethodOption::newton);

  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(60, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(200, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(50, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(200, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(40, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(30, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(200, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the third tuning phase.
  predictiveTuning.reset(iteration);

  // This tests the actual prediction.
  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
}

// Tests the first tuning iteration. There is no prediction and the whole searchSpace should be tested.
TEST_F(PredictiveTuningTest, testTuneFirstIteration) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, evidenceFirstPrediction, linePrediction);

  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(1, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
}

/*
 * Tests the selection of the right optimumSearchSpace in the regard of the relativeOptimum and tuning.
 * Three different configurations:
 *      - C08 is constant near the optimum (11).
 *      - Sliced is constant the optimum (10).
 *      - C01 is constant an not in the optimum range (20).
 * In the third iteration lc_c08 and sliced should be in _optimalSearchSpace and after the tuning phase sliced should be
 * the optimal configuration.
 */
TEST_F(PredictiveTuningTest, testTuningThreeIterations) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, evidenceFirstPrediction, linePrediction);

  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(11, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(11, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  // This test if a configuration near the optimum gets selected into _optimalSearchSpace.
  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(11, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);

  // This tests that the right optimum configuration gets selected at the end of a tuning phase.
  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
}

/*
 * Tests the selection of the right optimumSearchSpace in the regard of the number of iterations without a test and
 * tuning. Two different configurations:
 *      - C08 is constant out of the optimum range (20).
 *      - Sliced is constant the optimum (10).
 * In iteration three to six only sliced should be in _optimalSearchSpace.
 * In the seventh iteration lc_c08 and sliced should be in _optimalSearchSpace.
 */
TEST_F(PredictiveTuningTest, testTooLongNotTested) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning({autopas::ContainerOption::linkedCells}, {1.},
                                             {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_sliced},
                                             {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa},
                                             {autopas::Newton3Option::disabled}, relativeOptimumRange,
                                             maxTuningIterationsWithoutTest, evidenceFirstPrediction, linePrediction);

  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // Iterating through the maxNumberOfIterationsWithoutTest iterations where C08 should not be tested.
  for (int i = 0; i < maxTuningIterationsWithoutTest; i++) {
    // End of the (i + 2)-th tuning phase.
    predictiveTuning.reset(iteration);

    EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(10, iteration);
    ++iteration;

    predictiveTuning.tune();
    EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  }

  // End of the (maxTuningIterationsWithoutTest + 2)-th tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  // Tests that a configuration gets selected into _tooLongNotTested.
  predictiveTuning.tune();
  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
}

/*
 * Tests the recognition of an invalid _optimalSearchSpace in tuning and the reselection of _optimalSearchSpace
 * Three different configurations:
 *      - C08 is constant out of the optimum range (15).
 *      - Sliced is constant the optimum (10) and invalid in the third iteration.
 *      - C01 is constant out of the optimum range (20).
 * In the third iteration reselectOptimalSearchSpace should be called and C08 should be selected after the tuning phase.
 */
TEST_F(PredictiveTuningTest, testInvalidOptimalSearchSpaceOnce) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, evidenceFirstPrediction, linePrediction);

  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(15, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(15, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  ++iteration;

  predictiveTuning.tune(true);

  // Tests if a new _optimalSearchSpace gets selected if the first one is invalid.
  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(15, iteration);

  predictiveTuning.tune();
  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
}

/*
 * Tests the recognition of an invalid _optimalSearchSpace in tuning and the reselection of _optimalSearchSpace
 * Three different configurations:
 *      - C08 is constant out of the optimum range (15) and invalid in the third iteration.
 *      - Sliced is constant the optimum (10) and invalid in the third iteration.
 *      - C01 is constant out of the optimum range (20).
 * In the third iteration reselectOptimalSearchSpace should be called twice and C01 should be selected after the
 * tuning phase.
 */
TEST_F(PredictiveTuningTest, testInvalidOptimalSearchSpaceTwice) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest, evidenceFirstPrediction, linePrediction);

  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(15, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the first tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(15, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(10, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);
  ++iteration;

  predictiveTuning.tune();
  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());

  // End of the second tuning phase.
  predictiveTuning.reset(iteration);

  EXPECT_EQ(configurationSliced, predictiveTuning.getCurrentConfiguration());
  ++iteration;

  predictiveTuning.tune(true);

  // Tests if a new _optimalSearchSpace gets selected if the first one is invalid.
  EXPECT_EQ(configurationC08, predictiveTuning.getCurrentConfiguration());
  ++iteration;

  predictiveTuning.tune(true);

  // Tests if a new _optimalSearchSpace gets selected if the second one is invalid.
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
  predictiveTuning.addEvidence(20, iteration);

  predictiveTuning.tune();
  EXPECT_EQ(configurationC01, predictiveTuning.getCurrentConfiguration());
}
