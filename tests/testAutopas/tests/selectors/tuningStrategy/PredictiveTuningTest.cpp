/**
 * @file PredictiveTuningTest.cpp
 * @author Julian Pelloth
 * @date 01.04.2020
 */

#include "PredictiveTuningTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_F(PredictiveTuningTest, testSearchSpaceEmpty) {
  autopas::PredictiveTuning predictiveTuning({});
  EXPECT_TRUE(predictiveTuning.searchSpaceIsEmpty());
  EXPECT_FALSE(predictiveTuning.searchSpaceIsTrivial());
  EXPECT_THAT(predictiveTuning.getAllowedContainerOptions(), ::testing::IsEmpty());
}

TEST_F(PredictiveTuningTest, testSearchSpaceOneOption) {
  autopas::PredictiveTuning predictiveTuning(
      {autopas::Configuration(autopas::ContainerOption::directSum, 1., autopas::TraversalOption::directSumTraversal,
                              autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled)});
  EXPECT_FALSE(predictiveTuning.searchSpaceIsEmpty());
  EXPECT_TRUE(predictiveTuning.searchSpaceIsTrivial());
  EXPECT_THAT(predictiveTuning.getAllowedContainerOptions(),
              ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_F(PredictiveTuningTest, testSearchSpaceMoreOptions) {
  autopas::PredictiveTuning predictiveTuning({autopas::ContainerOption::linkedCells}, {1.},
                                             {autopas::TraversalOption::c08}, {autopas::DataLayoutOption::soa},
                                             {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled},
                                             relativeOptimumRange, maxTuningIterationsWithoutTest);
  EXPECT_FALSE(predictiveTuning.searchSpaceIsEmpty());
  EXPECT_FALSE(predictiveTuning.searchSpaceIsTrivial());
  EXPECT_THAT(predictiveTuning.getAllowedContainerOptions(),
              ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

TEST_F(PredictiveTuningTest, testRemoveN3OptionRemoveAll) {
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.}, {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos}, {autopas::Newton3Option::enabled},
      relativeOptimumRange, maxTuningIterationsWithoutTest);

  EXPECT_THROW(predictiveTuning.removeN3Option(autopas::Newton3Option::enabled),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_F(PredictiveTuningTest, testRemoveN3OptionRemoveSome) {
  autopas::PredictiveTuning predictiveTuning({autopas::ContainerOption::linkedCells}, {1.},
                                             {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
                                             {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
                                             {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled},
                                             relativeOptimumRange, maxTuningIterationsWithoutTest);

  EXPECT_NO_THROW(predictiveTuning.removeN3Option(autopas::Newton3Option::enabled));
  EXPECT_FALSE(predictiveTuning.searchSpaceIsEmpty());
  EXPECT_FALSE(predictiveTuning.searchSpaceIsTrivial());
}

/*
 * Tests the prediction and the selection of the right optimum.
 * Three different configurations:
 *      - C08 that is not optimal but it getting less expensive.
 *      - Sliced is optimal in the beginning but is getting more expensive.
 *      - C01 is constant an not in the optimum range.
 * In the third iteration c08 should be predicted to be the optimum.
 */
TEST_F(PredictiveTuningTest, testSelectPossibleConfigurations) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, relativeOptimumRange,
      maxTuningIterationsWithoutTest);

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

// Tests the first tuning iteration. There is no prediction and the whole searchSpace should be tested.
TEST_F(PredictiveTuningTest, testTuneFirstIteration) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, relativeOptimumRange,
      maxTuningIterationsWithoutTest);

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
 * In the third iteration c08 and sliced should be in _optimalSearchSpace and after the tuning phase sliced should be
 * the optimal configuration.
 */
TEST_F(PredictiveTuningTest, testTuningThreeIterations) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, relativeOptimumRange,
      maxTuningIterationsWithoutTest);

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
 * In the seventh iteration c08 and sliced should be in _optimalSearchSpace.
 */
TEST_F(PredictiveTuningTest, testTooLongNotTested) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning({autopas::ContainerOption::linkedCells}, {1.},
                                             {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
                                             {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
                                             relativeOptimumRange, maxTuningIterationsWithoutTest);

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
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, relativeOptimumRange,
      maxTuningIterationsWithoutTest);

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
 * tuning1.2, 5 phase.
 */
TEST_F(PredictiveTuningTest, testInvalidOptimalSearchSpaceTwice) {
  unsigned int iteration = 0;
  autopas::PredictiveTuning predictiveTuning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, relativeOptimumRange,
      maxTuningIterationsWithoutTest);

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
