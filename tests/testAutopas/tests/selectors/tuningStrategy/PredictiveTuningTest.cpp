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
    EXPECT_THAT(predictiveTuning.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_F(PredictiveTuningTest, testSearchSpaceMoreOptions) {
    autopas::PredictiveTuning predictiveTuning({autopas::ContainerOption::linkedCells}, {1.}, {autopas::TraversalOption::c08},
                                   {autopas::DataLayoutOption::soa},
                                   {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});
    EXPECT_FALSE(predictiveTuning.searchSpaceIsEmpty());
    EXPECT_FALSE(predictiveTuning.searchSpaceIsTrivial());
    EXPECT_THAT(predictiveTuning.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

TEST_F(PredictiveTuningTest, testRemoveN3OptionRemoveAll) {
    autopas::PredictiveTuning predictiveTuning(
            {autopas::ContainerOption::linkedCells}, {1.}, {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
            {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos}, {autopas::Newton3Option::enabled});

    EXPECT_THROW(predictiveTuning.removeN3Option(autopas::Newton3Option::enabled),
                 autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_F(PredictiveTuningTest, testRemoveN3OptionRemoveSome) {
    autopas::PredictiveTuning predictiveTuning({autopas::ContainerOption::linkedCells}, {1.},
                                   {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
                                   {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
                                   {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});

    EXPECT_NO_THROW(predictiveTuning.removeN3Option(autopas::Newton3Option::enabled));
    EXPECT_FALSE(predictiveTuning.searchSpaceIsEmpty());
    EXPECT_FALSE(predictiveTuning.searchSpaceIsTrivial());
}

TEST_F(PredictiveTuningTest, testSelectPossibleConfigurations) {
    autopas::PredictiveTuning predictiveTuning(
            {autopas::ContainerOption::linkedCells}, {1.},
            {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
            {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

    predictiveTuning.reset();

    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(4);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(1);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c01,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(20);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());

    predictiveTuning.reset();

    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(3);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(2);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c01,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(20);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());

    predictiveTuning.reset();

    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
}

TEST_F(PredictiveTuningTest, testTuneFirstIteration) {
    autopas::PredictiveTuning predictiveTuning(
            {autopas::ContainerOption::linkedCells}, {1.},
            {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
            {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

    predictiveTuning.reset();

    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(10);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(1);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c01,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(20);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
}

TEST_F(PredictiveTuningTest, testTuningThreeIterations) {
    autopas::PredictiveTuning predictiveTuning(
            {autopas::ContainerOption::linkedCells}, {1.},
            {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
            {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

    predictiveTuning.reset();

    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(11);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(10);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c01,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(20);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());

    predictiveTuning.reset();

    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(11);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(10);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c01,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(20);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());

    predictiveTuning.reset();

    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(11);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
    predictiveTuning.addEvidence(10);

    predictiveTuning.tune();
    EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                     autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
              predictiveTuning.getCurrentConfiguration());
}

