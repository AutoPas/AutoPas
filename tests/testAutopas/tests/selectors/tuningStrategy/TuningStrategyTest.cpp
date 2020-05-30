/**
 * @file TuningStrategyTest.cpp
 * @author Jan Nguyen
 * @date 27.04.20
 */

#include "TuningStrategyTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

/**
 * Generating a TuningStrategy without a valid configuration is expected to throw
 */
TEST_P(TuningStrategyTest, testSearchSpaceEmpty) {
  auto tuningStrategy = GetParam();
  auto noInterval = autopas::NumberSetFinite<double>({});
  EXPECT_THROW(
      autopas::TuningStrategyFactory::generateTuningStrategy(tuningStrategy, {}, noInterval, {}, {}, {}, {}, 42, 1.2, 5,
                                                             autopas::AcquisitionFunctionOption::expectedImprovement),
      autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_P(TuningStrategyTest, testSearchSpaceOneOption) {
  auto tuningStrategy = GetParam();
  auto oneInterval = autopas::NumberSetFinite<double>({1.});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      tuningStrategy, {autopas::ContainerOption::directSum}, oneInterval,
      {autopas::TraversalOption::directSumTraversal},{autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa},
      {autopas::Newton3Option::enabled}, 42, 1.2, 5, autopas::AcquisitionFunctionOption::expectedImprovement);

  EXPECT_FALSE(search->searchSpaceIsEmpty());
  EXPECT_TRUE(search->searchSpaceIsTrivial());
  EXPECT_THAT(search->getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_P(TuningStrategyTest, testSearchSpaceMoreOptions) {
  auto tuningStrategy = GetParam();
  auto oneInterval = autopas::NumberSetFinite<double>({1.});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      tuningStrategy, {autopas::ContainerOption::linkedCells}, oneInterval, {autopas::TraversalOption::c08},{autopas::LoadEstimatorOption::none},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled}, 1, 1.2, 5,
      autopas::AcquisitionFunctionOption::expectedImprovement);

  EXPECT_FALSE(search->searchSpaceIsEmpty());
  EXPECT_FALSE(search->searchSpaceIsTrivial());
  EXPECT_THAT(search->getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

/**
 * Initialize a search space and remove the only newton3 option.
 * It is expected that this throws.
 */
TEST_P(TuningStrategyTest, testRemoveN3OptionRemoveAll) {
  auto tuningStrategy = GetParam();
  auto oneInterval = autopas::NumberSetFinite<double>({1.});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      tuningStrategy, {autopas::ContainerOption::linkedCells}, oneInterval,
      {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},{autopas::LoadEstimatorOption::none},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos}, {autopas::Newton3Option::enabled}, 1, 1.2, 5,
      autopas::AcquisitionFunctionOption::expectedImprovement);

  EXPECT_THROW(search->removeN3Option(autopas::Newton3Option::enabled),
               autopas::utils::ExceptionHandler::AutoPasException);
}

/**
 * Initialize a search space and remove one of two newton3 options.
 * It is expected that this does not throw and the search space still contains elements.
 */
TEST_P(TuningStrategyTest, testRemoveN3OptionRemoveSome) {
  auto tuningStrategy = GetParam();
  auto oneInterval = autopas::NumberSetFinite<double>({1.});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      tuningStrategy, {autopas::ContainerOption::linkedCells}, oneInterval,
      {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},{autopas::LoadEstimatorOption::none},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
      {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled}, 1, 1.2, 5,
      autopas::AcquisitionFunctionOption::expectedImprovement);

  EXPECT_NO_THROW(search->removeN3Option(autopas::Newton3Option::enabled));
  EXPECT_FALSE(search->searchSpaceIsEmpty());
  EXPECT_FALSE(search->searchSpaceIsTrivial());
}

INSTANTIATE_TEST_SUITE_P(Generated, TuningStrategyTest,
                         testing::ValuesIn(autopas::TuningStrategyOption::getAllOptions()),
                         TuningStrategyTest::PrintToStringParamName());
