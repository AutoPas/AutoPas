/**
 * @file TuningStrategyTest.cpp
 * @author Jan Nguyen
 * @date 27.04.20
 */

#include "TuningStrategyTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_P(TuningStrategyTest, testSearchSpaceEmpty) {
  auto noInterval = autopas::NumberSetFinite<double>({});
  EXPECT_THROW(autopas::TuningStrategyFactory::generateTuningStrategy(
                   GetParam(), {}, noInterval, {}, {}, {}, 42, autopas::AcquisitionFunctionOption::expectedDecrease),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_P(TuningStrategyTest, testSearchSpaceOneOption) {
  auto oneInterval = autopas::NumberSetFinite<double>({1.});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      GetParam(), {autopas::ContainerOption::directSum}, oneInterval, {autopas::TraversalOption::directSumTraversal},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::enabled}, 42,
      autopas::AcquisitionFunctionOption::expectedDecrease);

  EXPECT_FALSE(search->searchSpaceIsEmpty());
  EXPECT_TRUE(search->searchSpaceIsTrivial());
  EXPECT_THAT(search->getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_P(TuningStrategyTest, testSearchSpaceMoreOptions) {
  auto oneInterval = autopas::NumberSetFinite<double>({1.});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      GetParam(), {autopas::ContainerOption::linkedCells}, oneInterval, {autopas::TraversalOption::c08},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled}, 1,
      autopas::AcquisitionFunctionOption::expectedDecrease);

  EXPECT_FALSE(search->searchSpaceIsEmpty());
  EXPECT_FALSE(search->searchSpaceIsTrivial());
  EXPECT_THAT(search->getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

TEST_P(TuningStrategyTest, testRemoveN3OptionRemoveAll) {
  auto oneInterval = autopas::NumberSetFinite<double>({1.});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      GetParam(), {autopas::ContainerOption::linkedCells}, oneInterval,
      {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos}, {autopas::Newton3Option::enabled}, 1,
      autopas::AcquisitionFunctionOption::expectedDecrease);

  EXPECT_THROW(search->removeN3Option(autopas::Newton3Option::enabled),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_P(TuningStrategyTest, testRemoveN3OptionRemoveSome) {
  auto oneInterval = autopas::NumberSetFinite<double>({1.});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      GetParam(), {autopas::ContainerOption::linkedCells}, oneInterval,
      {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
      {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled}, 1,
      autopas::AcquisitionFunctionOption::expectedDecrease);

  EXPECT_NO_THROW(search->removeN3Option(autopas::Newton3Option::enabled));
  EXPECT_FALSE(search->searchSpaceIsEmpty());
  EXPECT_FALSE(search->searchSpaceIsTrivial());
}

INSTANTIATE_TEST_SUITE_P(Generated, TuningStrategyTest,
                         testing::ValuesIn(autopas::TuningStrategyOption::getAllOptions()),
                         TuningStrategyTest::PrintToStringParamName());
