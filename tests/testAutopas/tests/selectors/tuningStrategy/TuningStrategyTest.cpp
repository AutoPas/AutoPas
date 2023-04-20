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
  auto noContainers = std::set<autopas::ContainerOption>({});
  auto noTraversals = std::set<autopas::TraversalOption>({});
  auto noLoadEstimators = std::set<autopas::LoadEstimatorOption>({});
  auto noDataLayouts = std::set<autopas::DataLayoutOption>({});
  auto noNewton3Options = std::set<autopas::Newton3Option>({});
  auto noFreq = autopas::NumberSetFinite<int>({});
  EXPECT_THROW(autopas::TuningStrategyFactory::generateTuningStrategy(
                   tuningStrategy, noContainers, noInterval, noTraversals, noLoadEstimators, noDataLayouts,
                   noNewton3Options, noFreq, 42, 1.2, 5, 0, 3, autopas::AcquisitionFunctionOption::expectedImprovement,
                   autopas::ExtrapolationMethodOption::linePrediction),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_P(TuningStrategyTest, testSearchSpaceOneOption) {
  auto tuningStrategy = GetParam();
  auto oneInterval = autopas::NumberSetFinite<double>({1.});
  auto oneContainer = std::set<autopas::ContainerOption>({autopas::ContainerOption::directSum});
  auto oneTraversal = std::set<autopas::TraversalOption>({autopas::TraversalOption::ds_sequential});
  auto oneLoadEstimator = std::set<autopas::LoadEstimatorOption>({autopas::LoadEstimatorOption::none});
  auto oneDataLayout = std::set<autopas::DataLayoutOption>({autopas::DataLayoutOption::soa});
  auto oneNewton3Option = std::set<autopas::Newton3Option>({autopas::Newton3Option::enabled});
  auto oneFreq = autopas::NumberSetFinite<int>({15});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      tuningStrategy, oneContainer, oneInterval, oneTraversal, oneLoadEstimator, oneDataLayout, oneNewton3Option,
      oneFreq, 42, 1.2, 5, 0, 3, autopas::AcquisitionFunctionOption::expectedImprovement,
      autopas::ExtrapolationMethodOption::linePrediction);

  EXPECT_FALSE(search->searchSpaceIsEmpty());
  EXPECT_TRUE(search->searchSpaceIsTrivial());
  EXPECT_THAT(search->getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_P(TuningStrategyTest, testSearchSpaceMoreOptions) {
  auto tuningStrategy = GetParam();
  auto oneInterval = autopas::NumberSetFinite<double>({1.});
  auto oneContainer = std::set<autopas::ContainerOption>({autopas::ContainerOption::linkedCells});
  auto oneTraversal = std::set<autopas::TraversalOption>({autopas::TraversalOption::lc_c08});
  auto oneLoadEstimator = std::set<autopas::LoadEstimatorOption>({autopas::LoadEstimatorOption::none});
  auto oneDataLayout = std::set<autopas::DataLayoutOption>({autopas::DataLayoutOption::soa});
  auto twoNewton3Options =
      std::set<autopas::Newton3Option>({autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});
  auto oneFreq = autopas::NumberSetFinite<int>({5});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      tuningStrategy, oneContainer, oneInterval, oneTraversal, oneLoadEstimator, oneDataLayout, twoNewton3Options,
      oneFreq, 1, 1.2, 5, 0, 3, autopas::AcquisitionFunctionOption::expectedImprovement,
      autopas::ExtrapolationMethodOption::linePrediction);

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
  auto oneContainer = std::set<autopas::ContainerOption>({autopas::ContainerOption::linkedCells});
  auto twoTraversals =
      std::set<autopas::TraversalOption>({autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_sliced});
  auto oneLoadEstimator = std::set<autopas::LoadEstimatorOption>({autopas::LoadEstimatorOption::none});
  auto twoDataLayouts =
      std::set<autopas::DataLayoutOption>({autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos});
  auto oneNewton3Option = std::set<autopas::Newton3Option>({autopas::Newton3Option::enabled});
  auto oneFreq = autopas::NumberSetFinite<int>({5});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      tuningStrategy, oneContainer, oneInterval, twoTraversals, oneLoadEstimator, twoDataLayouts, oneNewton3Option,
      oneFreq, 1, 1.2, 5, 0, 3, autopas::AcquisitionFunctionOption::expectedImprovement,
      autopas::ExtrapolationMethodOption::linePrediction);

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
  auto oneContainer = std::set<autopas::ContainerOption>({autopas::ContainerOption::linkedCells});
  auto twoTraversals =
      std::set<autopas::TraversalOption>({autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_sliced});
  auto oneLoadEstimator = std::set<autopas::LoadEstimatorOption>({autopas::LoadEstimatorOption::none});
  auto twoDataLayouts =
      std::set<autopas::DataLayoutOption>({autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos});
  auto twoNewton3Options =
      std::set<autopas::Newton3Option>({autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});
  auto oneFreq = autopas::NumberSetFinite<int>({5});
  auto search = autopas::TuningStrategyFactory::generateTuningStrategy(
      tuningStrategy, oneContainer, oneInterval, twoTraversals, oneLoadEstimator, twoDataLayouts, twoNewton3Options,
      oneFreq, 1, 1.2, 5, 0, 3, autopas::AcquisitionFunctionOption::expectedImprovement,
      autopas::ExtrapolationMethodOption::linePrediction);

  EXPECT_NO_THROW(search->removeN3Option(autopas::Newton3Option::enabled));
  EXPECT_FALSE(search->searchSpaceIsEmpty());
  EXPECT_FALSE(search->searchSpaceIsTrivial());
}

INSTANTIATE_TEST_SUITE_P(Generated, TuningStrategyTest,
                         testing::ValuesIn(autopas::TuningStrategyOption::getAllOptions()),
                         TuningStrategyTest::PrintToStringParamName());
