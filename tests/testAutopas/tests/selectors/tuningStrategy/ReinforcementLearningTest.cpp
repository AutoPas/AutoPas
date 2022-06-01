/**
 * @file ReinforcementLearningTest.cpp
 * @author L. Laumeyer
 * @date 30/5/22
 */

#include "ReinforcementLearningTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_F(ReinforcementLearningTest, testSearchSpaceEmpty) {
  autopas::ReinforcementLearning reinforcementLearning({});
  EXPECT_TRUE(reinforcementLearning.searchSpaceIsEmpty());
  EXPECT_FALSE(reinforcementLearning.searchSpaceIsTrivial());
  EXPECT_THAT(reinforcementLearning.getAllowedContainerOptions(), ::testing::IsEmpty());
}

TEST_F(ReinforcementLearningTest, testTuneSingleRun) {
  autopas::ReinforcementLearning reinforcementLearning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

  std::vector<autopas::Configuration> testedConfigs;
  std::vector<autopas::Configuration> allConfigs{
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled),
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled),
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled)};

  testedConfigs.emplace_back(reinforcementLearning.getCurrentConfiguration());
  reinforcementLearning.addEvidence(10, 0);

  reinforcementLearning.tune();
  testedConfigs.emplace_back(reinforcementLearning.getCurrentConfiguration());
  autopas::Configuration optimalConfig = reinforcementLearning.getCurrentConfiguration();
  reinforcementLearning.addEvidence(1, 0);

  reinforcementLearning.tune();
  testedConfigs.emplace_back(reinforcementLearning.getCurrentConfiguration());
  reinforcementLearning.addEvidence(20, 0);

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));
  reinforcementLearning.tune();
  EXPECT_EQ(optimalConfig, reinforcementLearning.getCurrentConfiguration());
}

TEST_F(ReinforcementLearningTest, testTuneDoubleRun) {
  autopas::ReinforcementLearning reinforcementLearning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

  std::vector<autopas::Configuration> testedConfigs;
  std::vector<autopas::Configuration> allConfigs{
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled),
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled),
      autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                             autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                             autopas::Newton3Option::disabled)};

  testedConfigs.emplace_back(reinforcementLearning.getCurrentConfiguration());
  reinforcementLearning.addEvidence(2, 0);

  reinforcementLearning.tune();
  testedConfigs.emplace_back(reinforcementLearning.getCurrentConfiguration());
  autopas::Configuration optimalConfig = reinforcementLearning.getCurrentConfiguration();
  reinforcementLearning.addEvidence(1, 0);

  reinforcementLearning.tune();
  testedConfigs.emplace_back(reinforcementLearning.getCurrentConfiguration());
  reinforcementLearning.addEvidence(20000, 0);

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));
  reinforcementLearning.tune();
  EXPECT_EQ(optimalConfig, reinforcementLearning.getCurrentConfiguration());

  //  Second Run
  optimalConfig = reinforcementLearning.getCurrentConfiguration();
  reinforcementLearning.addEvidence(5, 0);

  reinforcementLearning.tune();
  reinforcementLearning.addEvidence(100, 0);

  reinforcementLearning.tune();
  reinforcementLearning.addEvidence(4, 0);

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));
  reinforcementLearning.tune();
  EXPECT_EQ(optimalConfig, reinforcementLearning.getCurrentConfiguration());
}
