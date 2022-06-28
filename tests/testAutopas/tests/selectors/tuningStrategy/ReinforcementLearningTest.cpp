/**
 * @file ReinforcementLearningTest.cpp
 * @author L. Laumeyer
 * @date 30/5/22
 */

#include "ReinforcementLearningTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

/*
 * Test when there is no valid configuration
 */
TEST_F(ReinforcementLearningTest, testSearchSpaceEmpty) {
  autopas::ReinforcementLearning reinforcementLearning({});
  EXPECT_TRUE(reinforcementLearning.searchSpaceIsEmpty());
  EXPECT_FALSE(reinforcementLearning.searchSpaceIsTrivial());
  EXPECT_THAT(reinforcementLearning.getAllowedContainerOptions(), ::testing::IsEmpty());
}

/*
 * Test for selecting optimal configuration on the first tuning phase
 */
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

/*
 * Test for selectiing optimal configuration on the second tuning phase with switching configurations
 */
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
  reinforcementLearning.reset(4);

  reinforcementLearning.addEvidence(20000, 1);
  reinforcementLearning.tune();

  optimalConfig = reinforcementLearning.getCurrentConfiguration();
  reinforcementLearning.addEvidence(1, 1);
  reinforcementLearning.tune();

  reinforcementLearning.addEvidence(20000, 1);

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));
  EXPECT_THAT(false, reinforcementLearning.tune());
  EXPECT_EQ(optimalConfig, reinforcementLearning.getCurrentConfiguration());
}

/*
 * Test for correcting switching of the bools which should guarantee the correct sequence of if clauses in the tune
 * function
 */
TEST_F(ReinforcementLearningTest, testBools) {
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

  EXPECT_THAT(true, reinforcementLearning.getFirstTuningPhase());
  reinforcementLearning.addEvidence(2, 0);

  reinforcementLearning.tune();
  reinforcementLearning.addEvidence(1, 0);

  reinforcementLearning.tune();
  reinforcementLearning.addEvidence(20000, 0);

  reinforcementLearning.tune();

  //  Second Run
  reinforcementLearning.reset(4);

  EXPECT_THAT(true, reinforcementLearning.getStartTuningPhase());
  EXPECT_THAT(false, reinforcementLearning.getFirstTuningPhase());
  reinforcementLearning.addEvidence(5, 1);

  reinforcementLearning.tune();
  EXPECT_THAT(false, reinforcementLearning.getStartTuningPhase());
  reinforcementLearning.addEvidence(100, 1);

  reinforcementLearning.tune();
  reinforcementLearning.addEvidence(4, 1);

  reinforcementLearning.tune();
  EXPECT_THAT(true, reinforcementLearning.getStartTuningPhase());
  EXPECT_THAT(false, reinforcementLearning.getFirstTuningPhase());
}

/*
 * Test for the correct adding of new Timings
 */
TEST_F(ReinforcementLearningTest, testAddEvidence) {
  autopas::ReinforcementLearning reinforcementLearning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

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

  reinforcementLearning.addEvidence(10, 1);
  EXPECT_THAT(-10, reinforcementLearning.getState(reinforcementLearning.getCurrentConfiguration()));
  reinforcementLearning.tune();

  autopas::Configuration optimalConfig = reinforcementLearning.getCurrentConfiguration();
  reinforcementLearning.addEvidence(15, 2);
  EXPECT_THAT(-15, reinforcementLearning.getState(reinforcementLearning.getCurrentConfiguration()));
  reinforcementLearning.tune();

  reinforcementLearning.addEvidence(20, 3);
  EXPECT_THAT(-20, reinforcementLearning.getState(reinforcementLearning.getCurrentConfiguration()));
  reinforcementLearning.tune();

  //  Second Run
  reinforcementLearning.reset(4);

  reinforcementLearning.addEvidence(5, 4);
  reinforcementLearning.tune();

  reinforcementLearning.addEvidence(5, 5);
  reinforcementLearning.tune();

  reinforcementLearning.addEvidence(5, 6);
  EXPECT_THAT(-6, reinforcementLearning.getState(reinforcementLearning.getCurrentConfiguration()));
  reinforcementLearning.tune();
}

/*
 * Test for the correct State value after the first and second tuning phase
 */
TEST_F(ReinforcementLearningTest, testFindState) {
  autopas::ReinforcementLearning reinforcementLearning(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

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
  reinforcementLearning.addEvidence(10, 1);
  reinforcementLearning.tune();

  reinforcementLearning.addEvidence(15, 2);
  reinforcementLearning.tune();

  reinforcementLearning.addEvidence(20, 3);
  reinforcementLearning.tune();

  EXPECT_THAT(-20, reinforcementLearning._states.find(allConfigs.at(0))->second);

  reinforcementLearning.reset(4);

  reinforcementLearning.addEvidence(5, 4);
  reinforcementLearning.tune();

  reinforcementLearning.addEvidence(5, 5);
  reinforcementLearning.tune();

  reinforcementLearning.addEvidence(5, 6);
  reinforcementLearning.tune();
  auto config = reinforcementLearning.getCurrentConfiguration();
  double test = reinforcementLearning.getState(reinforcementLearning.getCurrentConfiguration());
  EXPECT_THAT(-5.5, reinforcementLearning.getState(reinforcementLearning.getCurrentConfiguration()));
}

/*
 * Test for selecting the correct timing in the first and second tuning phase
 */
TEST_F(ReinforcementLearningTest, testGetEvidence) {
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
  reinforcementLearning.addEvidence(22, 0);

  reinforcementLearning.tune();
  testedConfigs.emplace_back(reinforcementLearning.getCurrentConfiguration());
  reinforcementLearning.addEvidence(11, 1);

  reinforcementLearning.tune();
  testedConfigs.emplace_back(reinforcementLearning.getCurrentConfiguration());
  reinforcementLearning.addEvidence(33, 2);

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));
  EXPECT_THAT(33, reinforcementLearning.getEvidence(allConfigs[0]));
  EXPECT_THAT(22, reinforcementLearning.getEvidence(allConfigs[1]));
  EXPECT_THAT(11, reinforcementLearning.getEvidence(allConfigs[2]));
  reinforcementLearning.tune();

  //  Second Run
  reinforcementLearning.reset(4);
  reinforcementLearning.addEvidence(1, 4);
  EXPECT_THAT(allConfigs[1], reinforcementLearning.getCurrentConfiguration());
  reinforcementLearning.tune();

  EXPECT_THAT(allConfigs[2], reinforcementLearning.getCurrentConfiguration());
  reinforcementLearning.addEvidence(2, 5);
  reinforcementLearning.tune();

  EXPECT_THAT(allConfigs[0], reinforcementLearning.getCurrentConfiguration());
  reinforcementLearning.addEvidence(3, 6);

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));
  EXPECT_THAT(false, reinforcementLearning.tune());
  EXPECT_THAT(3, reinforcementLearning.getEvidence(allConfigs[0]));
  EXPECT_THAT(1, reinforcementLearning.getEvidence(allConfigs[1]));
  EXPECT_THAT(2, reinforcementLearning.getEvidence(allConfigs[2]));
}