/**
 * @file ActiveHarmonyTest.cpp
 * @author Jakob Englhauser
 * @date 14.11.19
 */

#include "ActiveHarmonyTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_F(ActiveHarmonyTest, testSearchSpaceEmpty) {
  autopas::ActiveHarmony activeHarmony(std::set<autopas::ContainerOption>({}), autopas::NumberInterval<double>(),
                                       std::set<autopas::TraversalOption>({}),
                                       std::set<autopas::LoadEstimatorOption>({}),
                                       std::set<autopas::DataLayoutOption>({}), std::set<autopas::Newton3Option>({}));

  EXPECT_TRUE(activeHarmony.searchSpaceIsEmpty());
  EXPECT_FALSE(activeHarmony.searchSpaceIsTrivial());
  EXPECT_THAT(activeHarmony.getAllowedContainerOptions(), ::testing::IsEmpty());
}

TEST_F(ActiveHarmonyTest, testSearchSpaceOneOption) {
  autopas::ActiveHarmony activeHarmony({autopas::ContainerOption::directSum}, autopas::NumberInterval<double>(),
                                       {autopas::TraversalOption::directSumTraversal},
                                       {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa},
                                       {autopas::Newton3Option::enabled});

  EXPECT_FALSE(activeHarmony.searchSpaceIsEmpty());
  EXPECT_TRUE(activeHarmony.searchSpaceIsTrivial());
  EXPECT_THAT(activeHarmony.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_F(ActiveHarmonyTest, testSearchSpaceMoreOptions) {
  autopas::ActiveHarmony activeHarmony({autopas::ContainerOption::linkedCells}, autopas::NumberInterval<double>(),
                                       {autopas::TraversalOption::c08}, {autopas::LoadEstimatorOption::none},
                                       {autopas::DataLayoutOption::soa}, autopas::Newton3Option::getAllOptions());

  EXPECT_FALSE(activeHarmony.searchSpaceIsEmpty());
  EXPECT_FALSE(activeHarmony.searchSpaceIsTrivial());
  EXPECT_THAT(activeHarmony.getAllowedContainerOptions(),
              ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

TEST_F(ActiveHarmonyTest, testRemoveN3OptionRemoveAll) {
  autopas::ActiveHarmony activeHarmony({autopas::ContainerOption::linkedCells}, autopas::NumberInterval<double>(),
                                       {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
                                       {autopas::LoadEstimatorOption::none}, autopas::DataLayoutOption::getAllOptions(),
                                       {autopas::Newton3Option::enabled});

  EXPECT_THROW(activeHarmony.removeN3Option(autopas::Newton3Option::enabled),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_F(ActiveHarmonyTest, testRemoveN3OptionRemoveSome) {
  autopas::ActiveHarmony activeHarmony({autopas::ContainerOption::linkedCells}, autopas::NumberInterval<double>(),
                                       {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
                                       {autopas::LoadEstimatorOption::none}, autopas::DataLayoutOption::getAllOptions(),
                                       autopas::Newton3Option::getAllOptions());

  EXPECT_NO_THROW(activeHarmony.removeN3Option(autopas::Newton3Option::enabled));
  EXPECT_FALSE(activeHarmony.searchSpaceIsEmpty());
  EXPECT_FALSE(activeHarmony.searchSpaceIsTrivial());
}
