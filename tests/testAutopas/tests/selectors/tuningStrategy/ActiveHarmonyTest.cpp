/**
 * @file ActiveHarmonyTest.cpp
 * @author Jakob Englhauser
 * @date 14.11.19
 */

#include "ActiveHarmonyTest.h"
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_F(ActiveHarmonyTest, testSearchSpaceEmpty) {
  autopas::ActiveHarmony activeHarmony(std::set<autopas::ContainerOption>({}), autopas::NumberInterval<double>(), std::set<autopas::TraversalOption>({}), std::set<autopas::DataLayoutOption>({}), std::set<autopas::Newton3Option>({}));
  EXPECT_TRUE(activeHarmony.searchSpaceIsEmpty());
  EXPECT_FALSE(activeHarmony.searchSpaceIsTrivial());
  EXPECT_THAT(activeHarmony.getAllowedContainerOptions(), ::testing::IsEmpty());
}