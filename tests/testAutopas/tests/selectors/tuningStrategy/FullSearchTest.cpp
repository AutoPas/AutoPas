/**
 * @file FullSearchTest.cpp
 * @author F. Gratl
 * @date 6/5/19
 */

#include "FullSearchTest.h"
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_F(FullSearchTest, testSearchSpaceEmpty) {
  autopas::FullSearch fullSearch({});
  EXPECT_TRUE(fullSearch.searchSpaceEmpty());
  EXPECT_FALSE(fullSearch.searchSpaceOneOption());
  EXPECT_THAT(fullSearch.getAllowedContainerOptions(), ::testing::IsEmpty());
}

TEST_F(FullSearchTest, testSearchSpaceOneOption) {
  autopas::FullSearch fullSearch(
      {autopas::Configuration(autopas::ContainerOption::directSum, autopas::TraversalOption::directSumTraversal,
                              autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled)});
  EXPECT_FALSE(fullSearch.searchSpaceEmpty());
  EXPECT_TRUE(fullSearch.searchSpaceOneOption());
  EXPECT_THAT(fullSearch.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_F(FullSearchTest, testSearchSpaceMoreOptions) {
  autopas::FullSearch fullSearch({autopas::ContainerOption::linkedCells}, {autopas::TraversalOption::c08},
                                 {autopas::DataLayoutOption::soa},
                                 {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});
  EXPECT_FALSE(fullSearch.searchSpaceEmpty());
  EXPECT_FALSE(fullSearch.searchSpaceOneOption());
  EXPECT_THAT(fullSearch.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

TEST_F(FullSearchTest, testRemoveN3OptionRemoveAll) {
  autopas::FullSearch fullSearch(
      {autopas::ContainerOption::linkedCells}, {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos}, {autopas::Newton3Option::enabled});

  EXPECT_THROW(fullSearch.removeN3Option(autopas::Newton3Option::enabled),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_F(FullSearchTest, testRemoveN3OptionRemoveSome) {
  autopas::FullSearch fullSearch({autopas::ContainerOption::linkedCells},
                                 {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
                                 {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
                                 {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});

  EXPECT_NO_THROW(fullSearch.removeN3Option(autopas::Newton3Option::enabled));
  EXPECT_FALSE(fullSearch.searchSpaceEmpty());
  EXPECT_FALSE(fullSearch.searchSpaceOneOption());
}

TEST_F(FullSearchTest, testTune) {
  autopas::FullSearch fullSearch(
      {autopas::ContainerOption::linkedCells},
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, autopas::TraversalOption::c08,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(10);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, autopas::TraversalOption::sliced,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(1);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, autopas::TraversalOption::c01,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(20);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, autopas::TraversalOption::sliced,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
}