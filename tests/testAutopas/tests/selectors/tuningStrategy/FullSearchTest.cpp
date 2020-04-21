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
  EXPECT_TRUE(fullSearch.searchSpaceIsEmpty());
  EXPECT_FALSE(fullSearch.searchSpaceIsTrivial());
  EXPECT_THAT(fullSearch.getAllowedContainerOptions(), ::testing::IsEmpty());
}

TEST_F(FullSearchTest, testSearchSpaceOneOption) {
  autopas::FullSearch fullSearch(
      {autopas::Configuration(autopas::ContainerOption::directSum, 1., autopas::TraversalOption::directSumTraversal,
                              autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled)});
  EXPECT_FALSE(fullSearch.searchSpaceIsEmpty());
  EXPECT_TRUE(fullSearch.searchSpaceIsTrivial());
  EXPECT_THAT(fullSearch.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_F(FullSearchTest, testSearchSpaceMoreOptions) {
  autopas::FullSearch fullSearch({autopas::ContainerOption::linkedCells}, {1.}, {autopas::TraversalOption::c08},
                                 {autopas::DataLayoutOption::soa},
                                 {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});
  EXPECT_FALSE(fullSearch.searchSpaceIsEmpty());
  EXPECT_FALSE(fullSearch.searchSpaceIsTrivial());
  EXPECT_THAT(fullSearch.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

TEST_F(FullSearchTest, testRemoveN3OptionRemoveAll) {
  autopas::FullSearch fullSearch(
      {autopas::ContainerOption::linkedCells}, {1.}, {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos}, {autopas::Newton3Option::enabled});

  EXPECT_THROW(fullSearch.removeN3Option(autopas::Newton3Option::enabled),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_F(FullSearchTest, testRemoveN3OptionRemoveSome) {
  autopas::FullSearch fullSearch({autopas::ContainerOption::linkedCells}, {1.},
                                 {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
                                 {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
                                 {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});

  EXPECT_NO_THROW(fullSearch.removeN3Option(autopas::Newton3Option::enabled));
  EXPECT_FALSE(fullSearch.searchSpaceIsEmpty());
  EXPECT_FALSE(fullSearch.searchSpaceIsTrivial());
}

TEST_F(FullSearchTest, testTune) {
  autopas::FullSearch fullSearch(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(10, 0);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(1, 0);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c01,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(20, 0);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
}
