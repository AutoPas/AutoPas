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

TEST_F(FullSearchTest, testTune) {
  autopas::FullSearch fullSearch(
      {autopas::ContainerOption::linkedCells}, {1.},
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(10);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(1);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c01,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(20);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::sliced,
                                   autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
}
