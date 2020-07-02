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
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled});

  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                                   autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                                   autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(10, 0);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
                                   autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                                   autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(1, 0);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01,
                                   autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                                   autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(20, 0);

  fullSearch.tune();
  EXPECT_EQ(autopas::Configuration(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_sliced,
                                   autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                                   autopas::Newton3Option::disabled),
            fullSearch.getCurrentConfiguration());
}
