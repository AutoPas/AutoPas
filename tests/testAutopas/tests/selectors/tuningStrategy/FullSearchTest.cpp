/**
 * @file FullSearchTest.cpp
 * @author F. Gratl
 * @date 6/5/19
 */

#include "FullSearchTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

namespace FullSearchTest {

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

  testedConfigs.emplace_back(fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(10, 0);

  fullSearch.tune();
  testedConfigs.emplace_back(fullSearch.getCurrentConfiguration());
  autopas::Configuration optimalConfig = fullSearch.getCurrentConfiguration();
  fullSearch.addEvidence(1, 0);

  fullSearch.tune();
  testedConfigs.emplace_back(fullSearch.getCurrentConfiguration());
  fullSearch.addEvidence(20, 0);

  EXPECT_THAT(allConfigs, testing::UnorderedElementsAreArray(testedConfigs));

  fullSearch.tune();
  EXPECT_EQ(optimalConfig, fullSearch.getCurrentConfiguration());
}

} // end namespace FullSearchTest
