/**
 * @file RandomSearchTest.cpp
 * @author Jan Nguyen
 * @date 10.07.19
 */

#include "RandomSearchTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_F(RandomSearchTest, testSearchSpaceEmpty) {
  autopas::RandomSearch randomSearch(std::set<autopas::ContainerOption>({}));
  EXPECT_TRUE(randomSearch.searchSpaceIsEmpty());
  EXPECT_FALSE(randomSearch.searchSpaceIsTrivial());
  EXPECT_THAT(randomSearch.getAllowedContainerOptions(), ::testing::IsEmpty());
}

TEST_F(RandomSearchTest, testSearchSpaceOneOption) {
  autopas::RandomSearch randomSearch({autopas::ContainerOption::directSum}, autopas::NumberSetFinite<double>({1.}),
                                     {autopas::TraversalOption::directSumTraversal},
                                     {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa},
                                     {autopas::Newton3Option::enabled});
  EXPECT_FALSE(randomSearch.searchSpaceIsEmpty());
  EXPECT_TRUE(randomSearch.searchSpaceIsTrivial());
  EXPECT_THAT(randomSearch.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_F(RandomSearchTest, testSearchSpaceMoreOptions) {
  autopas::RandomSearch randomSearch({autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1.}),
                                     {autopas::TraversalOption::c08}, {autopas::LoadEstimatorOption::none},
                                     {autopas::DataLayoutOption::soa},
                                     {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});
  EXPECT_FALSE(randomSearch.searchSpaceIsEmpty());
  EXPECT_FALSE(randomSearch.searchSpaceIsTrivial());
  EXPECT_THAT(randomSearch.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

TEST_F(RandomSearchTest, testRemoveN3OptionRemoveAll) {
  autopas::RandomSearch randomSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1.}),
      {autopas::TraversalOption::c08, autopas::TraversalOption::sliced}, {autopas::LoadEstimatorOption::none},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos}, {autopas::Newton3Option::enabled});

  EXPECT_THROW(randomSearch.removeN3Option(autopas::Newton3Option::enabled),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_F(RandomSearchTest, testRemoveN3OptionRemoveSome) {
  autopas::RandomSearch randomSearch({autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1.}),
                                     {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
                                     {autopas::LoadEstimatorOption::none},
                                     {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
                                     {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});

  EXPECT_NO_THROW(randomSearch.removeN3Option(autopas::Newton3Option::enabled));
  EXPECT_FALSE(randomSearch.searchSpaceIsEmpty());
  EXPECT_FALSE(randomSearch.searchSpaceIsTrivial());
}
