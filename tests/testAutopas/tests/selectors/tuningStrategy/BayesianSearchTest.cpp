/**
 * @file BayesianSearchTest.cpp
 * @author Jan Nguyen
 * @date 12.06.19
 */

#include "BayesianSearchTest.h"
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_F(BayesianSearchTest, testSearchSpaceEmpty) {
  autopas::BayesianSearch bayesianSearch(std::set<autopas::ContainerOption>({}));
  EXPECT_TRUE(bayesianSearch.searchSpaceIsEmpty());
  EXPECT_FALSE(bayesianSearch.searchSpaceIsTrivial());
  EXPECT_THAT(bayesianSearch.getAllowedContainerOptions(), ::testing::IsEmpty());
}

TEST_F(BayesianSearchTest, testSearchSpaceOneOption) {
  autopas::BayesianSearch bayesianSearch({autopas::ContainerOption::directSum}, autopas::NumberSetFinite<double>({1.}),
                                         {autopas::TraversalOption::directSumTraversal},
                                         {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::enabled});
  EXPECT_FALSE(bayesianSearch.searchSpaceIsEmpty());
  EXPECT_TRUE(bayesianSearch.searchSpaceIsTrivial());
  EXPECT_THAT(bayesianSearch.getAllowedContainerOptions(), ::testing::ElementsAre(autopas::ContainerOption::directSum));
}

TEST_F(BayesianSearchTest, testSearchSpaceMoreOptions) {
  autopas::BayesianSearch bayesianSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1.}), {autopas::TraversalOption::c08},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});
  EXPECT_FALSE(bayesianSearch.searchSpaceIsEmpty());
  EXPECT_FALSE(bayesianSearch.searchSpaceIsTrivial());
  EXPECT_THAT(bayesianSearch.getAllowedContainerOptions(),
              ::testing::ElementsAre(autopas::ContainerOption::linkedCells));
}

TEST_F(BayesianSearchTest, testRemoveN3OptionRemoveAll) {
  autopas::BayesianSearch bayesianSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1.}),
      {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos}, {autopas::Newton3Option::enabled});

  EXPECT_THROW(bayesianSearch.removeN3Option(autopas::Newton3Option::enabled),
               autopas::utils::ExceptionHandler::AutoPasException);
}

TEST_F(BayesianSearchTest, testRemoveN3OptionRemoveSome) {
  autopas::BayesianSearch bayesianSearch({autopas::ContainerOption::linkedCells},
                                         autopas::NumberSetFinite<double>({1.}),
                                         {autopas::TraversalOption::c08, autopas::TraversalOption::sliced},
                                         {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
                                         {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled});

  EXPECT_NO_THROW(bayesianSearch.removeN3Option(autopas::Newton3Option::enabled));
  EXPECT_FALSE(bayesianSearch.searchSpaceIsEmpty());
  EXPECT_FALSE(bayesianSearch.searchSpaceIsTrivial());
}
