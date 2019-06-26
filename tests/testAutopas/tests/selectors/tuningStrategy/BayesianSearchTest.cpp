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

TEST_F(BayesianSearchTest, testMaxEvidences) {
  size_t maxEvidences = 4;
  autopas::BayesianSearch bayesSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1}),
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, maxEvidences);

  // while #evidences < maxEvidences. tuning -> True
  for (size_t i = 1; i < maxEvidences; ++i) {
    bayesSearch.addEvidence(i);
    EXPECT_TRUE(bayesSearch.tune());
  }

  // #evidences == maxEvidences. tuning -> False
  bayesSearch.addEvidence(-1);
  EXPECT_FALSE(bayesSearch.tune());
}

TEST_F(BayesianSearchTest, testFindBest) {
  size_t maxEvidences = 4;
  autopas::BayesianSearch bayesSearch({autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1, 2}),
                                      {autopas::TraversalOption::c08, autopas::TraversalOption::c01},
                                      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
                                      {autopas::Newton3Option::disabled, autopas::Newton3Option::enabled},
                                      maxEvidences);

  // configuration to find
  autopas::FeatureVector best(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                              autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled);

  while (bayesSearch.tune()) {
    autopas::FeatureVector current(bayesSearch.getCurrentConfiguration());

    Eigen::VectorXd diff = static_cast<Eigen::VectorXd>(best - current);
    double distanceSquared = diff.array().square().sum();
    long dummyTime = static_cast<long>(654321 * distanceSquared);

    bayesSearch.addEvidence(dummyTime);
  }

  autopas::FeatureVector prediction(bayesSearch.getCurrentConfiguration());
  EXPECT_EQ(prediction, best);
}
