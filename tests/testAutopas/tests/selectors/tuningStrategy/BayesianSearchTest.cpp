/**
 * @file BayesianSearchTest.cpp
 * @author Jan Nguyen
 * @date 12.06.19
 */

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "BayesianSearchTest.h"
#include "autopas/options/Newton3Option.h"

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

TEST_F(BayesianSearchTest, testMaxEvidence) {
  size_t maxEvidence = 4;
  autopas::BayesianSearch bayesSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1}),
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, maxEvidence);

  // while #evidence < maxEvidence. tuning -> True
  for (size_t i = 1; i < maxEvidence; ++i) {
    bayesSearch.addEvidence(i);
    EXPECT_TRUE(bayesSearch.tune());
  }

  // #evidence == maxEvidence. tuning -> False
  bayesSearch.addEvidence(-1);
  EXPECT_FALSE(bayesSearch.tune());
}

TEST_F(BayesianSearchTest, testFindBest) {
  size_t maxEvidence = 16;
  unsigned long seed = 21;
  autopas::BayesianSearch bayesSearch({autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1, 2}),
                                      {autopas::TraversalOption::c08, autopas::TraversalOption::c01},
                                      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
                                      {autopas::Newton3Option::disabled, autopas::Newton3Option::enabled}, maxEvidence,
                                      autopas::AcquisitionFunctionOption::lowerConfidenceBound, 1000, seed);

  // configuration to find
  autopas::FeatureVector best(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                              autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled);

  while (bayesSearch.tune()) {
    autopas::FeatureVector current(bayesSearch.getCurrentConfiguration());

    Eigen::VectorXd diff = best - current;
    double distanceSquared = diff.array().square().sum();
    long dummyTime = static_cast<long>(654321 * distanceSquared);

    bayesSearch.addEvidence(dummyTime);
  }

  autopas::FeatureVector prediction(bayesSearch.getCurrentConfiguration());
  EXPECT_EQ(prediction, best);
}
