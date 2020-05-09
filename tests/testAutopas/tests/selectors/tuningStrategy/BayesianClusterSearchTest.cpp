/**
 * @file BayesianClusterSearchTest.cpp
 * @author Jan Nguyen
 * @date 12.05.20
 */

#include "BayesianClusterSearchTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

TEST_F(BayesianClusterSearchTest, testMaxEvidence) {
  size_t maxEvidence = 4;
  autopas::BayesianClusterSearch bayesClusterSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1}),
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01, autopas::TraversalOption::sliced},
      {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled}, maxEvidence);

  // while #evidence < maxEvidence. tuning -> True
  for (size_t i = 1; i < maxEvidence; ++i) {
    bayesClusterSearch.addEvidence(i);
    EXPECT_TRUE(bayesClusterSearch.tune());
  }

  // #evidence == maxEvidence. tuning -> False
  bayesClusterSearch.addEvidence(-1);
  EXPECT_FALSE(bayesClusterSearch.tune());
}

TEST_F(BayesianClusterSearchTest, testFindBest) {
  size_t maxEvidence = 7;
  unsigned long seed = 21;
  autopas::BayesianClusterSearch bayesClusterSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1, 2}),
      {autopas::TraversalOption::c08, autopas::TraversalOption::c01},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
      {autopas::Newton3Option::disabled, autopas::Newton3Option::enabled}, maxEvidence,
      autopas::AcquisitionFunctionOption::upperConfidenceBound, 50, seed);

  // configuration to find
  autopas::FeatureVector best(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::c08,
                              autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled);

  while (bayesClusterSearch.tune()) {
    autopas::FeatureVector current(bayesClusterSearch.getCurrentConfiguration());

    Eigen::VectorXd diff = best - current;
    double distanceSquared = diff.array().square().sum();
    long dummyTime = static_cast<long>(654321 * distanceSquared);

    bayesClusterSearch.addEvidence(dummyTime);
  }

  autopas::FeatureVector prediction(bayesClusterSearch.getCurrentConfiguration());
  EXPECT_EQ(prediction, best);
}
