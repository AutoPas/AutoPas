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
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      maxEvidence);

  // while #evidence < maxEvidence. tuning -> True
  for (size_t i = 1; i < maxEvidence; ++i) {
    bayesClusterSearch.addEvidence(i, 0);
    EXPECT_TRUE(bayesClusterSearch.tune());
  }

  // #evidence == maxEvidence. tuning -> False
  bayesClusterSearch.addEvidence(-1, 0);
  EXPECT_FALSE(bayesClusterSearch.tune());
}

TEST_F(BayesianClusterSearchTest, testFindBest) {
  size_t maxEvidence = 7;
  unsigned long seed = 21;
  autopas::BayesianClusterSearch bayesClusterSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberInterval<double>(1, 2),
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01}, {autopas::LoadEstimatorOption::none},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
      {autopas::Newton3Option::disabled, autopas::Newton3Option::enabled}, maxEvidence,
      autopas::AcquisitionFunctionOption::upperConfidenceBound, 50, seed);

  // configuration to find
  autopas::FeatureVector best(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                              autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                              autopas::Newton3Option::enabled);

  auto dummyTimeFun = [&best](autopas::FeatureVector target) -> long {
    Eigen::VectorXd diff = best - target;
    double distanceSquared = diff.squaredNorm();
    return static_cast<long>(654321 * distanceSquared);
  };

  while (bayesClusterSearch.tune()) {
    autopas::FeatureVector current(bayesClusterSearch.getCurrentConfiguration());
    bayesClusterSearch.addEvidence(dummyTimeFun(current), 0);
  }

  autopas::FeatureVector prediction(bayesClusterSearch.getCurrentConfiguration());
  long predictionTime = dummyTimeFun(prediction);
  long bestTime = dummyTimeFun(best);
  EXPECT_NEAR(predictionTime, bestTime, 1);
}
