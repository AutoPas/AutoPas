/**
 * @file BayesianSearchTest.cpp
 * @author Jan Nguyen
 * @date 12.06.19
 */

#include "BayesianSearchTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

namespace BayesianSearchTest {

TEST_F(BayesianSearchTest, testMaxEvidence) {
  size_t maxEvidence = 4;
  autopas::BayesianSearch bayesSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1}),
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      maxEvidence);

  // while #evidence < maxEvidence. tuning -> True
  for (size_t i = 1; i < maxEvidence; ++i) {
    bayesSearch.addEvidence(i, 0);
    EXPECT_TRUE(bayesSearch.tune());
  }

  // #evidence == maxEvidence. tuning -> False
  bayesSearch.addEvidence(-1, 0);
  EXPECT_FALSE(bayesSearch.tune());
}

TEST_F(BayesianSearchTest, testFindBest) {
  size_t maxEvidence = 7;
  unsigned long seed = 21;
  autopas::BayesianSearch bayesSearch({autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1, 2}),
                                      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01},
                                      {autopas::LoadEstimatorOption::none},
                                      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
                                      {autopas::Newton3Option::disabled, autopas::Newton3Option::enabled}, maxEvidence,
                                      autopas::AcquisitionFunctionOption::upperConfidenceBound, 50, seed);

  // configuration to find
  autopas::FeatureVector best(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                              autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                              autopas::Newton3Option::enabled);

  while (bayesSearch.tune()) {
    autopas::FeatureVector current(bayesSearch.getCurrentConfiguration());

    Eigen::VectorXd diff = best - current;
    double distanceSquared = diff.array().square().sum();
    long dummyTime = static_cast<long>(654321 * distanceSquared);

    bayesSearch.addEvidence(dummyTime, 0);
  }

  autopas::FeatureVector prediction(bayesSearch.getCurrentConfiguration());
  EXPECT_EQ(prediction, best);
}

}  // end namespace BayesianSearchTest
