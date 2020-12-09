/**
 * @file BayesianClusterSearchTest.cpp
 * @author Jan Nguyen
 * @date 12.05.20
 */

#include "BayesianClusterSearchTest.h"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

namespace BayesianClusterSearchTest {

TEST_F(BayesianClusterSearchTest, testMaxEvidence) {
  size_t maxEvidence = 4;
  autopas::BayesianClusterSearch bayesClusterSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberSetFinite<double>({1}),
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa}, {autopas::Newton3Option::disabled},
      maxEvidence);

  size_t iteration = 0;

  // fullSearch phase
  do {
    autopas::FeatureVector current(bayesClusterSearch.getCurrentConfiguration());
    bayesClusterSearch.addEvidence(0, iteration);
    ++iteration;
  } while (bayesClusterSearch.tune());

  // artificial time skip
  iteration += 100;
  bayesClusterSearch.reset(iteration);

  // while #evidence < maxEvidence. tuning -> True
  for (size_t i = 1; i < maxEvidence; ++i, ++iteration) {
    bayesClusterSearch.addEvidence(i, iteration);
    EXPECT_TRUE(bayesClusterSearch.tune());
  }

  // #evidence == maxEvidence. tuning -> False
  bayesClusterSearch.addEvidence(-1, maxEvidence);
  EXPECT_FALSE(bayesClusterSearch.tune());
}

/**
 * Find best configuration if configuration are similar through tuning phases.
 */
TEST_F(BayesianClusterSearchTest, testFindBestSimilar) {
  constexpr size_t maxEvidence = 10;
  constexpr unsigned long seed = 21;
  // we use a dummy time function which increases linearly with the squared distance to the chosen optimal solution
  constexpr double timePerDistanceSquared = 654321;
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
    return static_cast<long>(timePerDistanceSquared * distanceSquared);
  };

  size_t iteration = 0;

  // first tuning phase
  do {
    autopas::FeatureVector current(bayesClusterSearch.getCurrentConfiguration());
    bayesClusterSearch.addEvidence(dummyTimeFun(current), iteration);
    ++iteration;
  } while (bayesClusterSearch.tune());

  // artificial time skip
  iteration += 50;
  bayesClusterSearch.reset(iteration);

  // second tuning phase
  do {
    autopas::FeatureVector current(bayesClusterSearch.getCurrentConfiguration());
    bayesClusterSearch.addEvidence(dummyTimeFun(current), iteration);
    ++iteration;
  } while (bayesClusterSearch.tune());

  autopas::FeatureVector prediction(bayesClusterSearch.getCurrentConfiguration());
  long predictionTime = dummyTimeFun(prediction);
  long bestTime = dummyTimeFun(best);
  EXPECT_NEAR(predictionTime, bestTime, timePerDistanceSquared * 0.1);
}

/**
 * Find best configuration if optimal configuration changes over time.
 */
TEST_F(BayesianClusterSearchTest, testFindBestDifferent) {
  const size_t maxEvidence = 15;
  const unsigned long seed = 21;
  // we use a dummy time function which increases linearly with the squared distance to the chosen optimal solution
  constexpr double timePerDistanceSquared = 654321;

  autopas::BayesianClusterSearch bayesClusterSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberInterval<double>(1, 2),
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01}, {autopas::LoadEstimatorOption::none},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
      {autopas::Newton3Option::disabled, autopas::Newton3Option::enabled}, maxEvidence,
      autopas::AcquisitionFunctionOption::upperConfidenceBound, 50, seed);

  // optimal configuration in first tuning phase
  autopas::FeatureVector best1(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                               autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                               autopas::Newton3Option::enabled);

  auto dummyTimeFun1 = [&best1](autopas::FeatureVector target) -> long {
    Eigen::VectorXd diff = best1 - target;
    double distanceSquared = diff.squaredNorm();
    return static_cast<long>(timePerDistanceSquared * distanceSquared);
  };

  // optimal configuration in second tuning phase (only traversal changed)
  autopas::FeatureVector best2(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01,
                               autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                               autopas::Newton3Option::enabled);

  auto dummyTimeFun2 = [&best2](autopas::FeatureVector target) -> long {
    Eigen::VectorXd diff = best2 - target;
    double distanceSquared = diff.squaredNorm();
    return static_cast<long>(timePerDistanceSquared * distanceSquared);
  };

  size_t iteration = 0;

  // first tuning phase
  do {
    autopas::FeatureVector current(bayesClusterSearch.getCurrentConfiguration());
    bayesClusterSearch.addEvidence(dummyTimeFun1(current), iteration);
    ++iteration;
  } while (bayesClusterSearch.tune());

  // artificial time skip
  iteration += 50;
  bayesClusterSearch.reset(iteration);

  // second tuning phase
  do {
    autopas::FeatureVector current(bayesClusterSearch.getCurrentConfiguration());
    bayesClusterSearch.addEvidence(dummyTimeFun2(current), iteration);
    ++iteration;
  } while (bayesClusterSearch.tune());

  autopas::FeatureVector prediction(bayesClusterSearch.getCurrentConfiguration());
  long predictionTime = dummyTimeFun2(prediction);
  long bestTime = dummyTimeFun2(best2);
  EXPECT_NEAR(predictionTime, bestTime, timePerDistanceSquared * 0.25);
}

/**
 * Find best configuration if optimal configuration changes drastically over time.
 */
TEST_F(BayesianClusterSearchTest, testFindBestVeryDifferent) {
  const size_t maxEvidence = 20;
  const unsigned long seed = 21;
  // we use a dummy time function which increases linearly with the squared distance to the chosen optimal solution
  constexpr double timePerDistanceSquared = 654321;

  autopas::BayesianClusterSearch bayesClusterSearch(
      {autopas::ContainerOption::linkedCells}, autopas::NumberInterval<double>(1, 2),
      {autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01}, {autopas::LoadEstimatorOption::none},
      {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
      {autopas::Newton3Option::disabled, autopas::Newton3Option::enabled}, maxEvidence,
      autopas::AcquisitionFunctionOption::upperConfidenceBound, 50, seed);

  // optimal configuration in first tuning phase
  autopas::FeatureVector best1(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                               autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                               autopas::Newton3Option::enabled);

  auto dummyTimeFun1 = [&best1](autopas::FeatureVector target) -> long {
    Eigen::VectorXd diff = best1 - target;
    double distanceSquared = diff.squaredNorm();
    return static_cast<long>(timePerDistanceSquared * distanceSquared);
  };

  // optimal configuration in second tuning phase (every option changed)
  autopas::FeatureVector best2(autopas::ContainerOption::linkedCells, 2., autopas::TraversalOption::lc_c01,
                               autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                               autopas::Newton3Option::disabled);

  auto dummyTimeFun2 = [&best2](autopas::FeatureVector target) -> long {
    Eigen::VectorXd diff = best2 - target;
    double distanceSquared = diff.squaredNorm();
    return static_cast<long>(timePerDistanceSquared * distanceSquared);
  };

  size_t iteration = 0;

  // first tuning phase
  do {
    autopas::FeatureVector current(bayesClusterSearch.getCurrentConfiguration());
    bayesClusterSearch.addEvidence(dummyTimeFun1(current), iteration);
    ++iteration;
  } while (bayesClusterSearch.tune());

  // artificial time skip
  iteration += 50;
  bayesClusterSearch.reset(iteration);

  // second tuning phase
  do {
    autopas::FeatureVector current(bayesClusterSearch.getCurrentConfiguration());
    bayesClusterSearch.addEvidence(dummyTimeFun2(current), iteration);
    ++iteration;
  } while (bayesClusterSearch.tune());

  autopas::FeatureVector prediction(bayesClusterSearch.getCurrentConfiguration());
  long predictionTime = dummyTimeFun2(prediction);
  long bestTime = dummyTimeFun2(best2);
  EXPECT_NEAR(predictionTime, bestTime, timePerDistanceSquared * 0.25);
}

} // end namespace BayesianClusterSearchTest
