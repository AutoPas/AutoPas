/**
 * @file BayesianClusterSearchTest.cpp
 * @author Jan Nguyen
 * @date 12.05.20
 */

#include "BayesianClusterSearchTest.h"

#include <gmock/gmock-matchers.h>

#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "autopas/utils/ArrayUtils.h"

TEST_F(BayesianClusterSearchTest, testMaxEvidence) {
  const size_t maxEvidence = 3;

  const std::set<autopas::ContainerOption> containerOptions{autopas::ContainerOption::linkedCells};
  const std::set<autopas::TraversalOption> traversalOptions{
      autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced};
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions{autopas::LoadEstimatorOption::none};
  const std::set<autopas::DataLayoutOption> dataLayoutOptions{autopas::DataLayoutOption::aos,
                                                              autopas::DataLayoutOption::soa};
  const std::set<autopas::Newton3Option> newton3Options{autopas::Newton3Option::disabled};
  const autopas::NumberSetFinite<double> cellSizeFactors{1};
  const std::set<autopas::VectorizationPatternOption> vecPatternOptions{autopas::VectorizationPatternOption::p1xVec};

  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      vecPatternOptions, autopas::InteractionTypeOption::pairwise);
  autopas::BayesianClusterSearch bayesClusterSearch(autopas::InteractionTypeOption::pairwise, containerOptions,
                                                    cellSizeFactors, traversalOptions, loadEstimatorOptions,
                                                    dataLayoutOptions, newton3Options, vecPatternOptions, maxEvidence);

  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};

  ASSERT_GT(configQueue.size(), maxEvidence) << "The queue must be longer than the number of evidence that should be "
                                                "collected. Otherwise the cutoff mechanism is not tested.";

  autopas::EvidenceCollection evidenceCollection{};

  // The first tuning phase would have been only full search
  // so skip ahead to the start of the second phase.
  size_t iteration = searchSpace.size() + 100;
  bayesClusterSearch.reset(iteration, 1, configQueue, evidenceCollection);

  // Simulate relevant part of an iteration
  for (size_t i = 0; i < searchSpace.size() and not configQueue.empty(); ++i, ++iteration) {
    ASSERT_LT(i, maxEvidence) << "Tuning phase was not aborted after maxEvidence=" << maxEvidence << " iterations.\n"
                              << "configQueue: " << autopas::utils::ArrayUtils::to_string(configQueue);
    const autopas::Evidence evidence{iteration, 1ul, 42l};
    bayesClusterSearch.addEvidence(configQueue.back(), evidence);
    evidenceCollection.addEvidence(configQueue.back(), evidence);
    configQueue.pop_back();
    ASSERT_FALSE(configQueue.empty()) << "All configurations were popped out of the queue prematurely.";
    bayesClusterSearch.optimizeSuggestions(configQueue, evidenceCollection);
  }
  ASSERT_TRUE(configQueue.empty()) << "BayesianClusterSearch::optimizeSuggestions() should clear the queue once the "
                                      "maximum number of evidence is reached to signal the end of a tuning phase.\n"
                                      "configQueue: "
                                   << autopas::utils::ArrayUtils::to_string(configQueue);
}

/**
 * Find best configuration if configuration are similar through tuning phases.
 */
TEST_F(BayesianClusterSearchTest, testFindBestSimilar) {
  constexpr size_t maxEvidence = 10;
  constexpr unsigned long seed = 21;
  constexpr size_t predNumLHSamples = 50;
  // we use a dummy time function which increases linearly with the squared distance to the chosen optimal solution
  constexpr double timePerDistanceSquared = 654321;

  const std::set<autopas::ContainerOption> containerOptions{autopas::ContainerOption::linkedCells};
  const std::set<autopas::TraversalOption> traversalOptions{
      autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced};
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions{autopas::LoadEstimatorOption::none};
  const std::set<autopas::DataLayoutOption> dataLayoutOptions{autopas::DataLayoutOption::aos,
                                                              autopas::DataLayoutOption::soa};
  const std::set<autopas::Newton3Option> newton3Options{autopas::Newton3Option::disabled,
                                                        autopas::Newton3Option::enabled};
  const autopas::NumberSetFinite<double> cellSizeFactors{1.};
  const std::set<autopas::VectorizationPatternOption> vecPatternOptions{autopas::VectorizationPatternOption::p1xVec};

  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      vecPatternOptions, autopas::InteractionTypeOption::pairwise);
  autopas::BayesianClusterSearch bayesClusterSearch(
      autopas::InteractionTypeOption::pairwise, containerOptions, cellSizeFactors, traversalOptions,
      loadEstimatorOptions, dataLayoutOptions, newton3Options, vecPatternOptions, maxEvidence,
      autopas::AcquisitionFunctionOption::upperConfidenceBound, "", predNumLHSamples, seed);

  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
  autopas::EvidenceCollection evidenceCollection{};

  // configuration to find
  const autopas::FeatureVector best(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                                    autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                                    autopas::Newton3Option::enabled, autopas::InteractionTypeOption::pairwise,
                                    autopas::VectorizationPatternOption::p1xVec);

  auto dummyTimeFun = [&best](autopas::FeatureVector target) -> long {
    const Eigen::VectorXd diff = best - target;
    const double distanceSquared = diff.squaredNorm();
    return static_cast<long>(timePerDistanceSquared * distanceSquared);
  };

  size_t iteration = 0;
  size_t tuningPhase = 0;

  // first tuning phase (full search)
  bayesClusterSearch.reset(iteration, 0, configQueue, evidenceCollection);
  for (const auto conf : searchSpace) {
    const autopas::FeatureVector currentFeatureVec(conf);
    const autopas::Evidence evidence{iteration, tuningPhase, dummyTimeFun(currentFeatureVec)};
    bayesClusterSearch.addEvidence(conf, evidence);
    evidenceCollection.addEvidence(conf, evidence);
    ++iteration;
  }

  // artificial time skip
  iteration += 50;
  ++tuningPhase;
  bayesClusterSearch.reset(iteration, tuningPhase, configQueue, evidenceCollection);

  // second tuning phase
  while (not configQueue.empty()) {
    const autopas::FeatureVector currentFeatureVec(configQueue.back());
    const autopas::Evidence evidence{iteration, tuningPhase, dummyTimeFun(currentFeatureVec)};
    bayesClusterSearch.addEvidence(configQueue.back(), evidence);
    evidenceCollection.addEvidence(configQueue.back(), evidence);
    configQueue.pop_back();
    ++iteration;
    bayesClusterSearch.optimizeSuggestions(configQueue, evidenceCollection);
  }

  const auto [bestConf, bestEvidence] = evidenceCollection.getLatestOptimalConfiguration();
  const autopas::FeatureVector prediction(bestConf);
  const long predictionTime = dummyTimeFun(prediction);
  const long bestTime = dummyTimeFun(best);
  EXPECT_NEAR(predictionTime, bestTime, timePerDistanceSquared * 0.1);
}

/**
 * Find best configuration if optimal configuration changes over time.
 */
TEST_F(BayesianClusterSearchTest, testFindBestDifferent) {
  const size_t maxEvidence = 15;
  const unsigned long seed = 21;
  constexpr size_t predNumLHSamples = 50;
  // we use a dummy time function which increases linearly with the squared distance to the chosen optimal solution
  constexpr double timePerDistanceSquared = 654321;

  const std::set<autopas::ContainerOption> containerOptions{autopas::ContainerOption::linkedCells};
  const std::set<autopas::TraversalOption> traversalOptions{autopas::TraversalOption::lc_c08,
                                                            autopas::TraversalOption::lc_c01};
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions{autopas::LoadEstimatorOption::none};
  const std::set<autopas::DataLayoutOption> dataLayoutOptions{autopas::DataLayoutOption::aos,
                                                              autopas::DataLayoutOption::soa};
  const std::set<autopas::Newton3Option> newton3Options{autopas::Newton3Option::disabled};
  const autopas::NumberSetFinite<double> cellSizeFactors{1., 2.};
  const std::set<autopas::VectorizationPatternOption> vecPatternOptions{autopas::VectorizationPatternOption::p1xVec};

  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      vecPatternOptions, autopas::InteractionTypeOption::pairwise);
  autopas::BayesianClusterSearch bayesClusterSearch(
      autopas::InteractionTypeOption::pairwise, containerOptions, cellSizeFactors, traversalOptions,
      loadEstimatorOptions, dataLayoutOptions, newton3Options, vecPatternOptions, maxEvidence,
      autopas::AcquisitionFunctionOption::upperConfidenceBound, "", predNumLHSamples, seed);

  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
  autopas::EvidenceCollection evidenceCollection{};

  // optimal configuration in first tuning phase
  const autopas::FeatureVector best1(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                                     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                                     autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise,
                                     autopas::VectorizationPatternOption::p1xVec);

  auto dummyTimeFun1 = [&best1](autopas::FeatureVector target) -> long {
    const Eigen::VectorXd diff = best1 - target;
    const double distanceSquared = diff.squaredNorm();
    return static_cast<long>(timePerDistanceSquared * distanceSquared);
  };

  // optimal configuration in second tuning phase (only traversal changed)
  const autopas::FeatureVector best2(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c01,
                                     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                                     autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise,
                                     autopas::VectorizationPatternOption::p1xVec);

  auto dummyTimeFun2 = [&best2](autopas::FeatureVector target) -> long {
    const Eigen::VectorXd diff = best2 - target;
    const double distanceSquared = diff.squaredNorm();
    return static_cast<long>(timePerDistanceSquared * distanceSquared);
  };

  size_t iteration = 0;
  size_t tuningPhase = 0;

  // first tuning phase (full search)
  bayesClusterSearch.reset(iteration, 0, configQueue, evidenceCollection);
  for (const auto &conf : searchSpace) {
    const autopas::FeatureVector currentFeatureVec(conf);
    const autopas::Evidence evidence{iteration, tuningPhase, dummyTimeFun1(currentFeatureVec)};
    bayesClusterSearch.addEvidence(conf, evidence);
    evidenceCollection.addEvidence(conf, evidence);
    ++iteration;
  }

  // artificial time skip
  iteration += 50;
  ++tuningPhase;
  bayesClusterSearch.reset(iteration, tuningPhase, configQueue, evidenceCollection);

  // second tuning phase
  while (not configQueue.empty()) {
    const autopas::FeatureVector currentFeatureVec(configQueue.back());
    const autopas::Evidence evidence{iteration, tuningPhase, dummyTimeFun2(currentFeatureVec)};
    bayesClusterSearch.addEvidence(configQueue.back(), evidence);
    evidenceCollection.addEvidence(configQueue.back(), evidence);
    configQueue.pop_back();
    ++iteration;
    bayesClusterSearch.optimizeSuggestions(configQueue, evidenceCollection);
  }

  const auto [bestConf, bestEvidence] = evidenceCollection.getLatestOptimalConfiguration();
  const autopas::FeatureVector prediction(bestConf);
  const long predictionTime = dummyTimeFun2(prediction);
  const long bestTime = dummyTimeFun2(best2);
  EXPECT_NEAR(predictionTime, bestTime, timePerDistanceSquared * 0.25);
}

/**
 * Find best configuration if optimal configuration changes drastically over time.
 */
TEST_F(BayesianClusterSearchTest, testFindBestVeryDifferent) {
  const size_t maxEvidence = 20;
  const unsigned long seed = 21;
  constexpr size_t predNumLHSamples = 50;
  // we use a dummy time function which increases linearly with the squared distance to the chosen optimal solution
  constexpr double timePerDistanceSquared = 654321;

  const std::set<autopas::ContainerOption> containerOptions{autopas::ContainerOption::linkedCells};
  const std::set<autopas::TraversalOption> traversalOptions{autopas::TraversalOption::lc_c08,
                                                            autopas::TraversalOption::lc_c01};
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions{autopas::LoadEstimatorOption::none};
  const std::set<autopas::DataLayoutOption> dataLayoutOptions{autopas::DataLayoutOption::aos,
                                                              autopas::DataLayoutOption::soa};
  const std::set<autopas::Newton3Option> newton3Options{autopas::Newton3Option::disabled,
                                                        autopas::Newton3Option::enabled};
  const autopas::NumberSetFinite<double> cellSizeFactors{1., 2.};
  const std::set<autopas::VectorizationPatternOption> vecPatternOptions{autopas::VectorizationPatternOption::p1xVec};

  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      vecPatternOptions, autopas::InteractionTypeOption::pairwise);
  autopas::BayesianClusterSearch bayesClusterSearch(
      autopas::InteractionTypeOption::pairwise, containerOptions, cellSizeFactors, traversalOptions,
      loadEstimatorOptions, dataLayoutOptions, newton3Options, vecPatternOptions, maxEvidence,
      autopas::AcquisitionFunctionOption::upperConfidenceBound, "", predNumLHSamples, seed);

  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
  autopas::EvidenceCollection evidenceCollection{};

  // optimal configuration in first tuning phase
  const autopas::FeatureVector best1(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                                     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                                     autopas::Newton3Option::enabled, autopas::InteractionTypeOption::pairwise,
                                     autopas::VectorizationPatternOption::p1xVec);

  auto dummyTimeFun1 = [&best1](autopas::FeatureVector target) -> long {
    const Eigen::VectorXd diff = best1 - target;
    const double distanceSquared = diff.squaredNorm();
    return static_cast<long>(timePerDistanceSquared * distanceSquared);
  };

  // optimal configuration in second tuning phase (every option changed)
  const autopas::FeatureVector best2(autopas::ContainerOption::linkedCells, 2., autopas::TraversalOption::lc_c01,
                                     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                                     autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise,
                                     autopas::VectorizationPatternOption::p1xVec);

  auto dummyTimeFun2 = [&best2](autopas::FeatureVector target) -> long {
    const Eigen::VectorXd diff = best2 - target;
    const double distanceSquared = diff.squaredNorm();
    return static_cast<long>(timePerDistanceSquared * distanceSquared);
  };

  size_t iteration = 0;
  size_t tuningPhase = 0;

  // first tuning phase (full search)
  bayesClusterSearch.reset(iteration, 0, configQueue, evidenceCollection);
  for (const auto conf : searchSpace) {
    const autopas::FeatureVector currentFeatureVec(conf);
    const autopas::Evidence evidence{iteration, tuningPhase, dummyTimeFun1(currentFeatureVec)};
    bayesClusterSearch.addEvidence(conf, evidence);
    evidenceCollection.addEvidence(conf, evidence);
    ++iteration;
  }

  // artificial time skip
  iteration += 50;
  ++tuningPhase;
  bayesClusterSearch.reset(iteration, tuningPhase, configQueue, evidenceCollection);

  // second tuning phase
  while (not configQueue.empty()) {
    const autopas::FeatureVector currentFeatureVec(configQueue.back());
    const autopas::Evidence evidence{iteration, tuningPhase, dummyTimeFun2(currentFeatureVec)};
    bayesClusterSearch.addEvidence(configQueue.back(), evidence);
    evidenceCollection.addEvidence(configQueue.back(), evidence);
    configQueue.pop_back();
    ++iteration;
    bayesClusterSearch.optimizeSuggestions(configQueue, evidenceCollection);
  }

  const auto [bestConf, bestEvidence] = evidenceCollection.getLatestOptimalConfiguration();
  const autopas::FeatureVector prediction(bestConf);
  const long predictionTime = dummyTimeFun2(prediction);
  const long bestTime = dummyTimeFun2(best2);
  EXPECT_NEAR(predictionTime, bestTime, timePerDistanceSquared * 0.25);
}
