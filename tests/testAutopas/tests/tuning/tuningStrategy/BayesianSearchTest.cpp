/**
 * @file BayesianSearchTest.cpp
 * @author Jan Nguyen
 * @date 12.06.19
 */

#include "BayesianSearchTest.h"

#include <gmock/gmock-matchers.h>

#include "autopas/tuning/utils/SearchSpaceGenerators.h"

TEST_F(BayesianSearchTest, testMaxEvidence) {
  const size_t maxEvidence = 3;

  const std::set<autopas::ContainerOption> containerOptions{autopas::ContainerOption::linkedCells};
  const std::set<autopas::TraversalOption> traversalOptions{
      autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced};
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions{autopas::LoadEstimatorOption::none};
  const std::set<autopas::DataLayoutOption> dataLayoutOptions{autopas::DataLayoutOption::aos,
                                                              autopas::DataLayoutOption::soa};
  const std::set<autopas::Newton3Option> newton3Options{autopas::Newton3Option::disabled};
  const autopas::NumberSetFinite<double> cellSizeFactors{1};
  const autopas::NumberSetFinite<int> threadCounts{autopas::Configuration::ThreadCountNoTuning};

  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      autopas::InteractionTypeOption::pairwise, &threadCounts);
  autopas::BayesianSearch bayesSearch(autopas::InteractionTypeOption::pairwise, containerOptions, cellSizeFactors,
                                      traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, threadCounts,
                                      maxEvidence);

  std::vector<autopas::Configuration> configQueue{searchSpace.begin(), searchSpace.end()};
  autopas::EvidenceCollection evidenceCollection{};

  size_t iteration = searchSpace.size() + 100;
  bayesSearch.reset(iteration, 1, configQueue, evidenceCollection);

  for (size_t i = 0; i < searchSpace.size() and not configQueue.empty(); ++i, ++iteration) {
    ASSERT_LT(i, maxEvidence) << "Tuning phase was not aborted after maxEvidence=" << maxEvidence << " iterations.\n"
                              << "configQueue: " << autopas::utils::ArrayUtils::to_string(configQueue);
    const autopas::Evidence evidence{iteration, 1ul, 42l};
    bayesSearch.addEvidence(configQueue.back(), evidence);
    evidenceCollection.addEvidence(configQueue.back(), evidence);
    configQueue.pop_back();
    ASSERT_FALSE(configQueue.empty()) << "All configurations were popped out of the queue prematurely.";
    bayesSearch.optimizeSuggestions(configQueue, evidenceCollection);
  }
  ASSERT_TRUE(configQueue.empty()) << "BayesianSearch::optimizeSuggestions() should clear the queue once the "
                                      "maximum number of evidence is reached to signal the end of a tuning phase.\n"
                                      "configQueue: "
                                   << autopas::utils::ArrayUtils::to_string(configQueue);
}

TEST_F(BayesianSearchTest, testFindBest) {
  constexpr size_t maxEvidence = 7;
  constexpr unsigned long seed = 21;
  constexpr size_t predNumLHSamples = 50;

  const std::set<autopas::ContainerOption> containerOptions{autopas::ContainerOption::linkedCells};
  const std::set<autopas::TraversalOption> traversalOptions{autopas::TraversalOption::lc_c08,
                                                            autopas::TraversalOption::lc_c01};
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions{autopas::LoadEstimatorOption::none};
  const std::set<autopas::DataLayoutOption> dataLayoutOptions{autopas::DataLayoutOption::aos,
                                                              autopas::DataLayoutOption::soa};
  const std::set<autopas::Newton3Option> newton3Options{autopas::Newton3Option::disabled,
                                                        autopas::Newton3Option::enabled};
  const autopas::NumberSetFinite<double> cellSizeFactors{1., 2.};

  const auto threadCounts =
      std::make_unique<autopas::NumberSetFinite<int>>(std::set<int>{autopas::Configuration::ThreadCountNoTuning}).get();
  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      autopas::InteractionTypeOption::pairwise, threadCounts);
  autopas::BayesianSearch bayesSearch(autopas::InteractionTypeOption::pairwise, containerOptions, cellSizeFactors,
                                      traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, *threadCounts,
                                      maxEvidence, autopas::AcquisitionFunctionOption::upperConfidenceBound,
                                      predNumLHSamples, seed);
  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
  autopas::EvidenceCollection evidenceCollection{};

  // configuration to find
  const autopas::FeatureVector best(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                                    autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::soa,
                                    autopas::Newton3Option::enabled, autopas::InteractionTypeOption::pairwise,
                                    threadCounts->getMin());

  // artificial time skip
  size_t iteration = 0;
  size_t tuningPhase = 0;
  bayesSearch.reset(iteration, tuningPhase, configQueue, evidenceCollection);

  // second tuning phase
  while (not configQueue.empty()) {
    const autopas::FeatureVector currentFeatureVec(configQueue.back());
    Eigen::VectorXd diff = best - currentFeatureVec;
    const double distanceSquared = diff.array().square().sum();
    const long dummyTime = static_cast<long>(654321 * distanceSquared);
    const autopas::Evidence evidence{iteration, tuningPhase, dummyTime};
    bayesSearch.addEvidence(configQueue.back(), evidence);
    evidenceCollection.addEvidence(configQueue.back(), evidence);
    configQueue.pop_back();
    ++iteration;
    bayesSearch.optimizeSuggestions(configQueue, evidenceCollection);
  }

  const auto [bestConf, bestEvidence] = evidenceCollection.getLatestOptimalConfiguration();
  autopas::FeatureVector prediction(bestConf);
  EXPECT_EQ(prediction, best);
}
