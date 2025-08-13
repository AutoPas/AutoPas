/**
 * @file ReinforcementLearningTest.cpp
 * @author P. Metscher
 * @date 24.06.25
 */

#include "ReinforcementLearningTest.h"

#include "autopas/tuning/utils/SearchSpaceGenerators.h"

/**
 * Test if invalid configuration counts are handled correctly.
 */
TEST_F(ReinforcementLearningTest, invalidConfigurationCount) {
  // Set the exception policy to throw exceptions
  autopas::utils::ExceptionHandler::setBehavior(autopas::utils::ExceptionBehavior::throwException);

  // Create a search space with only one configuration
  std::set<autopas::Configuration> searchSpace = {};
  EXPECT_THROW(autopas::ReinforcementLearning(searchSpace, 0.8, 0.8),
               autopas::utils::ExceptionHandler::AutoPasException)
      << "The constructor must throw if the search space is empty.";

  // Create a search space with less configurations than the random explorations
  searchSpace = {
      autopas::Configuration(autopas::options::ContainerOption::linkedCells, 1.0,
                             autopas::options::TraversalOption::lc_c01_combined_SoA, autopas::LoadEstimatorOption::none,
                             autopas::options::DataLayoutOption::soa, autopas::Newton3Option::enabled,
                             autopas::options::InteractionTypeOption::pairwise),
      autopas::Configuration(autopas::options::ContainerOption::linkedCells, 1.0,
                             autopas::options::TraversalOption::lc_c01_combined_SoA, autopas::LoadEstimatorOption::none,
                             autopas::options::DataLayoutOption::soa, autopas::Newton3Option::disabled,
                             autopas::options::InteractionTypeOption::pairwise)};
  EXPECT_THROW(autopas::ReinforcementLearning(searchSpace, 0.8, 0.8),
               autopas::utils::ExceptionHandler::AutoPasException)
      << "The constructor must throw if the search space is smaller than the number of exploration samples.";

  const auto ns = autopas::NumberSetFinite<double>{1.0};
  // Create a search space with valid configurations
  searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      {autopas::ContainerOption::linkedCells, autopas::ContainerOption::verletLists,
       autopas::ContainerOption::pairwiseVerletLists},
      {autopas::TraversalOption::lc_c01_combined_SoA, autopas::TraversalOption::lc_c08,
       autopas::TraversalOption::vlp_c08, autopas::TraversalOption::vl_list_iteration},
      {autopas::LoadEstimatorOption::none}, {autopas::DataLayoutOption::soa, autopas::DataLayoutOption::aos},
      {autopas::Newton3Option::enabled, autopas::Newton3Option::disabled}, &ns,
      autopas::InteractionTypeOption::pairwise);

  // Test invalid learning rates
  EXPECT_THROW(autopas::ReinforcementLearning(searchSpace, -0.1, 0.8),
               autopas::utils::ExceptionHandler::AutoPasException)
      << "The constructor must throw when encountering invalid input arguments.";
  EXPECT_THROW(autopas::ReinforcementLearning(searchSpace, 1.1, 0.8),
               autopas::utils::ExceptionHandler::AutoPasException)
      << "The constructor must throw when encountering invalid input arguments.";

  // Test invalid discount factors
  EXPECT_THROW(autopas::ReinforcementLearning(searchSpace, 0.8, -0.1),
               autopas::utils::ExceptionHandler::AutoPasException)
      << "The constructor must throw when encountering invalid input arguments.";
  EXPECT_THROW(autopas::ReinforcementLearning(searchSpace, 0.8, 1.1),
               autopas::utils::ExceptionHandler::AutoPasException)
      << "The constructor must throw when encountering invalid input arguments.";

  // Test multiple valid configurations
  EXPECT_NO_THROW(autopas::ReinforcementLearning(searchSpace, 0.8, 0.8))
      << "The constructor must not fail when passing invalid arguments.";
  EXPECT_NO_THROW(autopas::ReinforcementLearning(searchSpace, 0.5, 0.5))
      << "The constructor must not fail when passing invalid arguments.";
  EXPECT_NO_THROW(autopas::ReinforcementLearning(searchSpace, 0.2, 0.8))
      << "The constructor must not fail when passing invalid arguments.";
  EXPECT_NO_THROW(autopas::ReinforcementLearning(searchSpace, 0.5, 0.1))
      << "The constructor must not fail when passing invalid arguments.";
}

/**
 * Test that the first search is a full search.
 */
TEST_F(ReinforcementLearningTest, firstSearchIsFullSearch) {
  // Create the core variables
  const std::set<autopas::ContainerOption> containerOptions{autopas::ContainerOption::linkedCells,
                                                            autopas::ContainerOption::verletLists,
                                                            autopas::ContainerOption::pairwiseVerletLists};
  const std::set<autopas::TraversalOption> traversalOptions{
      autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced,
      autopas::TraversalOption::vlp_c08, autopas::TraversalOption::vl_list_iteration};
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions{autopas::LoadEstimatorOption::none};
  const std::set<autopas::DataLayoutOption> dataLayoutOptions{autopas::DataLayoutOption::aos,
                                                              autopas::DataLayoutOption::soa};
  const std::set<autopas::Newton3Option> newton3Options{autopas::Newton3Option::disabled,
                                                        autopas::Newton3Option::enabled};
  const autopas::NumberSetFinite<double> cellSizeFactors{1};

  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      autopas::InteractionTypeOption::pairwise);

  autopas::ReinforcementLearning rl(searchSpace, 0.8, 0.8);

  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
  autopas::EvidenceCollection evidenceCollection{};

  // Test the initial full search
  std::vector<autopas::Configuration> configQueueResetCopy = configQueue;
  EXPECT_FALSE(rl.reset(0, 0, configQueue, evidenceCollection))
      << "ReinforcementLearning should not clear the queue ever intentionally when resetting.";

  EXPECT_EQ(configQueueResetCopy, configQueue)
      << "ReinforcementLearning should not change the config queue during the reset.";

  const autopas::Evidence evidence{0, 0, 100};
  while (true) {
    const std::vector<autopas::Configuration> configQueueCopy = configQueue;
    bool ret = rl.optimizeSuggestions(configQueue, evidenceCollection);
    rl.addEvidence(configQueue.back(), evidence);

    EXPECT_EQ(configQueue.size(), configQueueCopy.size())
        << "ReinforcementLearning should not change the size of the config queue.";

    EXPECT_EQ(configQueue, configQueueCopy)
        << "ReinforcementLearning should not change the config queue during the full search.";

    if (configQueue.empty()) {
      EXPECT_TRUE(ret) << "ReinforcementLearning should return true when the config queue is empty.";
      break;
    }

    configQueue.pop_back();
  }
}

/**
 * Test for selecting the optimal configuration on the second tuning phase with switching configurations.
 */
TEST_F(ReinforcementLearningTest, correctConfigExplorationPhase) {
  // Create the core variables
  const std::set<autopas::ContainerOption> containerOptions{autopas::ContainerOption::linkedCells,
                                                            autopas::ContainerOption::verletLists,
                                                            autopas::ContainerOption::pairwiseVerletLists};
  const std::set<autopas::TraversalOption> traversalOptions{
      autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced,
      autopas::TraversalOption::vlp_c08, autopas::TraversalOption::vl_list_iteration};
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions{autopas::LoadEstimatorOption::none};
  const std::set<autopas::DataLayoutOption> dataLayoutOptions{autopas::DataLayoutOption::aos,
                                                              autopas::DataLayoutOption::soa};
  const std::set<autopas::Newton3Option> newton3Options{autopas::Newton3Option::disabled,
                                                        autopas::Newton3Option::enabled};
  const autopas::NumberSetFinite<double> cellSizeFactors{1};

  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      autopas::InteractionTypeOption::pairwise);

  autopas::ReinforcementLearning rl(searchSpace, 0.8, 0.8, searchSpace.size() - 1);

  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
  autopas::EvidenceCollection evidenceCollection{};

  std::unordered_map<autopas::Configuration, autopas::Evidence, autopas::ConfigHash> evidenceMap;
  long i = 100;
  for (const auto &config : configQueue) {
    evidenceMap[config] = autopas::Evidence{0, 0, i};

    i += 20;
  }

  // Perform an initial full search
  rl.reset(0, 0, configQueue, evidenceCollection);

  while (!configQueue.empty()) {
    rl.optimizeSuggestions(configQueue, evidenceCollection);
    rl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);

    configQueue.pop_back();
  }

  // Update the evidence
  std::copy(searchSpace.begin(), searchSpace.end(), std::back_inserter(configQueue));
  autopas::Configuration bestConfig;
  i *= 2;
  for (const auto &config : configQueue) {
    evidenceMap[config] = autopas::Evidence{1, searchSpace.size() * 2, i};
    bestConfig = config;

    i -= 25;
  }

  evidenceMap[bestConfig].value = 1;

  // Test the second search
  std::vector<autopas::Configuration> configQueueResetCopy = configQueue;
  EXPECT_FALSE(rl.reset(1, searchSpace.size() * 2, configQueue, evidenceCollection))
      << "ReinforcementLearning should not clear the queue ever intentionally when resetting.";

  EXPECT_EQ(configQueueResetCopy, configQueue)
      << "ReinforcementLearning should not change the config queue during the reset.";

  autopas::Configuration lastConfig;
  size_t it = 0;
  while (true) {
    bool ret = rl.optimizeSuggestions(configQueue, evidenceCollection);

    if (configQueue.empty()) {
      EXPECT_TRUE(ret) << "ReinforcementLearning should return true when the config queue is empty.";
      break;
    }

    it++;
    rl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);

    lastConfig = configQueue.back();
    configQueue.pop_back();
  }

  EXPECT_EQ(bestConfig, lastConfig)
      << "ReinforcementLearning should do the exploitation search on the best configuration.";
  EXPECT_EQ(it, searchSpace.size() + 1)
      << "ReinforcementLearning should search all configurations once, and the best configuration twice.";
}

/**
 * Test that the correct number of random samples gets sampled.
 */
TEST_F(ReinforcementLearningTest, correctNumberOfRandomSamples) {
  // Create the core variables
  const std::set<autopas::ContainerOption> containerOptions{autopas::ContainerOption::linkedCells,
                                                            autopas::ContainerOption::verletLists,
                                                            autopas::ContainerOption::pairwiseVerletLists};
  const std::set<autopas::TraversalOption> traversalOptions{
      autopas::TraversalOption::lc_c08, autopas::TraversalOption::lc_c01, autopas::TraversalOption::lc_sliced,
      autopas::TraversalOption::vlp_c08, autopas::TraversalOption::vl_list_iteration};
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions{autopas::LoadEstimatorOption::none};
  const std::set<autopas::DataLayoutOption> dataLayoutOptions{autopas::DataLayoutOption::aos,
                                                              autopas::DataLayoutOption::soa};
  const std::set<autopas::Newton3Option> newton3Options{autopas::Newton3Option::disabled,
                                                        autopas::Newton3Option::enabled};
  const autopas::NumberSetFinite<double> cellSizeFactors{0.5, 1, 2.0};

  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      autopas::InteractionTypeOption::pairwise);

  autopas::ReinforcementLearning rl(searchSpace, 0.8, 0.8, 3);
  ASSERT_GT(searchSpace.size(), 5);

  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
  autopas::EvidenceCollection evidenceCollection{};
  autopas::Evidence evidence{0, 0, 100};

  // Perform an initial full search
  rl.reset(0, 0, configQueue, evidenceCollection);

  while (!configQueue.empty()) {
    rl.optimizeSuggestions(configQueue, evidenceCollection);
    rl.addEvidence(configQueue.back(), evidence);
    configQueue.pop_back();
  }

  // Test the second search
  std::copy(searchSpace.begin(), searchSpace.end(), std::back_inserter(configQueue));
  evidence.tuningPhase = 1;
  std::vector<autopas::Configuration> configQueueResetCopy2 = configQueue;
  EXPECT_FALSE(rl.reset(1, searchSpace.size() * 2, configQueue, evidenceCollection))
      << "ReinforcementLearning should not clear the queue ever intentionally when resetting.";

  EXPECT_EQ(configQueueResetCopy2, configQueue)
      << "ReinforcementLearning should not change the config queue during the reset.";

  size_t it = 0;
  while (true) {
    bool ret = rl.optimizeSuggestions(configQueue, evidenceCollection);

    if (configQueue.empty()) {
      EXPECT_TRUE(ret) << "ReinforcementLearning should return true when the config queue is empty.";
      break;
    }

    it++;
    rl.addEvidence(configQueue.back(), evidence);

    configQueue.pop_back();
  }

  EXPECT_EQ(it, 5) << "ReinforcementLearning should search the number of random explorations, the previous best and "
                      "the new best configuration.";

  // Test the third search
  std::copy(searchSpace.begin(), searchSpace.end(), std::back_inserter(configQueue));
  evidence.tuningPhase = 2;
  std::vector<autopas::Configuration> configQueueResetCopy3 = configQueue;
  EXPECT_FALSE(rl.reset(2, searchSpace.size() * 4, configQueue, evidenceCollection))
      << "ReinforcementLearning should not clear the queue ever intentionally when resetting.";

  EXPECT_EQ(configQueueResetCopy3, configQueue)
      << "ReinforcementLearning should not change the config queue during the reset.";

  it = 0;
  while (true) {
    bool ret = rl.optimizeSuggestions(configQueue, evidenceCollection);

    if (configQueue.empty()) {
      EXPECT_TRUE(ret) << "ReinforcementLearning should return true when the config queue is empty.";
      break;
    }

    it++;
    rl.addEvidence(configQueue.back(), evidence);

    configQueue.pop_back();
  }

  EXPECT_EQ(it, 5) << "ReinforcementLearning should search the number of random explorations, the previous best and "
                      "the new best configuration.";
}
