/**
 * @file ReinforcementLearningTest.cpp
 * @author P. Metscher
 * @date 24.06.25
 */

#include "ReinforcementLearningTest.h"

#include "autopas/tuning/utils/SearchSpaceGenerators.h"

/**
 * Test if invalid configuration counts are handled correctly by the constructor.
 */
TEST_F(ReinforcementLearningTest, invalidConfigurationConstructor) {
  // Set the exception policy to throw exceptions
  autopas::utils::ExceptionHandler::setBehavior(autopas::utils::ExceptionBehavior::throwException);

  const auto ns = autopas::NumberSetFinite<double>{1.0};
  // Create a search space with valid configurations
  auto searchSpace = configQueueBase;

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
  const auto searchSpace = configQueueBase;

  autopas::ReinforcementLearning rl(searchSpace, 0.8, 0.8);

  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
  autopas::EvidenceCollection evidenceCollection{};

  // Test the initial full search
  std::vector<autopas::Configuration> configQueueResetCopy = configQueue;
  EXPECT_TRUE(rl.reset(0, 0, configQueue, evidenceCollection))
      << "ReinforcementLearning should allow empty queues in the initial full search.";

  EXPECT_EQ(configQueueResetCopy, configQueue)
      << "ReinforcementLearning should not change the config queue during the full search.";

  const autopas::Evidence evidence{0, 0, 100};
  while (true) {
    const std::vector<autopas::Configuration> configQueueCopy = configQueue;
    bool ret = rl.optimizeSuggestions(configQueue, evidenceCollection);

    EXPECT_EQ(configQueue.size(), configQueueCopy.size())
        << "ReinforcementLearning should not change the size of the config queue.";

    EXPECT_EQ(configQueue, configQueueCopy)
        << "ReinforcementLearning should not change the config queue during the full search.";

    if (configQueue.empty()) {
      EXPECT_TRUE(ret) << "ReinforcementLearning should return true when the config queue is empty.";
      break;
    }

    rl.addEvidence(configQueue.back(), evidence);
    configQueue.pop_back();
  }
}

/**
 * Test for selecting the optimal configuration on the second tuning phase with switching configurations.
 *
 * If the best configuration changes from the first to the second tuning phase, the RL strategy should still find this
 * change, as it chooses all configurations as random exploration samples.
 */
TEST_F(ReinforcementLearningTest, correctConfigExplorationPhase) {
  // Create the core variables
  const auto searchSpace = configQueueBase;

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
    rl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);
    configQueue.pop_back();

    rl.optimizeSuggestions(configQueue, evidenceCollection);
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
  EXPECT_FALSE(rl.reset(1, searchSpace.size() * 2, configQueue, evidenceCollection))
      << "ReinforcementLearning should not clear the queue ever intentionally when resetting.";

  EXPECT_EQ(configQueue.size(), searchSpace.size())
      << "ReinforcementLearning should choose the exploration samples during the reset.";

  autopas::Configuration lastConfig;
  size_t it = 0;

  it++;
  rl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);

  lastConfig = configQueue.back();
  configQueue.pop_back();

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
  const auto searchSpace = configQueueBase;

  autopas::ReinforcementLearning rl(searchSpace, 0.8, 0.8, 3);
  ASSERT_GT(searchSpace.size(), 5);

  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
  autopas::EvidenceCollection evidenceCollection{};
  autopas::Evidence evidence{0, 0, 100};

  // Perform an initial full search
  rl.reset(0, 0, configQueue, evidenceCollection);
  rl.optimizeSuggestions(configQueue, evidenceCollection);

  while (!configQueue.empty()) {
    rl.addEvidence(configQueue.back(), evidence);
    configQueue.pop_back();
    rl.optimizeSuggestions(configQueue, evidenceCollection);
  }

  // Test the second search
  size_t it = 0;
  std::copy(searchSpace.begin(), searchSpace.end(), std::back_inserter(configQueue));
  evidence.tuningPhase = 1;

  EXPECT_FALSE(rl.reset(1, searchSpace.size() * 2, configQueue, evidenceCollection))
      << "ReinforcementLearning should not clear the queue ever intentionally when resetting.";

  EXPECT_EQ(configQueue.size(), 4)
      << "ReinforcementLearning should not clear the queue ever intentionally when resetting.";

  it++;
  rl.addEvidence(configQueue.back(), evidence);
  configQueue.pop_back();

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

  EXPECT_FALSE(rl.reset(1, searchSpace.size() * 2, configQueue, evidenceCollection))
      << "ReinforcementLearning should not clear the queue ever intentionally when resetting.";

  EXPECT_EQ(configQueue.size(), 4)
      << "ReinforcementLearning should never clear the queue intentionally when resetting.";

  it = 1;
  rl.addEvidence(configQueue.back(), evidence);
  configQueue.pop_back();

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

/**
 * Test that if there is a single best configuration (i.e., the configuration with the lowest evidence value)
 * throughout all tuning phases, the ReinforcementLearning strategy always includes this configuration in its
 * sampling set during each exploration phase, ensuring it is never omitted from evaluation.
 */
TEST_F(ReinforcementLearningTest, singleBestConfigAlwaysSampled) {
  // Create the core variables
  const auto searchSpace = configQueueBase;

  autopas::ReinforcementLearning rl(searchSpace, 0.8, 0.8, 3);
  ASSERT_GT(searchSpace.size(), 5);

  autopas::Configuration targetConfig = *(searchSpace.begin()++);
  std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
  autopas::EvidenceCollection evidenceCollection{};
  std::unordered_map<autopas::Configuration, autopas::Evidence, autopas::ConfigHash> evidenceMap;
  for (const auto &config : configQueue) {
    evidenceMap[config] = autopas::Evidence{0, 0, config == targetConfig ? 100 : 110};
  }

  // Perform an initial full search
  rl.reset(0, 0, configQueue, evidenceCollection);
  rl.optimizeSuggestions(configQueue, evidenceCollection);

  while (!configQueue.empty()) {
    rl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);
    configQueue.pop_back();
    rl.optimizeSuggestions(configQueue, evidenceCollection);
  }

  // Perform the other searches
  for (size_t i = 0; i < 200; i++) {
    std::copy(searchSpace.begin(), searchSpace.end(), std::back_inserter(configQueue));
    int occurrences = 0;

    EXPECT_FALSE(rl.reset(1, searchSpace.size() * 2, configQueue, evidenceCollection))
        << "ReinforcementLearning should not clear the queue ever intentionally when resetting.";

    EXPECT_EQ(configQueue.size(), 4) << "ReinforcementLearning should choose the exploration samples during the reset.";

    autopas::Configuration lastConfig;
    rl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);

    if (configQueue.back() == targetConfig) {
      occurrences++;
    }

    lastConfig = configQueue.back();
    configQueue.pop_back();

    while (true) {
      bool ret = rl.optimizeSuggestions(configQueue, evidenceCollection);

      if (configQueue.empty()) {
        EXPECT_TRUE(ret) << "ReinforcementLearning should return true when the config queue is empty.";
        break;
      }

      if (configQueue.back() == targetConfig) {
        occurrences++;
      }

      rl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);

      lastConfig = configQueue.back();
      configQueue.pop_back();
    }

    EXPECT_GE(occurrences, 1) << "ReinforcementLearning should always sample the best configuration.";
  }
}
