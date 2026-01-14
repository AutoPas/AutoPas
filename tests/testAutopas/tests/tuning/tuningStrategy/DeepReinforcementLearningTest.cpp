/**
 * @file DeepReinforcementLearningTest.cpp
 * @author P. Metscher
 * @date 24.06.25
 */

#include "DeepReinforcementLearningTest.h"

/**
 * Test if valid constructors do not throw.
 */
TEST_F(DeepReinforcementLearningTest, validConstructor) {
  // Set the exception policy to throw exceptions
  autopas::utils::ExceptionHandler::setBehavior(autopas::utils::ExceptionBehavior::throwException);

  EXPECT_NO_THROW(
      autopas::DeepReinforcementLearning(true, 10, autopas::DeepReinforcementLearning::ExplorationMethod::polynomial))
      << "DeepReinforcementLearning valid constructor should not throw.";
  EXPECT_NO_THROW(
      autopas::DeepReinforcementLearning(false, 12, autopas::DeepReinforcementLearning::ExplorationMethod::polynomial))
      << "DeepReinforcementLearning valid constructor should not throw.";
  EXPECT_NO_THROW(
      autopas::DeepReinforcementLearning(true, 9, autopas::DeepReinforcementLearning::ExplorationMethod::random))
      << "DeepReinforcementLearning valid constructor should not throw.";
  EXPECT_NO_THROW(
      autopas::DeepReinforcementLearning(false, 11, autopas::DeepReinforcementLearning::ExplorationMethod::random))
      << "DeepReinforcementLearning valid constructor should not throw.";
  EXPECT_NO_THROW(
      autopas::DeepReinforcementLearning(true, 2, autopas::DeepReinforcementLearning::ExplorationMethod::longestAgo))
      << "DeepReinforcementLearning valid constructor should not throw.";
}

/**
 * Test that invalid constructors do throw.
 */
TEST_F(DeepReinforcementLearningTest, invalidConstructors) {
  // Set the exception policy to throw exceptions
  autopas::utils::ExceptionHandler::setBehavior(autopas::utils::ExceptionBehavior::throwException);

  for (const bool train : {true, false}) {
    for (const autopas::DeepReinforcementLearning::ExplorationMethod explorationMethod :
         {autopas::DeepReinforcementLearning::ExplorationMethod::polynomial,
          autopas::DeepReinforcementLearning::ExplorationMethod::random,
          autopas::DeepReinforcementLearning::ExplorationMethod::longestAgo}) {
      EXPECT_THROW(autopas::DeepReinforcementLearning(train, 0, explorationMethod),
                   autopas::utils::ExceptionHandler::AutoPasException)
          << "DeepReinforcementLearning should reject zero exploration samples.";
      EXPECT_THROW(autopas::DeepReinforcementLearning(train, 1, explorationMethod),
                   autopas::utils::ExceptionHandler::AutoPasException)
          << "DeepReinforcementLearning should reject one exploration sample.";
    }
  }
}

/**
 * Test that the first search is a full search.
 */
TEST_F(DeepReinforcementLearningTest, firstSearchIsFullSearch) {
  // Set the exception policy to throw exceptions
  autopas::utils::ExceptionHandler::setBehavior(autopas::utils::ExceptionBehavior::throwException);

  // Create the core variables
  const auto searchSpace = configQueueBase;

  for (const bool train : {true, false}) {
    for (size_t numExplorationSamples = 2; numExplorationSamples <= 8; numExplorationSamples++) {
      for (const autopas::DeepReinforcementLearning::ExplorationMethod explorationMethod :
           {autopas::DeepReinforcementLearning::ExplorationMethod::polynomial,
            autopas::DeepReinforcementLearning::ExplorationMethod::random,
            autopas::DeepReinforcementLearning::ExplorationMethod::longestAgo}) {
        autopas::DeepReinforcementLearning drl(train, numExplorationSamples, explorationMethod);
        std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
        autopas::EvidenceCollection evidenceCollection{};

        // Test the initial full search
        std::vector<autopas::Configuration> configQueueResetCopy = configQueue;
        EXPECT_FALSE(drl.reset(0, 0, configQueue, evidenceCollection))
            << "DeepReinforcementLearning should not clear the queue ever intentionally when resetting.";

        EXPECT_EQ(configQueueResetCopy, configQueue)
            << "DeepReinforcementLearning should not change the config queue during the reset for the full search.";

        const autopas::Evidence evidence{0, 0, 100};

        drl.addEvidence(configQueue.back(), evidence);
        configQueue.pop_back();

        while (true) {
          const std::vector<autopas::Configuration> configQueueCopy = configQueue;
          bool ret = drl.optimizeSuggestions(configQueue, evidenceCollection);

          EXPECT_EQ(configQueue.size(), configQueueCopy.size())
              << "DeepReinforcementLearning should not change the size of the config queue.";

          EXPECT_EQ(configQueue, configQueueCopy)
              << "DeepReinforcementLearning should not change the config queue during the full search.";

          if (configQueue.empty()) {
            EXPECT_TRUE(ret) << "DeepReinforcementLearning should return true when the config queue is empty.";
            break;
          }

          drl.addEvidence(configQueue.back(), evidence);
          configQueue.pop_back();
        }
      }
    }
  }
}

/**
 * Test that the remaining searches search the correct number of exploration and exploitation samples.
 */
TEST_F(DeepReinforcementLearningTest, validExplorationExploitationSize) {
  // Set the exception policy to throw exceptions
  autopas::utils::ExceptionHandler::setBehavior(autopas::utils::ExceptionBehavior::throwException);

  // Create the core variables
  const auto searchSpace = configQueueBase;

  for (const bool train : {true, false}) {
    for (size_t numExplorationSamples = 2; numExplorationSamples <= 8; numExplorationSamples++) {
      for (const autopas::DeepReinforcementLearning::ExplorationMethod explorationMethod :
           {autopas::DeepReinforcementLearning::ExplorationMethod::polynomial,
            autopas::DeepReinforcementLearning::ExplorationMethod::random,
            autopas::DeepReinforcementLearning::ExplorationMethod::longestAgo}) {
        autopas::DeepReinforcementLearning drl(train, numExplorationSamples, explorationMethod);
        std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
        autopas::EvidenceCollection evidenceCollection{};

        // Perform the initial full search
        drl.reset(0, 0, configQueue, evidenceCollection);
        const autopas::Evidence evidence{0, 0, 100};

        drl.addEvidence(configQueue.back(), evidence);
        configQueue.pop_back();

        while (true) {
          if (configQueue.empty()) {
            break;
          }

          drl.addEvidence(configQueue.back(), evidence);
          configQueue.pop_back();

          const std::vector<autopas::Configuration> configQueueCopy = configQueue;
          drl.optimizeSuggestions(configQueue, evidenceCollection);
        }

        // Perform the future searches
        for (size_t tuningPhase = 1; tuningPhase < 100; tuningPhase++) {
          configQueue = {searchSpace.begin(), searchSpace.end()};

          // Reset
          std::vector<autopas::Configuration> configQueueResetCopy = configQueue;
          EXPECT_FALSE(drl.reset(tuningPhase, tuningPhase * 1000, configQueue, evidenceCollection))
              << "DeepReinforcementLearning should not clear the queue ever intentionally when resetting.";
          EXPECT_EQ(configQueue.size(), numExplorationSamples)
              << "DeepReinforcementLearning should shrink the config queue to the number of exploration samples.";

          drl.addEvidence(configQueue.back(), evidence);
          configQueue.pop_back();

          if (configQueue.empty()) {
            break;
          }

          // Search all exploration samples
          while (true) {
            drl.addEvidence(configQueue.back(), evidence);
            configQueue.pop_back();

            if (configQueue.empty()) {
              break;
            }

            const std::vector<autopas::Configuration> configQueueCopy = configQueue;
            EXPECT_FALSE(drl.optimizeSuggestions(configQueue, evidenceCollection))
                << "DeepReinforcementLearning should never return 0 during the exploration search.";

            EXPECT_EQ(configQueue.size(), configQueueCopy.size())
                << "DeepReinforcementLearning should not change the size of the config queue.";

            EXPECT_EQ(configQueue, configQueueCopy)
                << "DeepReinforcementLearning should not change the config queue during the exploration search.";
          }

          // Search all exploitation samples
          EXPECT_TRUE(drl.optimizeSuggestions(configQueue, evidenceCollection))
              << "DeepReinforcementLearning should terminate after the exploitation samples.";

          EXPECT_LE(configQueue.size(), 1)
              << "DeepReinforcementLearning should not choose more than 1 exploitation sample.";

          while (!configQueue.empty()) {
            drl.addEvidence(configQueue.back(), evidence);
            configQueue.pop_back();

            const std::vector<autopas::Configuration> configQueueCopy = configQueue;
            EXPECT_TRUE(drl.optimizeSuggestions(configQueue, evidenceCollection))
                << "DeepReinforcementLearning should terminate after the exploitation samples.";

            EXPECT_EQ(configQueue.size(), configQueueCopy.size())
                << "DeepReinforcementLearning should not change the size of the config queue.";

            EXPECT_EQ(configQueue, configQueueCopy)
                << "DeepReinforcementLearning should not change the config queue during the exploitation search.";
          }
        }
      }
    }
  }
}

/**
 * Test that a bad configuration does not starve.
 */
TEST_F(DeepReinforcementLearningTest, badConfigurationDoesNotStarve) {
  // Set the exception policy to throw exceptions
  autopas::utils::ExceptionHandler::setBehavior(autopas::utils::ExceptionBehavior::throwException);

  // Create the core variables
  const auto searchSpace = configQueueBase;

  // This test only applies to the polynomial exploration method. A lower number of exploration sample options is
  // searched, for the test to be faster.
  for (const bool train : {true, false}) {
    for (size_t numExplorationSamples = 3; numExplorationSamples <= 5; numExplorationSamples++) {
      for (const autopas::DeepReinforcementLearning::ExplorationMethod explorationMethod :
           {autopas::DeepReinforcementLearning::ExplorationMethod::polynomial}) {
        autopas::DeepReinforcementLearning drl(train, numExplorationSamples, explorationMethod);
        std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
        autopas::EvidenceCollection evidenceCollection{};
        autopas::Configuration targetConfig = *(searchSpace.begin()++);
        std::unordered_map<autopas::Configuration, autopas::Evidence, autopas::ConfigHash> evidenceMap;

        for (const auto &config : configQueue) {
          evidenceMap[config] = autopas::Evidence{0, 0, config == targetConfig ? 224 : 112};
        }

        // Perform the initial full search
        drl.reset(0, 0, configQueue, evidenceCollection);

        drl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);
        configQueue.pop_back();

        while (true) {
          if (configQueue.empty()) {
            break;
          }

          drl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);
          configQueue.pop_back();

          const std::vector<autopas::Configuration> configQueueCopy = configQueue;
          drl.optimizeSuggestions(configQueue, evidenceCollection);
        }

        int occurrences = 0;

        // Perform the future searches
        for (size_t tuningPhase = 1; tuningPhase < 10000; tuningPhase++) {
          configQueue = {searchSpace.begin(), searchSpace.end()};

          // Reset
          std::vector<autopas::Configuration> configQueueResetCopy = configQueue;
          drl.reset(tuningPhase, tuningPhase * 1000, configQueue, evidenceCollection);

          if (configQueue.back() == targetConfig) {
            occurrences++;
          }

          drl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);

          configQueue.pop_back();

          if (configQueue.empty()) {
            break;
          }

          // Search all exploration samples
          while (true) {
            if (configQueue.back() == targetConfig) {
              occurrences++;
            }

            drl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);
            configQueue.pop_back();

            if (configQueue.empty()) {
              break;
            }

            const std::vector<autopas::Configuration> configQueueCopy = configQueue;
            drl.optimizeSuggestions(configQueue, evidenceCollection);
          }

          // Search all exploitation samples
          drl.optimizeSuggestions(configQueue, evidenceCollection);

          while (!configQueue.empty()) {
            if (configQueue.back() == targetConfig) {
              occurrences++;
            }

            drl.addEvidence(configQueue.back(), evidenceMap[configQueue.back()]);
            configQueue.pop_back();

            const std::vector<autopas::Configuration> configQueueCopy = configQueue;
            drl.optimizeSuggestions(configQueue, evidenceCollection);
          }
        }

        EXPECT_GT(occurrences, 0) << "The worst configuration must not starve.";
      }
    }
  }
}
