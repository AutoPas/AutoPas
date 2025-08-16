/**
 * @file DeepReinforcementLearningTest.cpp
 * @author P. Metscher
 * @date 24.06.25
 */

#include "DeepReinforcementLearningTest.h"

#include "autopas/tuning/utils/SearchSpaceGenerators.h"

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

  for (const bool train : {true, false}) {
    for (size_t explorationSamples = 2; explorationSamples <= 8; explorationSamples++) {
      for (const autopas::DeepReinforcementLearning::ExplorationMethod explorationMethod :
           {autopas::DeepReinforcementLearning::ExplorationMethod::polynomial,
            autopas::DeepReinforcementLearning::ExplorationMethod::random,
            autopas::DeepReinforcementLearning::ExplorationMethod::longestAgo}) {
        autopas::DeepReinforcementLearning drl(train, explorationSamples, explorationMethod);
        std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
        autopas::EvidenceCollection evidenceCollection{};

        // Test the initial full search
        std::vector<autopas::Configuration> configQueueResetCopy = configQueue;
        EXPECT_FALSE(drl.reset(0, 0, configQueue, evidenceCollection))
            << "DeepReinforcementLearning should not clear the queue ever intentionally when resetting.";

        EXPECT_EQ(configQueueResetCopy, configQueue)
            << "DeepReinforcementLearning should not change the config queue during the reset.";

        const autopas::Evidence evidence{0, 0, 100};
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

  for (const bool train : {true, false}) {
    for (size_t explorationSamples = 2; explorationSamples <= 8; explorationSamples++) {
      for (const autopas::DeepReinforcementLearning::ExplorationMethod explorationMethod :
           {autopas::DeepReinforcementLearning::ExplorationMethod::polynomial,
            autopas::DeepReinforcementLearning::ExplorationMethod::random,
            autopas::DeepReinforcementLearning::ExplorationMethod::longestAgo}) {
        autopas::DeepReinforcementLearning drl(train, explorationSamples, explorationMethod);
        std::vector<autopas::Configuration> configQueue{searchSpace.rbegin(), searchSpace.rend()};
        autopas::EvidenceCollection evidenceCollection{};

        // Perform the initial full search
        drl.reset(0, 0, configQueue, evidenceCollection);
        const autopas::Evidence evidence{0, 0, 100};

        while (true) {
          const std::vector<autopas::Configuration> configQueueCopy = configQueue;
          drl.optimizeSuggestions(configQueue, evidenceCollection);

          if (configQueue.empty()) {
            break;
          }

          drl.addEvidence(configQueue.back(), evidence);
          configQueue.pop_back();
        }

        // Perform the future searches
        for (size_t tuningPhase = 1; tuningPhase < 100; tuningPhase++) {
          configQueue = {searchSpace.begin(), searchSpace.end()};

          // Reset
          std::vector<autopas::Configuration> configQueueResetCopy = configQueue;
          EXPECT_FALSE(drl.reset(tuningPhase, tuningPhase * 1000, configQueue, evidenceCollection))
              << "DeepReinforcementLearning should not clear the queue ever intentionally when resetting.";

          EXPECT_EQ(configQueueResetCopy, configQueue)
              << "DeepReinforcementLearning should not change the config queue during the reset.";

          // Choose Exploration samples
          const std::vector<autopas::Configuration> initialConfigQueueCopy = configQueue;
          drl.optimizeSuggestions(configQueue, evidenceCollection);

          EXPECT_EQ(configQueue.size(), explorationSamples)
              << "DeepReinforcementLearning should shrink the config queue to the number of exploration samples.";

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
