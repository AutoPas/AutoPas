/**
 * @file AutoTunerTest.cpp
 * @author F. Gratl
 * @date 8/10/18
 */

#include "AutoTunerTest.h"

#include <vector>

#include "autopas/LogicHandler.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/SlowConfigFilter.h"
#include "autopas/tuning/tuningStrategy/SortByName.h"
#include "autopas/tuning/utils/AutoTunerInfo.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * NOTICE: This class uses always the MockFunctor, even when the mock functionalities are not needed,
 * in order to keep the number of template instantiations of AutoTuner to a minimum.
 */

using ::testing::_;

/**
 * Generates no configurations.
 */
TEST_F(AutoTunerTest, testNoConfig) {
  // wrap experiment into lambda to make simpler EXPECT_THROW expression
  auto experiment = []() {
    constexpr unsigned int rebuildFrequency = 20;
    constexpr autopas::AutoTunerInfo autoTunerInfo{};
    autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
    autopas::AutoTuner::SearchSpaceType searchSpace = {};
    autopas::AutoTuner autoTuner(tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, "");
  };

  EXPECT_THROW(experiment(), autopas::utils::ExceptionHandler::AutoPasException)
      << "The Constructor should catch empty search spaces.";
}

/**
 * Iterate with two configs.
 * First has short rebuild and long non-rebuild iterations
 * Second has long rebuild and short non-rebuild iterations
 * Expect to choose the first because the second one is worse on average.
 */
TEST_F(AutoTunerTest, testBuildNotBuildTimeEstimation) {
  constexpr unsigned int rebuildFrequency = 20;
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .selectorStrategy = autopas::options::SelectorStrategyOption::fastestMean,
      .tuningInterval = 1000,
      .maxSamples = 3};
  // Use configurations with N3, otherwise there are more calls to AoSFunctor
  const auto searchSpace = {_confLc_c08_N3, _confDs_seq_N3};

  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""};

  size_t iteration = 0;
  size_t tuningPhase = 0;

  autopas::LiveInfo info{};

  // Iteration 0, Config 0, Rebuilding Neighbor List
  autoTuner.tuneConfiguration(iteration, tuningPhase, /*isStartOfTuningPhase*/ true);
  const auto config0a = autoTuner.getCurrentConfig();
  EXPECT_TRUE(autoTuner.inTuningPhase()) << "AutoTuner should be in tuning state after just one iteration.";
  autoTuner.addMeasurement(200000, 25000, true, iteration, tuningPhase);
  iteration++;

  // Iteration 1, Config 0, Not Rebuilding
  autoTuner.tuneConfiguration(iteration, tuningPhase, /*isStartOfTuningPhase*/ false);
  const auto config0b = autoTuner.getCurrentConfig();
  // Sanity check that configuration didn't change
  ASSERT_EQ(config0a, config0b);
  EXPECT_TRUE(autoTuner.inTuningPhase()) << "AutoTuner should be in tuning state after second iteration.";
  autoTuner.addMeasurement(0, 30000, false, iteration, tuningPhase);
  iteration++;

  // Iteration 2, Config 0, Not Rebuilding
  autoTuner.tuneConfiguration(iteration, tuningPhase, /*isStartOfTuningPhase*/ false);
  const auto config0c = autoTuner.getCurrentConfig();
  // Sanity check that configuration didn't change
  ASSERT_EQ(config0a, config0c);
  EXPECT_TRUE(autoTuner.inTuningPhase()) << "AutoTuner should be in tuning state after third iteration.";
  autoTuner.addMeasurement(0, 35000, false, iteration, tuningPhase);
  iteration++;

  // Iteration 3, Config 1, Rebuilding Neighbor List
  autoTuner.tuneConfiguration(iteration, tuningPhase, /*isStartOfTuningPhase*/ false);
  const auto config1a = autoTuner.getCurrentConfig();
  // Sanity check that configuration did change
  ASSERT_NE(config0a, config1a);
  EXPECT_TRUE(autoTuner.inTuningPhase()) << "AutoTuner should be in tuning state after fourth iteration.";
  autoTuner.addMeasurement(5500000, 25000, true, iteration, tuningPhase);
  iteration++;

  // Iteration 4, Config 1, Not Rebuilding
  autoTuner.tuneConfiguration(iteration, tuningPhase, /*isStartOfTuningPhase*/ false);
  const auto config1b = autoTuner.getCurrentConfig();
  // Sanity check that configuration didn't change
  ASSERT_EQ(config1a, config1b);
  EXPECT_TRUE(autoTuner.inTuningPhase()) << "AutoTuner should be in tuning state after fifth iteration.";
  autoTuner.addMeasurement(0, 25000, false, iteration, tuningPhase);
  iteration++;

  // Iteration 5, Config 1, Not Rebuilding
  autoTuner.tuneConfiguration(iteration, tuningPhase, /*isStartOfTuningPhase*/ false);
  const auto config1c = autoTuner.getCurrentConfig();
  // Sanity check that configuration didn't change
  ASSERT_EQ(config1a, config1c);
  // This was the last tuning iteration. The AutoTuner stays in tuning mode until the next call to tuneConfiguration()
  EXPECT_TRUE(autoTuner.inTuningPhase()) << "AutoTuner should be in tuning state after sixth iteration.";
  autoTuner.addMeasurement(0, 15000, false, iteration, tuningPhase);
  iteration++;

  autoTuner.tuneConfiguration(iteration, tuningPhase, /*isStartOfTuningPhase*/ false);
  EXPECT_FALSE(autoTuner.inTuningPhase()) << "AutoTuner should no longer be in tuning state.";

  // Expected Weighted Averages
  // Config 1: (200000 + 0 + 0) / 3 / 20 + (25000 + 30000 + 35000) / 3 = 40000
  // Config 2: (5500000 + 0 + 0) / 3 / 20 + (25000 + 25000 + 15000) / 3  = 113333.33
  // => Config 1 is optimal
  EXPECT_EQ(autoTuner.getCurrentConfig(), config0a);
  EXPECT_NE(autoTuner.getCurrentConfig(), config1a);
}

/**
 *  Add fewer measurements than the rebuild frequency and check if the weighted average for the evidence is correct.
 */
TEST_F(AutoTunerTest, testSampleWeightingOneRebuild) {
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  const autopas::AutoTuner::SearchSpaceType searchSpace{_confLc_c08_noN3, _confLc_c01_noN3};
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .maxSamples = 3,
  };
  constexpr size_t rebuildFrequency = 10;
  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""};

  autoTuner.tuneConfiguration(0, 0, /*isStartOfTuningPhase*/ true);
  const auto config = autoTuner.getCurrentConfig();

  constexpr long sampleRebuild = 10;
  constexpr long sampleNonRebuild = 2;
  autoTuner.addMeasurement(sampleRebuild, sampleNonRebuild, true, 0, 0);
  autoTuner.addMeasurement(0, sampleNonRebuild, false, 1, 0);
  autoTuner.addMeasurement(0, sampleNonRebuild, false, 2, 0);

  constexpr long expectedEvidence = sampleRebuild / rebuildFrequency + sampleNonRebuild;
  EXPECT_EQ(expectedEvidence, autoTuner.getEvidenceCollection().getEvidence(config)->front().reducedValue);
}

/**
 *  Add more measurements than the rebuild frequency and check if the weighted average for the evidence is correct.
 *  Version with two rebuilds during sampling.
 */
TEST_F(AutoTunerTest, testSampleWeightingTwoRebuild) {
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  const autopas::AutoTuner::SearchSpaceType searchSpace{_confLc_c08_noN3, _confLc_c01_noN3};
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .maxSamples = 5,
  };
  constexpr size_t rebuildFrequency = 3;
  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""};

  autoTuner.tuneConfiguration(0, 0, /*isStartOfTuningPhase*/ true);
  const auto config = autoTuner.getCurrentConfig();

  constexpr long sampleRebuild = 8;
  constexpr long sampleNonRebuild = 2;
  autoTuner.addMeasurement(sampleRebuild, sampleNonRebuild, true, 0, 0);
  autoTuner.addMeasurement(0, sampleNonRebuild, false, 1, 0);
  autoTuner.addMeasurement(0, sampleNonRebuild, false, 2, 0);
  autoTuner.addMeasurement(sampleRebuild, sampleNonRebuild, true, 3, 0);
  autoTuner.addMeasurement(0, sampleNonRebuild, false, 4, 0);

  constexpr long expectedEvidence = sampleRebuild / rebuildFrequency + sampleNonRebuild;
  EXPECT_EQ(expectedEvidence, autoTuner.getEvidenceCollection().getEvidence(config)->front().reducedValue);
}

/**
 * Test that if a tuning strategy wipes the whole config queue it is not applied.
 */
TEST_F(AutoTunerTest, testRestoreAfterWipe) {
  // Create a tuning strategy that will always reject everything
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  // This strategy will throw out all configurations that are slower than 50% of the fastest
  // (incl. fastest, yes this is bug abuse)
  tuningStrategies.emplace_back(std::make_unique<autopas::SlowConfigFilter>(0.5));
  tuningStrategies.emplace_back(std::make_unique<autopas::SortByName>());

  // Set up the tuner
  const autopas::AutoTuner::SearchSpaceType searchSpace{_confLc_c08_noN3, _confLc_c01_noN3};
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .maxSamples = 1,
  };
  constexpr size_t rebuildFrequency = 3;
  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""};

  size_t iteration = 0;
  // Fill the search space with random data so the slow config filter can work
  for (const auto conf : searchSpace) {
    autoTuner.tuneConfiguration(iteration, 0, iteration == 0);
    autoTuner.addMeasurement(10, 42, true, iteration, 0);
    iteration++;
  }

  // Trigger the tuning process with evidence. Here the slow config filter would wipe out everything
  autoTuner.tuneConfiguration(iteration, 0, false);

  // The slow config filter should have been ignored
  // But the second strategy should still have been applied reversing the order of configs
  EXPECT_EQ(autoTuner.getConfigQueue()[0], _confLc_c01_noN3);
  EXPECT_EQ(autoTuner.getConfigQueue()[1], _confLc_c08_noN3);
}

void AutoTunerTest::testEndingTuningPhaseWithRejectedConfig(bool rejectIndefinitely) const {
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  constexpr autopas::AutoTunerInfo autoTunerInfo{.tuningInterval = 100, .maxSamples = 1};
  constexpr unsigned int rebuildFrequency = 1;

  const autopas::AutoTuner::SearchSpaceType searchSpace{
      _confDs_seq_noN3,
      _confDs_seq_N3,
      _confLc_c08_N3,
      _confLc_c18_N3,
  };

  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, "2B"};

  size_t iteration = 0;
  size_t tuningPhase = 0;

  // 1) Sample first config (accepted).
  autoTuner.tuneConfiguration(iteration, 0, true);
  const auto config1 = autoTuner.getCurrentConfig();
  autoTuner.addMeasurement(180, 20, /*neighborListRebuilt*/ true, iteration, tuningPhase);
  iteration++;

  // 2) Second config is rejected, should immediately get the next one.
  autoTuner.tuneConfiguration(iteration, 0, false);
  const auto config2 = autoTuner.getCurrentConfig();
  const auto configAfterReject2 =
      autoTuner.rejectConfig(config2, /*indefinitely*/ rejectIndefinitely, iteration, tuningPhase);
  // EXPECT_TRUE(stillTuningAfterReject2);

  // 3) Sample third config (accepted). Make it faster than config1 so it becomes the optimum.
  autoTuner.addMeasurement(90, 10, /*neighborListRebuilt*/ true, iteration, tuningPhase);
  iteration++;

  // 4) Final config is rejected.
  autoTuner.tuneConfiguration(iteration, 0, false);
  const auto config4 = autoTuner.getCurrentConfig();
  const auto configAfterReject4 =
      autoTuner.rejectConfig(config4, /*indefinitely*/ rejectIndefinitely, iteration, tuningPhase);

  const auto expectedBest = _confLc_c08_N3;
  EXPECT_EQ(configAfterReject4, expectedBest);

  iteration++;
  autoTuner.tuneConfiguration(iteration, 0, false);
  // We should be able to later call getCurrentConfig and receive the expectedBest
  EXPECT_EQ(autoTuner.getCurrentConfig(), expectedBest);
}

/**
 * Tests that ending a tuning phase with a reject configuration is handled correctly.
 *
 * Mimics an AutoTuner trialling four configurations. The second and the final configurations will be rejected. After
 * the tuning phase is completed, the test should get the correct configuration upon calling `getNextConfig`
 */
TEST_F(AutoTunerTest, testEndingTuningPhaseWithRejectedConfig) {
  testEndingTuningPhaseWithRejectedConfig(true);
  testEndingTuningPhaseWithRejectedConfig(false);
}

/**
 * Test that forcing an optimal configuration immediately ends the tuning phase,
 * clears the config queue, and prevents further sampling.
 */
TEST_F(AutoTunerTest, testForceOptimalConfiguration) {
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  constexpr autopas::AutoTunerInfo autoTunerInfo{.maxSamples = 3};

  // Provide a search space with multiple options
  const auto searchSpace = {_confLc_c08_noN3, _confDs_seq_noN3, _confLc_c18_noN3};
  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, 20, ""};

  autoTuner.tuneConfiguration(0, 0, true);
  EXPECT_TRUE(autoTuner.inTuningPhase());

  // The TuningManager steps in and forces a specific config
  autoTuner.forceOptimalConfiguration(_confDs_seq_noN3);

  // Assert the state machine correctly aborted tuning and locked in the choice
  EXPECT_FALSE(autoTuner.inTuningPhase()) << "Tuner should not be tuning after an optimal config is forced.";
  EXPECT_EQ(autoTuner.getCurrentConfig(), _confDs_seq_noN3) << "Current config does not match the forced config.";

  // Assert that further calls to tune do not mistakenly restart it
  const bool stillTuning = autoTuner.tuneConfiguration(1, 0, false);
  EXPECT_FALSE(stillTuning) << "Tuner should refuse to continue tuning.";
}

/**
 * Test the early stopping mechanic. If a configuration is significantly slower
 * than the currently established optimum, its remaining samples should be skipped.
 */
TEST_F(AutoTunerTest, testEarlyStopping) {
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .maxSamples = 3,
      .earlyStoppingFactor = 2.0  // Stop if a config is > 2x slower than the current best
  };
  const auto searchSpace = {_confLc_c08_noN3, _confDs_seq_noN3};
  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, 20, ""};

  size_t iteration = 0;
  size_t tuningPhase = 0;

  // --- Evaluate First Config (Establish a fast baseline) ---
  autoTuner.tuneConfiguration(iteration, tuningPhase, true);
  const auto fastConfig = autoTuner.getCurrentConfig();

  // Add 3 fast measurements (e.g., 100 ticks)
  autoTuner.addMeasurement(10, 100, true, iteration, tuningPhase);
  iteration++;
  autoTuner.tuneConfiguration(iteration, tuningPhase, false);
  autoTuner.addMeasurement(0, 100, false, iteration, tuningPhase);
  iteration++;
  autoTuner.tuneConfiguration(iteration, tuningPhase, false);
  autoTuner.addMeasurement(0, 100, false, iteration, tuningPhase);
  iteration++;

  // --- Evaluate Second Config (The slow one) ---
  autoTuner.tuneConfiguration(iteration, tuningPhase, false);  // Transition to next config
  const auto slowConfig = autoTuner.getCurrentConfig();
  EXPECT_NE(fastConfig, slowConfig);
  EXPECT_TRUE(autoTuner.inTuningPhase());

  // Add ONE terrible measurement (500 ticks, which is 5x slower than 100)
  autoTuner.addMeasurement(50, 500, true, iteration, tuningPhase);
  iteration++;

  // The next call to tuneConfiguration will trigger checkEarlyStoppingCondition()
  const bool stillTuning = autoTuner.tuneConfiguration(iteration, tuningPhase, false);

  EXPECT_FALSE(stillTuning) << "Tuning should be aborted due to the second config being to slow.";

  // Assert that the tuner instantly gave up on the slow config and ended the phase
  // without asking for the 2nd and 3rd samples!
  EXPECT_FALSE(autoTuner.inTuningPhase()) << "Early stopping failed to abort the tuning phase.";

  // Assert it correctly fell back to the fastest known configuration
  EXPECT_EQ(autoTuner.getCurrentConfig(), fastConfig);
}

/**
 * Test that an AutoTuner initialized with a trivial search space (size == 1)
 * never enters a tuning phase and simply acts as a static configuration provider.
 */
TEST_F(AutoTunerTest, testTrivialSearchSpace) {
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  constexpr autopas::AutoTunerInfo autoTunerInfo{.maxSamples = 3};

  // Provide a search space with EXACTLY ONE option
  const auto searchSpace = {_confLc_c08_noN3};
  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, 20, ""};

  // Assert that it immediately recognizes it has nothing to tune
  EXPECT_FALSE(autoTuner.inTuningPhase()) << "Tuner with size 1 should never report being in a tuning phase.";

  // Assert that explicitly commanding it to start a phase is safely ignored
  const bool stillTuning = autoTuner.tuneConfiguration(0, 0, true);
  EXPECT_FALSE(stillTuning) << "Tuner should refuse to start tuning when search space is trivial.";

  // Assert it returns the only valid configuration
  EXPECT_EQ(autoTuner.getCurrentConfig(), _confLc_c08_noN3);

  EXPECT_FALSE(autoTuner.inTuningPhase())
      << "Tuner should not have entered a tuning phase when search space is trivial.";
}