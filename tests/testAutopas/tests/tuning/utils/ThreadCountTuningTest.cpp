/**
 * @file ThreadCountTuningTest.cpp
 * @author R. Horn
 * @date 06/03/2026
 */

#include "ThreadCountTuningTest.h"

#include "autopas/LogicHandler.h"
#include "autopas/LogicHandlerInfo.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/utils/AutoTunerInfo.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "generators/src/GridGenerator.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;

void ThreadCountTuningTest::testThreadCountTuningWithBoxMax(const size_t boxMax,
                                                            const std::set<int> &threadCountOptions,
                                                            int expectedSelectedThreadCount) const {
  const std::set<autopas::ContainerOption> containerOptions({autopas::ContainerOption::linkedCells});
  const std::set<autopas::TraversalOption> traversalOptions({autopas::TraversalOption::lc_c01});
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions({autopas::LoadEstimatorOption::none});
  const std::set<autopas::DataLayoutOption> dataLayoutOptions({autopas::DataLayoutOption::soa});
  const std::set<autopas::Newton3Option> newton3Options({autopas::Newton3Option::disabled});
  const autopas::NumberSetFinite<double> cellSizeFactors({1});
  const autopas::NumberSetFinite<int> threadCounts(threadCountOptions);
  const std::set<autopas::VectorizationPatternOption> vecPatternOptions({autopas::VectorizationPatternOption::p1xVec});
  const unsigned int verletRebuildFrequency = 20;
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{static_cast<double>(boxMax), static_cast<double>(boxMax), static_cast<double>(boxMax)},
      .cutoff = 1.8,
      .verletSkin = .2,
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 10,  // Ensure results are not too flaky
  };

  mdLib::LJFunctor<Molecule> functor(logicHandlerInfo.cutoff);

  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      containerOptions, traversalOptions, loadEstimatorOptions, dataLayoutOptions, newton3Options, &cellSizeFactors,
      &threadCounts, vecPatternOptions, autopas::InteractionTypeOption::pairwise);
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  auto tunerManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tunerManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);
  autopas::LogicHandler<Molecule> logicHandler(tunerManager, logicHandlerInfo, verletRebuildFrequency, "");

  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::off);
  //  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  bool stillTuning = true;
  autopasTools::generators::GridGenerator::fillWithParticles(logicHandler.getContainer(), {boxMax, boxMax, boxMax},
                                                             Molecule());
  const size_t numInsertedMolecules = logicHandler.getContainer().size();

  int iterationsAfterTuning = 0;
  while (stillTuning and iterationsAfterTuning < 1) {
    // Should not have any leaving molecules in this test
    auto dummyMoleculesVec = logicHandler.updateContainer();
    stillTuning = logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
    if (not stillTuning) iterationsAfterTuning++;
  }

  EXPECT_EQ(numInsertedMolecules,
            logicHandler.getContainer().size());  // Should not have any leaving molecules in this test
  int selectedThreadCount = autopas::autopas_get_tuned_num_threads();
  if (expectedSelectedThreadCount == 0) {
    // Tuning should be disabled (Only check if all threads are used)
    expectedSelectedThreadCount = autopas::autopas_get_max_threads();
  } else {
    // Check if the actual number of threads is in the set of options
    EXPECT_TRUE(threadCountOptions.contains(selectedThreadCount));
  }
  // Check the actual number of threads to be used as set by the current configuration
  if (expectedSelectedThreadCount > 0) {
    EXPECT_EQ(expectedSelectedThreadCount, selectedThreadCount);
  }
}

/**
 * Tests: * very small scenario (8 particles) -> lowest number of threads
 */
TEST_F(ThreadCountTuningTest, testThreadCountTuningSmall) {
  const int maxThreads = std::min(autopas::autopas_get_max_threads(), 4);
  testThreadCountTuningWithBoxMax(2, {1, maxThreads}, 1);
}

/**
 * Tests: Only select valid options
 * excluding possible default options (1, max threads)
 */
TEST_F(ThreadCountTuningTest, testThreadCountTuningOptions) {
  const int maxThreads = autopas::autopas_get_max_threads();
  if (maxThreads < 3) {
    GTEST_SKIP() << "Too few physical threads";
  }
  std::set<int> threadCountOptions{2, maxThreads - 1};
  testThreadCountTuningWithBoxMax(2, threadCountOptions);
}

/**
 * Tests: Tuning disabled -> use all threads
 */
TEST_F(ThreadCountTuningTest, testThreadCountTuningDisabled) { testThreadCountTuningWithBoxMax(2, {0}, 0); }

/**
 * Tests: larger scenario (33k particles)   -> highest number of threads
 * Sensitive to shared CPU usage! (Skipped in CI.)
 * To run the test explicitly, use:
 *   ./build/tests/testAutopas/runTests \
 *     --gtest_filter=ThreadCountTuningTest.DISABLED_testThreadCountTuningLarge \
 *     --gtest_also_run_disabled_tests
 */
TEST_F(ThreadCountTuningTest, DISABLED_testThreadCountTuningLarge) {
  const int maxThreads = std::min(autopas::autopas_get_max_threads(), 4);
  testThreadCountTuningWithBoxMax(32, {1, maxThreads}, maxThreads);
}
