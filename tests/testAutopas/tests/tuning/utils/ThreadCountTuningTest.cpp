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
#include "autopasTools/generators/GridGenerator.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;

void ThreadCountTuningTest::testThreadCountTuningWithBoxMax(const size_t boxMax,
                                                            const std::set<int> &threadCountOptions,
                                                            const size_t expectedSelectedThreadCount) const {
  const std::set<autopas::ContainerOption> containerOptions({autopas::ContainerOption::linkedCells});
  const std::set<autopas::TraversalOption> traversalOptions({autopas::TraversalOption::lc_c01});
  const std::set<autopas::LoadEstimatorOption> loadEstimatorOptions({autopas::LoadEstimatorOption::none});
  const std::set<autopas::DataLayoutOption> dataLayoutOptions({autopas::DataLayoutOption::soa});
  const std::set<autopas::Newton3Option> newton3Options({autopas::Newton3Option::disabled});
  const autopas::NumberSetFinite<double> cellSizeFactors({1});
  const autopas::NumberSetFinite<int> threadCounts(threadCountOptions);
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
      &threadCounts, autopas::InteractionTypeOption::pairwise);
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::off);
  //  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  bool stillTuning = true;
  autopas::Configuration currentConfig;
  autopasTools::generators::GridGenerator::fillWithParticles(logicHandler.getContainer(), {boxMax, boxMax, boxMax},
                                                             Molecule());
  const size_t numInsertedMolecules = logicHandler.getContainer().size();

  int iterationsAfterTuning = 0;
  while (stillTuning and iterationsAfterTuning < 1) {
    // Should not have any leaving molecules in this test
    auto dummyMoleculesVec = logicHandler.updateContainer();
    stillTuning = logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
    currentConfig = tunerMap[autopas::InteractionTypeOption::pairwise]->getCurrentConfig();
    if (not stillTuning) iterationsAfterTuning++;
  }

  EXPECT_EQ(numInsertedMolecules,
            logicHandler.getContainer().size());  // Should not have any leaving molecules in this test
  // NOTE: currentConfig.threadCount does not return the actual number of threads when thread count tuning is turned off
  // Instead, check the actual number of threads to be used as set by the current configuration
  EXPECT_EQ(expectedSelectedThreadCount, autopas::autopas_get_preferred_num_threads());
}

/**
 * Tests thread count tuning:
 * vary small scenario (8 particles) -> lowest number of threads
 * larger scenario (33k particles)   -> highest number of threads
 */
TEST_F(ThreadCountTuningTest, testThreadCountTuning) {
  const int maxThreads = std::min(autopas::autopas_get_max_threads(), 4);
  testThreadCountTuningWithBoxMax(2, {1, maxThreads}, 1);
  testThreadCountTuningWithBoxMax(32, {1, maxThreads}, maxThreads);
}
