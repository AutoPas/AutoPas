/**
 * @file AutoTunerTest.cpp
 * @author F. Gratl
 * @date 8/10/18
 */

#include "AutoTunerTest.h"

#include <cstddef>
#include <vector>

#include "autopas/AutoPasDecl.h"
#include "autopas/LogicHandler.h"
#include "autopas/LogicHandlerInfo.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/SlowConfigFilter.h"
#include "autopas/tuning/tuningStrategy/SortByName.h"
#include "autopas/tuning/utils/AutoTunerInfo.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/checkFunctorType.h"
#include "autopasTools/generators/GridGenerator.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * NOTICE: This class uses always the MockFunctor, even when the mock functionalities are not needed,
 * in order to keep the number of template instantiations of AutoTuner to a minimum.
 */

using ::testing::_;

TEST_F(AutoTunerTest, testAllConfigurations) {
  const autopas::NumberSetFinite<double> cellSizeFactors({1});
  const double verletSkinPerTimestep = 0;
  const unsigned int verletRebuildFrequency = 20;
  const unsigned int verletClusterSize = 64;
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{8., 8., 8.},
      .cutoff = 1.8,
      .verletSkinPerTimestep = .01,
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 2,
  };

  // the NiceMock wrapper suppresses warnings from uninteresting function calls
  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  // Need to resize cells during loading, otherwise we get exceptions in SoAFunctors
  EXPECT_CALL(functor, SoALoader(::testing::Matcher<autopas::ReferenceParticleCell<Molecule> &>(_), _, _, _))
      .Times(testing::AtLeast(1))
      .WillRepeatedly(
          testing::WithArgs<0, 1>(testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.size()); })));
  EXPECT_CALL(functor, SoALoader(::testing::Matcher<FMCell &>(_), _, _, _))
      .Times(testing::AtLeast(1))
      .WillRepeatedly(
          testing::WithArgs<0, 1>(testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.size()); })));
  const auto searchSpace = autopas::SearchSpaceGenerators::cartesianProduct(
      autopas::ContainerOption::getAllOptions(), autopas::TraversalOption::getAllOptions(),
      autopas::LoadEstimatorOption::getAllOptions(), autopas::DataLayoutOption::getAllOptions(),
      autopas::Newton3Option::getAllOptions(), &cellSizeFactors, autopas::InteractionTypeOption::pairwise);
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::off);
  //  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  bool stillTuning = true;
  auto prevConfig = autopas::Configuration();

  std::map<autopas::ContainerOption, size_t> configsPerContainer;

  // number of configs manually counted:
  //
  // Direct Sum:            ds_sequential               (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  configsPerContainer[autopas::ContainerOption::directSum] = 4;
  // LinkedCells:           lc_c08                      (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        lc_sliced                   (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        lc_sliced_balanced          (AoS <=> SoA, newton3 <=> noNewton3, 2 heuristics)   = 8
  //                        lc_sliced_c02               (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        lc_c18                      (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        lc_c01                      (AoS <=> SoA, noNewton3)                             = 2
  //                        lc_c01_combined_SoA         (SoA, noNewton3)                                     = 1
  //                        lc_c04                      (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        lc_c04_combined_SoA         (SoA, newton3 <=> noNewton3)                         = 2
  //                        lc_c04_HCP                  (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  configsPerContainer[autopas::ContainerOption::linkedCells] = 37;
  // same as linked Cells but load estimator stuff is currently missing
  configsPerContainer[autopas::ContainerOption::linkedCellsReferences] =
      configsPerContainer[autopas::ContainerOption::linkedCells] - 4;
  // VerletLists:           vl_list_iteration           (AoS <=> SoA, noNewton3)                             = 2
  configsPerContainer[autopas::ContainerOption::verletLists] = 2;
  // VerletListsCells:      vlc_sliced                  (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        vlc_sliced_balanced         (AoS <=> SoA, newton3 <=> noNewton3, 3 heuristics)   = 12
  //                        vlc_sliced_colored          (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        vlc_c18                     (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        vlc_c01                     (AoS <=> SoA, noNewton3)                             = 2
  configsPerContainer[autopas::ContainerOption::verletListsCells] = 26;
  // VerletClusterLists:    vcl_cluster_iteration       (AoS <=> SoA, noNewton3)                             = 2
  //                        vcl_c06                     (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        vcl_c01_balanced            (AoS <=> SoA, noNewton3)                             = 2
  //                        vcl_sliced                  (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        vcl_sliced_c02              (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        vcl_sliced_balanced         (AoS <=> SoA, newton3 <=> noNewton3, 2 heuristics)   = 8
  configsPerContainer[autopas::ContainerOption::verletClusterLists] = 24;
  // VarVerletListsAsBuild: vvl_as_built                (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  configsPerContainer[autopas::ContainerOption::varVerletListsAsBuild] = 4;

  // PairwiseVerletLists:   vlp_sliced                  (AoS <=> SoA, newton3 <=> noNewton3)                         = 4
  //                        vlp_sliced_balanced         (AoS <=> SoA, newton3 <=> noNewton3)                         = 4
  //                        vlp_sliced_colored          (AoS <=> SoA, newton3 <=> noNewton3)                         = 4
  //                        vlp_c18                     (AoS <=> SoA, newton3 <=> noNewton3)                         = 4
  //                        vlp_c01                     (AoS <=> SoA, noNewton3)                                     = 2
  //                        vlp_c08                     (AoS <=> SoA, newton3 <=> noNewton3)                         = 4
  configsPerContainer[autopas::ContainerOption::pairwiseVerletLists] = 22;

  // Octree:                ot_c01                      (AoS <=> SoA, noNewton3)                             = 2
  //                        ot_c18                      (AoS <=> SoA, newton3)                               = 2
  configsPerContainer[autopas::ContainerOption::octree] = 4;

  // check that there is an entry for every container.
  ASSERT_EQ(configsPerContainer.size(), autopas::ContainerOption::getAllOptions().size());

  const size_t numberOfConfigs = std::accumulate(configsPerContainer.begin(), configsPerContainer.end(), 0ul,
                                                 [](auto acc, auto &pair) { return acc + pair.second; });
  ASSERT_EQ(numberOfConfigs, searchSpace.size())
      << "The calculated number of configurations is not equal to the cross product search space!";
  // total number of possible configurations * number of samples + last iteration after tuning
  const size_t expectedNumberOfIterations = numberOfConfigs * autoTunerInfo.maxSamples + 1;

  int collectedSamples = 0;
  int iterations = 0;
  while (stillTuning) {
    if (collectedSamples == autoTunerInfo.maxSamples) {
      collectedSamples = 0;
      logicHandler.getContainer().deleteAllParticles();
      // add particles, so VerletClusterLists uses more than one tower, otherwise its traversals are invalid.
      if (logicHandler.getContainer().getContainerType() == autopas::ContainerOption::verletClusterLists) {
        const std::array<size_t, 3> particlesPerDim = {8, 16, 8};
        const std::array<double, 3> spacing = {0.25, 0.25, 0.25};
        const std::array<double, 3> offset = {0.125, 0.125, 0.125};
        Molecule defaultParticle{};
        autopasTools::generators::GridGenerator::fillWithParticles(logicHandler.getContainer(), particlesPerDim,
                                                                   defaultParticle, spacing, offset);
      }
    }
    stillTuning =
        logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);
    ++iterations;
    ++collectedSamples;
    const auto currentConfig = tunerMap[autopas::InteractionTypeOption::pairwise]->getCurrentConfig();
    if (stillTuning) {
      if (collectedSamples == 1) {
        EXPECT_NE(currentConfig, prevConfig)
            << "Iteration: " << iterations << " collectedSamples: " << collectedSamples;
      } else {
        EXPECT_EQ(currentConfig, prevConfig)
            << "Iteration: " << iterations << " collectedSamples: " << collectedSamples;
      }
    }
    prevConfig = currentConfig;
  }

  EXPECT_EQ(expectedNumberOfIterations, iterations);
}

TEST_F(AutoTunerTest, testWillRebuildDDL) {
  // also check if rebuild is detected if next config is invalid

  const double cellSizeFactor = 1.;
  const unsigned int verletRebuildFrequency = 20;
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 2,
  };
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const autopas::AutoTuner::SearchSpaceType searchSpace{
      _confDs_seq_noN3,
      _confDs_seq_N3,
      _confLc_c08_noN3,
  };

  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  auto &autoTuner = *tunerMap[autopas::InteractionTypeOption::pairwise];

  EXPECT_EQ(*(searchSpace.rbegin()), autoTuner.getCurrentConfig());

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "Expect rebuild for first iteration.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // DS NoN3
  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // DS NoN3
  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "Expect rebuild because we change config.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // DS N3
  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // DS N3
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "Expect rebuild because we change config.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // LC NoN3
  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // LC NoN3
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "Expect rebuild because reached end of tuning phase.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // optimum
  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "Expect no rebuild because not tuning.";
}

/**
 * This test simulates that the next config (which is checked by willRebuildNeighborLists) is invalid.
 */
TEST_F(AutoTunerTest, testWillRebuildDDLOneConfigKicked) {
  // also check if rebuild is detected if next config is invalid

  const double cellSizeFactor = 1.;
  const unsigned int verletRebuildFrequency = 20;
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 2,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const autopas::AutoTuner::SearchSpaceType searchSpace{
      _confDs_seq_noN3,
      _confDs_seq_N3,
      _confLc_c08_N3,
  };

  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  auto &autoTuner = *tunerMap[autopas::InteractionTypeOption::pairwise];

  EXPECT_EQ(*(searchSpace.rbegin()), autoTuner.getCurrentConfig());

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(false));

  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "Expect rebuild for first iteration.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // DS N3
  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // DS N3
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "Expect rebuild because we change config.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // LC N3
  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // LC N3
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "Expect rebuild because reached end of tuning phase.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // optimum
  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "Expect no rebuild because not tuning.";
}

TEST_F(AutoTunerTest, testWillRebuildDL) {
  // also check if rebuild is detected if next config is invalid

  const double cellSizeFactor = 1.;
  const unsigned int verletRebuildFrequency = 20;
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 2,
  };
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  const autopas::AutoTuner::SearchSpaceType searchSpace{
      _confDs_seq_noN3,
      _confLc_c08_noN3,
  };

  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  auto &autoTuner = *tunerMap[autopas::InteractionTypeOption::pairwise];

  EXPECT_EQ(*(searchSpace.rbegin()), autoTuner.getCurrentConfig());

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "Expect rebuild for first iteration.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // DS NoN3
  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // DS NoN3
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "Expect rebuild because we change config.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // LC NoN3
  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // LC NoN3
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "Expect rebuild because reached end of tuning phase.";
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                              autopas::InteractionTypeOption::pairwise);  // optimum
  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "Expect no rebuild because not tuning.";
}

/**
 * Initialize a tuner, do a tuning phase, then instead of waiting for the full tuning interval restart the tuning
 * earlier via forceRetune.
 */
TEST_F(AutoTunerTest, testForceRetuneBetweenPhases) {
  const unsigned int verletRebuildFrequency = 20;
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{2, 4, 2},
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };

  autopas::AutoTuner::SearchSpaceType searchSpace{_confLc_c01_noN3, _confLc_c18_noN3, _confLc_c08_noN3};
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  auto &autoTuner = *tunerMap[autopas::InteractionTypeOption::pairwise];

  const size_t numExpectedTuningIterations = searchSpace.size() * autoTunerInfo.maxSamples;
  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  // expect a full tuning phase
  for (size_t i = 0; i < numExpectedTuningIterations; ++i) {
    // since we don't actually do anything doRebuild can always be false.
    EXPECT_TRUE((logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                                             autopas::InteractionTypeOption::pairwise)))
        << "Tuner should still be tuning in iteration " << i;
  }
  // first iteration after tuning phase
  EXPECT_FALSE(
      (logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise)))
      << "Tuner should be done be tuning in the first iteration after the tuning phase.";

  EXPECT_FALSE(autoTuner.willRebuildNeighborLists()) << "No rebuilding expected here.";
  // instead of waiting the full tuning interval restart tuning immediately
  autoTuner.forceRetune();
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "willRebuildNeighborLists() does not recognize forceRetune()";

  // expect a full tuning phase
  for (size_t i = 0; i < numExpectedTuningIterations; ++i) {
    // since we don't actually do anything doRebuild can always be false.
    EXPECT_TRUE((logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                                             autopas::InteractionTypeOption::pairwise)))
        << "Tuner should still be tuning.";
  }
  // first iteration after tuning phase
  EXPECT_FALSE(
      (logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise)))
      << "Tuner should be done be tuning.";
}

TEST_F(AutoTunerTest, testForceRetuneInPhase) {
  const double cellSizeFactor = 1;
  const unsigned int verletRebuildFrequency = 20;
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{2, 4, 2},
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const auto searchSpace = {_confLc_c01_noN3, _confLc_c18_noN3, _confLc_c08_noN3};

  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  auto &autoTuner = *tunerMap[autopas::InteractionTypeOption::pairwise];

  size_t numExpectedTuningIterations = searchSpace.size() * autoTunerInfo.maxSamples;
  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  // Do part of the tuning phase. After the loop we should be in the middle of sampling the second configuration.
  ASSERT_GT(autoTunerInfo.maxSamples, 1);
  ASSERT_GT(searchSpace.size(), 1);
  size_t iteration = 0;
  for (; iteration < autoTunerInfo.maxSamples + 1; ++iteration) {
    // since we don't actually do anything doRebuild can always be false.
    EXPECT_TRUE((logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                                             autopas::InteractionTypeOption::pairwise)))
        << "Tuner should still be tuning.\n"
           "Phase 1\n"
           "Iteration "
        << iteration;
  }
  // restart the full tuning phase
  autoTuner.forceRetune();
  EXPECT_TRUE(autoTuner.willRebuildNeighborLists()) << "willRebuildNeighborLists() does not recognize forceRetune()";

  // expect a full tuning phase
  for (size_t i = 0; i < numExpectedTuningIterations; ++i, ++iteration) {
    // since we don't actually do anything doRebuild can always be false.
    EXPECT_TRUE((logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor,
                                                                             autopas::InteractionTypeOption::pairwise)))
        << "Tuner should still be tuning.\n"
           "Phase 2\n"
           "Iteration "
        << iteration;
  }
  // first iteration after tuning phase
  EXPECT_FALSE(
      (logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise)))
      << "Tuner should be done be tuning.\n"
         "Iteration "
      << iteration;
}

/**
 * Generates no configurations.
 */
TEST_F(AutoTunerTest, testNoConfig) {
  // wrap experiment into lambda to make simpler EXPECT_THROW expression
  auto experiment = []() {
    const unsigned int verletRebuildFrequency = 20;
    const autopas::AutoTunerInfo autoTunerInfo{};
    autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

    autopas::AutoTuner::SearchSpaceType searchSpace = {};
    autopas::AutoTuner autoTuner(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, "");
  };

  EXPECT_THROW(experiment(), autopas::utils::ExceptionHandler::AutoPasException)
      << "The Constructor should catch empty search spaces.";
}

/**
 * Generates exactly one valid configuration.
 */
TEST_F(AutoTunerTest, testOneConfig) {
  const double cellSizeFactor = 1.;
  const unsigned int verletRebuildFrequency = 20;
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  const auto searchSpace = {_confLc_c08_noN3};
  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  auto &autoTuner = *tunerMap[autopas::InteractionTypeOption::pairwise];

  EXPECT_EQ(_confLc_c08_noN3, autoTuner.getCurrentConfig());

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  size_t numSamples = 0;
  for (int i = 0; i < 5; ++i) {
    if (numSamples == autoTunerInfo.maxSamples) {
      numSamples = 0;
    }
    logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);
    ++numSamples;
    EXPECT_EQ(_confLc_c08_noN3, autoTuner.getCurrentConfig());
  }
}

/**
 * Generates exactly one valid and one invalid configuration.
 */
TEST_F(AutoTunerTest, testConfigSecondInvalid) {
  const double cellSizeFactor = 1.;
  const unsigned int verletRebuildFrequency = 20;
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  const auto searchSpace = {_confLc_c08_noN3, _confLc_c08_N3};
  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  auto &autoTuner = *tunerMap[autopas::InteractionTypeOption::pairwise];

  EXPECT_EQ(*(std::rbegin(searchSpace)), autoTuner.getCurrentConfig());

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;

  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(false));

  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);
  EXPECT_EQ(_confLc_c08_N3, autoTuner.getCurrentConfig());
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);
  EXPECT_EQ(_confLc_c08_N3, autoTuner.getCurrentConfig());
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);
  EXPECT_EQ(_confLc_c08_N3, autoTuner.getCurrentConfig());
}

/**
 * All generated configurations are thrown out at runtime.
 */
TEST_F(AutoTunerTest, testLastConfigThrownOut) {
  const unsigned int verletRebuildFrequency = 20;
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const auto searchSpace = {_confLc_c08_noN3, _confLc_c18_noN3};
  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  auto &autoTuner = *tunerMap[autopas::InteractionTypeOption::pairwise];

  EXPECT_EQ(*std::rbegin(searchSpace), autoTuner.getCurrentConfig());

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(false));

  EXPECT_THROW(
      (logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise)),
      autopas::utils::ExceptionHandler::AutoPasException);
}

/**
 * Iterate with two configs.
 * First has short rebuild and long non-rebuild iterations
 * Second has long rebuild and short non-rebuild iterations
 * Expect to choose the first because the second one is worse on average.
 */
TEST_F(AutoTunerTest, testBuildNotBuildTimeEstimation) {
  const unsigned int verletRebuildFrequency = 20;
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 2,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  // Use configurations with N3, otherwise there are more calls to AoSFunctor
  const auto searchSpace = {_confLc_c08_N3, _confDs_seq_N3};
  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<Molecule> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
  auto &autoTuner = *tunerMap[autopas::InteractionTypeOption::pairwise];

  using ::testing::_;
  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, SoALoader(::testing::Matcher<autopas::FullParticleCell<Molecule> &>(_), _, _, _))
      .Times(testing::AtLeast(0));
  EXPECT_CALL(functor, SoAExtractor(::testing::Matcher<autopas::FullParticleCell<Molecule> &>(_), _, _))
      .Times(testing::AtLeast(0));
  EXPECT_CALL(functor, SoAFunctorPair(_, _, _)).Times(testing::AtLeast(0));
  EXPECT_CALL(functor, SoAFunctorSingle(_, _)).Times(testing::AtLeast(0));
  logicHandler.getContainer().addParticle((Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0}));
  logicHandler.getContainer().addParticle((Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0}));

  using namespace std::literals;

  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(100ms); }));
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);

  auto firstConfig = autoTuner.getCurrentConfig();

  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(30ms); }));
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);

  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(30ms); }));
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);

  // Here, second config will start to be tuned

  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(300ms); }));
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);

  auto secondConfig = autoTuner.getCurrentConfig();

  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(25ms); }));
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);
  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(25ms); }));
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);

  // Here, tuning should be finished and first should have been chosen (100 + 2 * 30 = 160 < 350 = 300 + 2 * 25)
  EXPECT_CALL(functor, AoSFunctor).Times(1);
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);

  EXPECT_EQ(autoTuner.getCurrentConfig(), firstConfig);
  EXPECT_NE(autoTuner.getCurrentConfig(), secondConfig);
}

/**
 *  Add less measurements than the rebuild frequency and check if the weighted average for the evidence is correct.
 */
TEST_F(AutoTunerTest, testSampleWeightingOneRebuild) {
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  autopas::AutoTuner::SearchSpaceType searchSpace{_confLc_c08_noN3, _confLc_c01_noN3};
  const autopas::AutoTunerInfo autoTunerInfo{
      .maxSamples = 3,
  };
  constexpr size_t rebuildFrequency = 10;
  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""};

  const auto [config, _] = autoTuner.getNextConfig();

  constexpr long sampleWithRebuild = 10;
  constexpr long sampleWithoutRebuild = 2;
  autoTuner.addMeasurement(sampleWithRebuild, true);
  autoTuner.addMeasurement(sampleWithoutRebuild, false);
  autoTuner.addMeasurement(sampleWithoutRebuild, false);

  constexpr long expectedEvidence =
      (sampleWithRebuild + (rebuildFrequency - 1) * sampleWithoutRebuild) / rebuildFrequency;
  EXPECT_EQ(expectedEvidence, autoTuner.getEvidenceCollection().getEvidence(config)->front().value);
}

/**
 *  Add more measurements than the rebuild frequency and check if the weighted average for the evidence is correct.
 *  Version with two rebuilds during sampling.
 */
TEST_F(AutoTunerTest, testSampleWeightingTwoRebuild) {
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  autopas::AutoTuner::SearchSpaceType searchSpace{_confLc_c08_noN3, _confLc_c01_noN3};
  const autopas::AutoTunerInfo autoTunerInfo{
      .maxSamples = 5,
  };
  constexpr size_t rebuildFrequency = 3;
  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""};

  const auto [config, _] = autoTuner.getNextConfig();

  constexpr long sampleWithRebuild = 10;
  constexpr long sampleWithoutRebuild = 2;
  autoTuner.addMeasurement(sampleWithRebuild, true);
  autoTuner.addMeasurement(sampleWithoutRebuild, false);
  autoTuner.addMeasurement(sampleWithoutRebuild, false);
  autoTuner.addMeasurement(sampleWithRebuild, true);
  autoTuner.addMeasurement(sampleWithoutRebuild, false);

  constexpr long expectedEvidence =
      (sampleWithRebuild + (rebuildFrequency - 1) * sampleWithoutRebuild) / rebuildFrequency;
  EXPECT_EQ(expectedEvidence, autoTuner.getEvidenceCollection().getEvidence(config)->front().value);
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
  const autopas::AutoTunerInfo autoTunerInfo{
      .maxSamples = 1,
  };
  constexpr size_t rebuildFrequency = 3;
  autopas::AutoTuner autoTuner{tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""};

  // Fill the search space with random data so the slow config filter can work
  for (const auto conf : searchSpace) {
    const auto iDontCare = autoTuner.getNextConfig();
    autoTuner.addMeasurement(42, true);
  }

  // Trigger the tuning process with evidence. Here the slow config filter would wipe out everything
  const auto iDontCare = autoTuner.getNextConfig();

  // The slow config filter should have been ignored
  // But the second strategy should still have been applied reversing the order of configs
  EXPECT_EQ(autoTuner.getConfigQueue()[0], _confLc_c01_noN3);
  EXPECT_EQ(autoTuner.getConfigQueue()[1], _confLc_c08_noN3);
}