/**
 * @file AutoTunerTest.cpp
 * @author F. Gratl
 * @date 8/10/18
 */

#include "AutoTunerTest.h"

#include "autopas/selectors/AutoTuner.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/selectors/tuningStrategy/FullSearch.h"
#include "autopasTools/generators/GridGenerator.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * NOTICE: This class uses always the MockFunctor, even when the mock functionalities are not needed,
 * in order to keep the number of template instantiations of AutoTuner to a minimum.
 */

using ::testing::_;

TEST_F(AutoTunerTest, testAllConfigurations) {
  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {2, 4, 2};
  // adaptive domain size so sliced is always applicable.
  // bBoxMax[1] = autopas::autopas::autopas_get_max_threads() * 2;
  const double cutoff = 1;
  const double cellSizeFactor = 1;
  const double verletSkinPerTimestep = 0;
  const unsigned verletRebuildFrequency = 20;
  const unsigned int verletClusterSize = 64;
  const double mpiTuningMaxDifferenceForBucket = 0.3;
  const double mpiTuningWeightForMaxDensity = 0.0;
  const unsigned int maxSamples = 2;
  // the NiceMock wrapper suppresses warnings from uninteresting function calls
  testing::NiceMock<MockFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  // Need to resize cells during loading, otherwise we get exceptions in SoAFunctors
  EXPECT_CALL(functor, SoALoader(::testing::Matcher<autopas::ReferenceParticleCell<Molecule> &>(_), _, _))
      .Times(testing::AtLeast(1))
      .WillRepeatedly(testing::WithArgs<0, 1>(
          testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.numParticles()); })));
  EXPECT_CALL(functor, SoALoader(::testing::Matcher<FMCell &>(_), _, _))
      .Times(testing::AtLeast(1))
      .WillRepeatedly(testing::WithArgs<0, 1>(
          testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.numParticles()); })));
  auto tuningStrategy = std::make_unique<autopas::FullSearch>(
      autopas::ContainerOption::getAllOptions(), std::set<double>({cellSizeFactor}),
      autopas::TraversalOption::getAllOptions(), autopas::LoadEstimatorOption::getAllOptions(),
      autopas::DataLayoutOption::getAllOptions(), autopas::Newton3Option::getAllOptions());
  autopas::AutoTuner<Molecule> autoTuner(bBoxMin, bBoxMax, cutoff, verletSkinPerTimestep, verletClusterSize,
                                         std::move(tuningStrategy), mpiTuningMaxDifferenceForBucket,
                                         mpiTuningWeightForMaxDensity, autopas::SelectorStrategyOption::fastestAbs, 100,
                                         maxSamples, verletRebuildFrequency);

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

  size_t numberOfConfigs = std::accumulate(configsPerContainer.begin(), configsPerContainer.end(), 0ul,
                                           [](auto acc, auto &pair) { return acc + pair.second; });

  // total number of possible configurations * number of samples + last iteration after tuning
  const size_t expectedNumberOfIterations = numberOfConfigs * maxSamples + 1;

  int collectedSamples = 0;
  int iterations = 0;
  bool doRebuild = true;  // defines whether the verlet lists should be rebuild.
  while (stillTuning) {
    if (collectedSamples == maxSamples) {
      collectedSamples = 0;
      doRebuild = true;
      // add particles, so VerletClusterLists uses more than one tower.
      autoTuner.getContainer().get()->deleteAllParticles();
      const std::array<size_t, 3> particlesPerDim = {8, 16, 8};
      const std::array<double, 3> spacing = {0.25, 0.25, 0.25};
      const std::array<double, 3> offset = {0.125, 0.125, 0.125};
      Molecule defaultParticle{};
      autopasTools::generators::GridGenerator::fillWithParticles(*(autoTuner.getContainer().get()), particlesPerDim,
                                                                 defaultParticle, spacing, offset);
    }
    std::vector<std::vector<Molecule>> emptyVec(autopas::autopas_get_max_threads());
    stillTuning = autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);
    doRebuild = false;
    ++iterations;
    ++collectedSamples;
    auto currentConfig = autoTuner.getCurrentConfig();
    if (stillTuning) {
      if (collectedSamples == 1) {
        EXPECT_NE(currentConfig, prevConfig)
            << "current:" << currentConfig.toString() << ", previous: " << prevConfig.toString() << std::endl;
      } else {
        EXPECT_EQ(currentConfig, prevConfig)
            << "current:" << currentConfig.toString() << ", previous: " << prevConfig.toString() << std::endl;
      }
    }
    prevConfig = currentConfig;
  }

  EXPECT_EQ(expectedNumberOfIterations, iterations);
}

TEST_F(AutoTunerTest, testWillRebuildDDL) {
  // also check if rebuild is detected if next config is invalid

  double cellSizeFactor = 1.;
  std::set<autopas::Configuration> configs;
  configs.emplace(autopas::ContainerOption::directSum, cellSizeFactor, autopas::TraversalOption::ds_sequential,
                  autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled);
  configs.emplace(autopas::ContainerOption::directSum, cellSizeFactor, autopas::TraversalOption::ds_sequential,
                  autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled);
  configs.emplace(autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
                  autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled);

  auto tuningStrategy = std::make_unique<autopas::FullSearch>(configs);
  autopas::AutoTuner<Molecule> autoTuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy), 0.3, 0.0,
                                         autopas::SelectorStrategyOption::fastestAbs, 1000, 2, 20);

  EXPECT_EQ(*(configs.begin()), autoTuner.getCurrentConfig());

  MockFunctor<Molecule> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  std::vector<std::vector<Molecule>> emptyVec(autopas::autopas_get_max_threads());

  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild for first iteration.";
  bool doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // DS NoN3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // DS NoN3
  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because we change config.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // DS N3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // DS N3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because we change config.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // LC NoN3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // LC NoN3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because reached end of tuning phase.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // optimum
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because not tuning.";
}

/**
 * This test simulates that the next config (which is checked by willRebuild) is invalid.
 */
TEST_F(AutoTunerTest, testWillRebuildDDLOneConfigKicked) {
  // also check if rebuild is detected if next config is invalid

  double cellSizeFactor = 1.;
  std::set<autopas::Configuration> configs;
  configs.emplace(autopas::ContainerOption::directSum, cellSizeFactor, autopas::TraversalOption::ds_sequential,
                  autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled);
  configs.emplace(autopas::ContainerOption::directSum, cellSizeFactor, autopas::TraversalOption::ds_sequential,
                  autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled);
  configs.emplace(autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
                  autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled);

  auto tuningStrategy = std::make_unique<autopas::FullSearch>(configs);
  autopas::AutoTuner<Molecule> autoTuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy), 0.3, 0.0,
                                         autopas::SelectorStrategyOption::fastestAbs, 1000, 2, 20);

  EXPECT_EQ(*(configs.begin()), autoTuner.getCurrentConfig());

  MockFunctor<Molecule> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(false));

  std::vector<std::vector<Molecule>> emptyVec(autopas::autopas_get_max_threads());

  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild for first iteration.";
  bool doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // DS N3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // DS N3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because we change config.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // LC N3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // LC N3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because reached end of tuning phase.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // optimum
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because not tuning.";
}

TEST_F(AutoTunerTest, testWillRebuildDL) {
  // also check if rebuild is detected if next config is invalid

  double cellSizeFactor = 1.;
  std::set<autopas::Configuration> configs;
  configs.emplace(autopas::ContainerOption::directSum, cellSizeFactor, autopas::TraversalOption::ds_sequential,
                  autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled);
  configs.emplace(autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
                  autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled);

  auto tuningStrategy = std::make_unique<autopas::FullSearch>(configs);
  autopas::AutoTuner<Molecule> autoTuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy), 0.3, 0.0,
                                         autopas::SelectorStrategyOption::fastestAbs, 1000, 2, 20);

  EXPECT_EQ(*(configs.begin()), autoTuner.getCurrentConfig());

  MockFunctor<Molecule> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  std::vector<std::vector<Molecule>> emptyVec(autopas::autopas_get_max_threads());

  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild for first iteration.";
  bool doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // DS NoN3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // DS NoN3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because we change config.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // LC NoN3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // LC NoN3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because reached end of tuning phase.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);  // optimum
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because not tuning.";
}

/**
 * Initialize a tuner, do a tuning phase, then instead of waiting for the full tuning interval restart the tuning
 * earlier via forceRetune.
 */
TEST_F(AutoTunerTest, testForceRetuneBetweenPhases) {
  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {2, 4, 2};
  const double cutoff = 1;
  const double verletSkinPerTimestep = 0;
  const unsigned int verletRebuildFrequency = 20;
  const unsigned int verletClusterSize = 4;
  const double mpiTuningMaxDifferenceForBucket = 0.3;
  const double mpiTuningWeightForMaxDensity = 0.0;
  const unsigned int maxSamples = 3;

  auto configsList = {_confLc_c01, _confLc_c04, _confLc_c08};

  auto tuningStrategy = std::make_unique<autopas::FullSearch>(configsList);
  autopas::AutoTuner<Molecule> autoTuner(bBoxMin, bBoxMax, cutoff, verletSkinPerTimestep, verletClusterSize,
                                         std::move(tuningStrategy), mpiTuningMaxDifferenceForBucket,
                                         mpiTuningWeightForMaxDensity, autopas::SelectorStrategyOption::fastestAbs, 100,
                                         maxSamples, verletRebuildFrequency);

  size_t numExpectedTuningIterations = configsList.size() * maxSamples;
  MockFunctor<Molecule> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  std::vector<std::vector<Molecule>> emptyVec(autopas::autopas_get_max_threads());

  // expect a full tuning phase
  for (size_t i = 0; i < numExpectedTuningIterations; ++i) {
    // since we don't actually do anything doRebuild can always be false.
    EXPECT_TRUE(autoTuner.iteratePairwise(&functor, false, emptyVec, emptyVec)) << "Tuner should still be tuning.";
  }
  // first iteration after tuning phase
  EXPECT_FALSE(autoTuner.iteratePairwise(&functor, false, emptyVec, emptyVec)) << "Tuner should be done be tuning.";

  EXPECT_FALSE(autoTuner.willRebuild()) << "No rebuilding expected here.";
  // instead of waiting the full tuning interval restart tuning immediately
  autoTuner.forceRetune();
  EXPECT_TRUE(autoTuner.willRebuild()) << "willRebuild() does not recognize forceRetune()";

  // expect a full tuning phase
  for (size_t i = 0; i < numExpectedTuningIterations; ++i) {
    // since we don't actually do anything doRebuild can always be false.
    EXPECT_TRUE(autoTuner.iteratePairwise(&functor, false, emptyVec, emptyVec)) << "Tuner should still be tuning.";
  }
  // first iteration after tuning phase
  EXPECT_FALSE(autoTuner.iteratePairwise(&functor, false, emptyVec, emptyVec)) << "Tuner should be done be tuning.";
}

TEST_F(AutoTunerTest, testForceRetuneInPhase) {
  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {2, 4, 2};
  const double cutoff = 1;
  const double cellSizeFactor = 1;
  const double verletSkinPerTimestep = 0;
  const unsigned int verletRebuildFrequency = 20;
  const unsigned int verletClusterSize = 4;
  const double mpiTuningMaxDifferenceForBucket = 0.3;
  const double mpiTuningWeightForMaxDensity = 0.0;
  const unsigned int maxSamples = 3;

  autopas::Configuration confLc_c01(autopas::ContainerOption::linkedCells, cellSizeFactor,
                                    autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
                                    autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled);
  autopas::Configuration confLc_c04(autopas::ContainerOption::linkedCells, cellSizeFactor,
                                    autopas::TraversalOption::lc_c04, autopas::LoadEstimatorOption::none,
                                    autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled);
  autopas::Configuration confLc_c08(autopas::ContainerOption::linkedCells, cellSizeFactor,
                                    autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
                                    autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled);

  auto configsList = {confLc_c01, confLc_c04, confLc_c08};

  auto tuningStrategy = std::make_unique<autopas::FullSearch>(configsList);
  autopas::AutoTuner<Molecule> autoTuner(bBoxMin, bBoxMax, cutoff, verletSkinPerTimestep, verletClusterSize,
                                         std::move(tuningStrategy), mpiTuningMaxDifferenceForBucket,
                                         mpiTuningWeightForMaxDensity, autopas::SelectorStrategyOption::fastestAbs, 100,
                                         maxSamples, verletRebuildFrequency);

  size_t numExpectedTuningIterations = configsList.size() * maxSamples;
  MockFunctor<Molecule> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  std::vector<std::vector<Molecule>> emptyVec(autopas::autopas_get_max_threads());

  // Do part of the tuning phase. After the loop we should be in the middle of sampling the second configuration.
  ASSERT_GT(maxSamples, 1);
  ASSERT_GT(configsList.size(), 1);
  size_t iteration = 0;
  for (; iteration < maxSamples + 1; ++iteration) {
    // since we don't actually do anything doRebuild can always be false.
    EXPECT_TRUE(autoTuner.iteratePairwise(&functor, false, emptyVec, emptyVec)) << "Tuner should still be tuning.\n"
                                                                                   "Phase 1\n"
                                                                                   "Iteration "
                                                                                << iteration;
  }
  // restart the full tuning phase
  autoTuner.forceRetune();
  EXPECT_TRUE(autoTuner.willRebuild()) << "willRebuild() does not recognize forceRetune()";

  // expect a full tuning phase
  for (size_t i = 0; i < numExpectedTuningIterations; ++i, ++iteration) {
    // since we don't actually do anything doRebuild can always be false.
    EXPECT_TRUE(autoTuner.iteratePairwise(&functor, false, emptyVec, emptyVec)) << "Tuner should still be tuning.\n"
                                                                                   "Phase 2\n"
                                                                                   "Iteration "
                                                                                << iteration;
  }
  // first iteration after tuning phase
  EXPECT_FALSE(autoTuner.iteratePairwise(&functor, false, emptyVec, emptyVec)) << "Tuner should be done be tuning.\n"
                                                                                  "Iteration "
                                                                               << iteration;
}

/**
 * Generates no configurations.
 */
TEST_F(AutoTunerTest, testNoConfig) {
  // wrap constructor call into lambda to avoid parser errors
  auto exp1 = []() {
    std::set<autopas::Configuration> configsList = {};
    auto tuningStrategy = std::make_unique<autopas::FullSearch>(configsList);
    autopas::AutoTuner<Molecule> autoTuner({0, 0, 0}, {10, 10, 10}, 0.05, 0, 64, std::move(tuningStrategy), 0.3, 0.0,
                                           autopas::SelectorStrategyOption::fastestAbs, 1000, 3, 20);
  };

  EXPECT_THROW(exp1(), autopas::utils::ExceptionHandler::AutoPasException) << "Constructor with given configs";

  // wrap constructor call into lambda to avoid parser errors
  auto exp2 = []() {
    std::set<autopas::ContainerOption> co = {};
    std::set<double> csf = {};
    std::set<autopas::TraversalOption> tr = {};
    std::set<autopas::LoadEstimatorOption> le = {};
    std::set<autopas::DataLayoutOption> dl = {};
    std::set<autopas::Newton3Option> n3 = {};
    auto tuningStrategy = std::make_unique<autopas::FullSearch>(co, csf, tr, le, dl, n3);
    autopas::AutoTuner<Molecule> autoTuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy), 0.3, 0.0,
                                           autopas::SelectorStrategyOption::fastestAbs, 1000, 3, 20);
  };

  EXPECT_THROW(exp2(), autopas::utils::ExceptionHandler::AutoPasException) << "Constructor which generates configs";
}

/**
 * Generates exactly one valid configuration.
 */
TEST_F(AutoTunerTest, testOneConfig) {
  auto configsList = {_confLc_c08};
  auto tuningStrategy = std::make_unique<autopas::FullSearch>(configsList);
  size_t maxSamples = 3;
  autopas::AutoTuner<Molecule> tuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy), 0.3, 0.0,
                                     autopas::SelectorStrategyOption::fastestAbs, 1000, maxSamples, 20);

  EXPECT_EQ(_confLc_c08, tuner.getCurrentConfig());

  MockFunctor<Molecule> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));

  std::vector<std::vector<Molecule>> emptyVec(autopas::autopas_get_max_threads());

  bool doRebuild = true;
  size_t numSamples = 0;
  for (int i = 0; i < 5; ++i) {
    if (numSamples == maxSamples) {
      numSamples = 0;
      doRebuild = true;
    }
    tuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);
    doRebuild = false;
    ++numSamples;
    EXPECT_EQ(_confLc_c08, tuner.getCurrentConfig());
  }
}

/**
 * Generates exactly one valid and one invalid configuration.
 */
TEST_F(AutoTunerTest, testConfigSecondInvalid) {
  double cellSizeFactor = 1.;
  autopas::Configuration confN3(autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
                                autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                                autopas::Newton3Option::enabled);
  autopas::Configuration confNoN3(autopas::ContainerOption::linkedCells, cellSizeFactor,
                                  autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
                                  autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled);

  auto configsList = {confNoN3, confN3};
  auto tuningStrategy = std::make_unique<autopas::FullSearch>(configsList);
  autopas::AutoTuner<Molecule> tuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy), 0.3, 0.0,
                                     autopas::SelectorStrategyOption::fastestAbs, 1000, 3, 20);

  EXPECT_EQ(confNoN3, tuner.getCurrentConfig());

  std::vector<std::vector<Molecule>> emptyVec(autopas::autopas_get_max_threads());

  MockFunctor<Molecule> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(false));
  bool doRebuild = true;
  tuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);
  EXPECT_EQ(confN3, tuner.getCurrentConfig());
  doRebuild = false;
  tuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);
  EXPECT_EQ(confN3, tuner.getCurrentConfig());
  doRebuild = false;
  tuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);
  EXPECT_EQ(confN3, tuner.getCurrentConfig());
}

/**
 * All generated configurations are thrown out at runtime.
 */
TEST_F(AutoTunerTest, testLastConfigThrownOut) {
  double cellSizeFactor = 1.;
  autopas::Configuration confN3(autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
                                autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                                autopas::Newton3Option::enabled);
  autopas::Configuration confNoN3(autopas::ContainerOption::linkedCells, cellSizeFactor,
                                  autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
                                  autopas::DataLayoutOption::soa, autopas::Newton3Option::enabled);

  auto configsList = {confN3, confNoN3};
  auto tuningStrategy = std::make_unique<autopas::FullSearch>(configsList);
  autopas::AutoTuner<Molecule> tuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy), 0.3, 0.0,
                                     autopas::SelectorStrategyOption::fastestAbs, 1000, 3, 20);

  EXPECT_EQ(confN3, tuner.getCurrentConfig());

  MockFunctor<Molecule> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(false));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  std::vector<std::vector<Molecule>> emptyVec(autopas::autopas_get_max_threads());

  bool doRebuild = true;
  EXPECT_THROW(tuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec),
               autopas::utils::ExceptionHandler::AutoPasException);
}

/**
 * Iterate with two configs.
 * First has short rebuild and long non-rebuild iterations
 * Second has long rebuild and short non-rebuild iterations
 * Expect to choose the first because the second one is worse on average.
 */
TEST_F(AutoTunerTest, testBuildNotBuildTimeEstimation) {
  const double cellSizeFactor = 1.;
  autopas::Configuration confA(autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
                               autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                               autopas::Newton3Option::enabled);
  autopas::Configuration confB(autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c18,
                               autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                               autopas::Newton3Option::enabled);

  auto rebuildFrequency = 3;

  auto configsList = {confA, confB};
  auto tuningStrategy = std::make_unique<autopas::FullSearch>(configsList);
  autopas::AutoTuner<Molecule> tuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy), 0.3, 0.0,
                                     autopas::SelectorStrategyOption::fastestAbs, 1000, 3, rebuildFrequency);

  MockFunctor<Molecule> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));

  std::vector<std::vector<Molecule>> emptyVec(autopas::autopas_get_max_threads());
  std::vector<std::vector<Molecule>> twoParticles(autopas::autopas_get_max_threads());
  twoParticles[0] = {Molecule{}, Molecule{}};

  using namespace std::literals;

  bool doRebuild = true;
  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(100ms); }));
  tuner.iteratePairwise(&functor, doRebuild, twoParticles, emptyVec);

  auto firstConfig = tuner.getCurrentConfig();

  doRebuild = false;
  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(30ms); }));
  tuner.iteratePairwise(&functor, doRebuild, twoParticles, emptyVec);
  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(30ms); }));
  tuner.iteratePairwise(&functor, doRebuild, twoParticles, emptyVec);

  // Here, second config will start to be tuned

  doRebuild = true;
  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(300ms); }));
  tuner.iteratePairwise(&functor, doRebuild, twoParticles, emptyVec);

  auto secondConfig = tuner.getCurrentConfig();

  doRebuild = false;
  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(25ms); }));
  tuner.iteratePairwise(&functor, doRebuild, twoParticles, emptyVec);
  EXPECT_CALL(functor, AoSFunctor).WillOnce(::testing::Invoke([]() { std::this_thread::sleep_for(25ms); }));
  tuner.iteratePairwise(&functor, doRebuild, twoParticles, emptyVec);

  // Here, tuning should be finished and first should have been chosen (100 + 2 * 30 = 160 < 350 = 300 + 2 * 25)
  tuner.iteratePairwise(&functor, doRebuild, emptyVec, emptyVec);

  EXPECT_EQ(tuner.getCurrentConfig(), firstConfig);
  EXPECT_NE(tuner.getCurrentConfig(), secondConfig);
}