/**
 * @file AutoTunerTest.cpp
 * @author F. Gratl
 * @date 8/10/18
 */

#include "AutoTunerTest.h"

#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/selectors/AutoTuner.h"
#include "autopas/selectors/tuningStrategy/FullSearch.h"

using ::testing::_;

TEST_F(AutoTunerTest, testAllConfigurations) {
  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 42};
  // adaptive domain size so sliced is always applicable.
  bBoxMax[2] = autopas::autopas_get_max_threads() * 2;
  const double cutoff = 1;
  const double cellSizeFactor = 1;
  const double verletSkin = 0;
  const unsigned int verletClusterSize = 64;
  const unsigned int maxSamples = 2;
  autopas::LJFunctor<Molecule, FMCell> functor(cutoff);
  auto tuningStrategy = std::make_unique<autopas::FullSearch>(
      autopas::ContainerOption::getAllOptions(), std::set<double>({cellSizeFactor}),
      autopas::TraversalOption::getAllOptions(), autopas::LoadEstimatorOption::getAllOptions(),
      autopas::DataLayoutOption::getAllOptions(), autopas::Newton3Option::getAllOptions());
  autopas::AutoTuner<Molecule, FMCell> autoTuner(bBoxMin, bBoxMax, cutoff, verletSkin, verletClusterSize,
                                                 std::move(tuningStrategy), autopas::SelectorStrategyOption::fastestAbs,
                                                 100, maxSamples);

  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::off);
  //  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  bool stillTuning = true;
  auto prevConfig = autopas::Configuration();

  // total number of possible configurations * number of samples + last iteration after tuning
  // number of configs manually counted:
  //
  // Direct Sum:            directSum traversal         (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  // LinkedCells:           lc_c08 traversal            (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        lc_sliced                   (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        balanced-sliced             (AoS <=> SoA, newton3 <=> noNewton3, 2 heuristics)   = 8
  //                        lc_c18                      (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        lc_c01                      (AoS <=> SoA, noNewton3)                             = 2
  //                        c01-combined-SoA            (SoA, noNewton3)                                     = 1
  //                        lc_c04                      (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        lc_c04_combined_SoA         (SoA, newton3 <=> noNewton3)                         = 2
  //                        lc_c04_HCP                  (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  // VerletLists:           verlet-lists                (AoS <=> SoA, noNewton3)                             = 2
  // VerletListsCells:      verlet-sliced               (AoS, newton3 <=> noNewton3)                         = 2
  //                        balanced-verlet-sliced      (AoS, newton3 <=> noNewton3, 3 heuristics)           = 6
  //                        verlet-c18                  (AoS, newton3 <=> noNewton3)                         = 2
  //                        verlet-c01                  (AoS, noNewton3)                                     = 1
  // VerletClusterLists:    verlet-clusters             (AoS <=> SoA, noNewton3)                             = 2
  //                        verlet-clusters-coloring    (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                        verlet-clusters-static      (AoS <=> SoA, noNewton3)                             = 2
  // VarVerletListsAsBuild: var-verlet-lists-as-build   (AoS <=> SoA, newton3 <=> noNewton3)                 = 4
  //                                                                                                    --------
  //                                                                                                          62
  // Additional with cuda
  // Direct Sum:            directSum traversal            (Cuda, newton3 <=> noNewton3)                     = 2
  // LinkedCells:           lc_c01_cuda traversal          (Cuda, newton3 <=> noNewton3)                     = 2
  // VerletClusterCells:    verlet-cluster-cells traversal (Cuda, newton3 <=> noNewton3)                     = 2
  //                                                                                                    --------
  //                                                                                                          68
  //
  // currently disabled:
  // NORMAL:
  //                                                                                                    --------
  // TOTAL:                                                                                                   68
  //
  // CUDA:
  // C01CudaTraversal for enabled N3, see #420                                                                -1
  //                                                                                                    --------
  // TOTAL:                                                                                                   67

#ifndef AUTOPAS_CUDA
  const size_t expectedNumberOfIterations = 62 * maxSamples + 1;
#else
  const size_t expectedNumberOfIterations = 67 * maxSamples + 1;
#endif

  int collectedSamples = 0;
  int iterations = 0;
  bool doRebuild = true;  // defines whether the verlet lists should be rebuild.
  while (stillTuning) {
    if (collectedSamples == maxSamples) {
      collectedSamples = 0;
      doRebuild = true;
    }
    stillTuning = autoTuner.iteratePairwise(&functor, doRebuild);
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
  autopas::AutoTuner<Particle, FPCell> autoTuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy),
                                                 autopas::SelectorStrategyOption::fastestAbs, 1000, 2);

  EXPECT_EQ(*(configs.begin()), autoTuner.getCurrentConfig());

  MockFunctor<Particle, FPCell> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild for first iteration.";
  bool doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild);  // DS NoN3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild);  // DS NoN3
  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because we change config.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild);  // DS N3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild);  // DS N3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because we change config.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild);  // LC NoN3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild);  // LC NoN3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because reached end of tuning phase.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild);  // optimum
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
  autopas::AutoTuner<Particle, FPCell> autoTuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy),
                                                 autopas::SelectorStrategyOption::fastestAbs, 1000, 2);

  EXPECT_EQ(*(configs.begin()), autoTuner.getCurrentConfig());

  MockFunctor<Particle, FPCell> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(false));

  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild for first iteration.";
  bool doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild);  // DS N3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild);  // DS N3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because we change config.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild);  // LC N3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild);  // LC N3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because reached end of tuning phase.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild);  // optimum
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
  autopas::AutoTuner<Particle, FPCell> autoTuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy),
                                                 autopas::SelectorStrategyOption::fastestAbs, 1000, 2);

  EXPECT_EQ(*(configs.begin()), autoTuner.getCurrentConfig());

  MockFunctor<Particle, FPCell> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  // Intended false positive
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild for first iteration.";
  bool doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild);  // DS NoN3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild);  // DS NoN3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because we change config.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild);  // LC NoN3
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because more samples needed.";
  doRebuild = false;
  autoTuner.iteratePairwise(&functor, doRebuild);  // LC NoN3
  EXPECT_TRUE(autoTuner.willRebuild()) << "Expect rebuild because reached end of tuning phase.";
  doRebuild = true;
  autoTuner.iteratePairwise(&functor, doRebuild);  // optimum
  EXPECT_FALSE(autoTuner.willRebuild()) << "Expect no rebuild because not tuning.";
}

/**
 * Generates no configurations.
 */
TEST_F(AutoTunerTest, testNoConfig) {
  // wrap constructor call into lambda to avoid parser errors
  auto exp1 = []() {
    std::set<autopas::Configuration> configsList = {};
    auto tuningStrategy = std::make_unique<autopas::FullSearch>(configsList);
    autopas::AutoTuner<Particle, FPCell> autoTuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy),
                                                   autopas::SelectorStrategyOption::fastestAbs, 1000, 3);
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
    autopas::AutoTuner<Particle, FPCell> autoTuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy),
                                                   autopas::SelectorStrategyOption::fastestAbs, 1000, 3);
  };

  EXPECT_THROW(exp2(), autopas::utils::ExceptionHandler::AutoPasException) << "Constructor which generates configs";
}

/**
 * Generates exactly one valid configuration.
 */
TEST_F(AutoTunerTest, testOneConfig) {
  autopas::Configuration conf(autopas::ContainerOption::linkedCells, 1., autopas::TraversalOption::lc_c08,
                              autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos,
                              autopas::Newton3Option::enabled);

  auto configsList = {conf};
  auto tuningStrategy = std::make_unique<autopas::FullSearch>(configsList);
  size_t maxSamples = 3;
  autopas::AutoTuner<Particle, FPCell> tuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy),
                                             autopas::SelectorStrategyOption::fastestAbs, 1000, maxSamples);

  EXPECT_EQ(conf, tuner.getCurrentConfig());

  MFunctor functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));

  bool doRebuild = true;
  size_t numSamples = 0;
  for (int i = 0; i < 5; ++i) {
    if (numSamples == maxSamples) {
      numSamples = 0;
      doRebuild = true;
    }
    tuner.iteratePairwise(&functor, doRebuild);
    doRebuild = false;
    ++numSamples;
    EXPECT_EQ(conf, tuner.getCurrentConfig());
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
  autopas::AutoTuner<Particle, FPCell> tuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy),
                                             autopas::SelectorStrategyOption::fastestAbs, 1000, 3);

  EXPECT_EQ(confNoN3, tuner.getCurrentConfig());

  MFunctor functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(false));
  bool doRebuild = true;
  tuner.iteratePairwise(&functor, doRebuild);
  EXPECT_EQ(confN3, tuner.getCurrentConfig());
  doRebuild = false;
  tuner.iteratePairwise(&functor, doRebuild);
  EXPECT_EQ(confN3, tuner.getCurrentConfig());
  doRebuild = false;
  tuner.iteratePairwise(&functor, doRebuild);
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
  autopas::AutoTuner<Particle, FPCell> tuner({0, 0, 0}, {10, 10, 10}, 1, 0, 64, std::move(tuningStrategy),
                                             autopas::SelectorStrategyOption::fastestAbs, 1000, 3);

  EXPECT_EQ(confN3, tuner.getCurrentConfig());

  MFunctor functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(false));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));

  bool doRebuild = true;
  EXPECT_THROW(tuner.iteratePairwise(&functor, doRebuild), autopas::utils::ExceptionHandler::AutoPasException);
}
