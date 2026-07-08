/**
 * @file TuningManagerTest.cpp
 * @author muehlhaeusser
 * @date 18.04.2026
 */

#include "TuningManagerTest.h"

#include <ranges>

#include "autopas/LogicHandler.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/TuningManager.h"
#include "autopas/tuning/tuningStrategy/RandomSearch.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyFactory.h"
#include "autopas/tuning/tuningStrategy/fuzzyTuning/FuzzyTuning.h"
#include "autopas/tuning/utils/AutoTunerInfo.h"
#include "autopas/tuning/utils/SearchSpaceGenerators.h"
#include "testingHelpers/GenerateValidConfigurations.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;

/**
 * Test that the TuningManager strictly adheres to the fixed tuning interval.
 */
TEST_F(TuningManagerTest, testTuningIntervalIsFixed) {
  constexpr unsigned int verletRebuildFrequency = 20;
  constexpr unsigned int tuningInterval = 1000;
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = tuningInterval,
      .maxSamples = 3,
  };
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const autopas::AutoTuner::SearchSpaceType searchSpace{
      _confDs_seq_noN3,
      _confDs_seq_N3,
      _confLc_c08_noN3,
  };

  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, verletRebuildFrequency, "");

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));

  constexpr size_t numIterations = 5 * tuningInterval;
  bool lastWasTuningInteration = false;
  size_t iterations = 0;

  while (iterations < numIterations) {
    [[maybe_unused]] auto emigrants = logicHandler.updateContainer();
    const auto stillTuning =
        logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);

    if (not lastWasTuningInteration and stillTuning) {
      EXPECT_TRUE(iterations % tuningInterval == 0)
          << "Tuning phase did not start at fixed iteration number. Iteration: " << iterations;
    }
    ++iterations;
    lastWasTuningInteration = stillTuning;
  }
}

/**
 * Test tuning with two autotuners in combination by checking the return value of `bool stillTuning =
 * LogicHandler::computeInteractionsPipeline()`. The tuners have differently sized search spaces and therefore finish
 * tuning at a different time step. Yet, they should start the second tuning phase in the same time step.
 */
TEST_F(TuningManagerTest, testMultipleTuners) {
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  constexpr size_t rebuildFrequency = 3;
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 10,
      .maxSamples = 2,
  };
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };

  const auto pairwiseSearchSpace = {_confDs_seq_N3, _confLc_c18_N3, _confLc_c08_N3};
  const auto triwiseSearchSpace = {_confDs_3b_N3, _confLc_c01_3b_noN3};

  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(std::make_unique<autopas::AutoTuner>(tuningStrategies, pairwiseSearchSpace, autoTunerInfo,
                                                                   rebuildFrequency, "2B"),
                              autopas::InteractionTypeOption::pairwise);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, triwiseSearchSpace, autoTunerInfo, rebuildFrequency, "3B"),
      autopas::InteractionTypeOption::triwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, rebuildFrequency, "");

  testing::NiceMock<MockPairwiseFunctor<Molecule>> pairFunctor;
  testing::NiceMock<MockTriwiseFunctor<Molecule>> triFunctor;

  EXPECT_CALL(pairFunctor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(triFunctor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(pairFunctor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(triFunctor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(pairFunctor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(triFunctor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(pairFunctor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(triFunctor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));

  // Add three particles into one (linked cells) cell
  logicHandler.getContainer().addParticle((Molecule{{0.1, 0.1, 0.1}, {0., 0., 0.}, 0, 0}));
  logicHandler.getContainer().addParticle((Molecule{{0.2, 0.1, 0.1}, {0., 0., 0.}, 1, 0}));
  logicHandler.getContainer().addParticle((Molecule{{0.1, 0.2, 0.1}, {0., 0., 0.}, 2, 0}));

  // Beginning of the first tuning phase. Both tuners are still tuning.
  EXPECT_CALL(pairFunctor, AoSFunctor).Times(4 * 3);  // 3 AoS functor calls per timestep (All configs with N3)
  EXPECT_CALL(triFunctor, AoSFunctor).Times(4 * 3);   // 3 AoS functor calls per timestep (All configs without N3)
  for (int i = 0; i < 4; i++) {
    logicHandler.updateContainer();
    EXPECT_TRUE(logicHandler.computeInteractionsPipeline(&pairFunctor, autopas::InteractionTypeOption::pairwise));
    EXPECT_TRUE(logicHandler.computeInteractionsPipeline(&triFunctor, autopas::InteractionTypeOption::triwise));
  }

  // End of first tuning phase. The triwise tuner has finished tuning and is waiting for the pairwise tuner.
  EXPECT_CALL(pairFunctor, AoSFunctor).Times(2 * 3);
  EXPECT_CALL(triFunctor, AoSFunctor).Times(2 * 3);
  for (int i = 0; i < 2; i++) {
    logicHandler.updateContainer();
    EXPECT_TRUE(logicHandler.computeInteractionsPipeline(&pairFunctor, autopas::InteractionTypeOption::pairwise));
    EXPECT_TRUE(logicHandler.computeInteractionsPipeline(&triFunctor, autopas::InteractionTypeOption::triwise));
  }

  // Outside the tuning phase. Both tuners run with their best configuration. 4 more iterations until next tuning phase.
  EXPECT_CALL(pairFunctor, AoSFunctor).Times(4 * 3);
  EXPECT_CALL(triFunctor, AoSFunctor).Times(4 * 3);
  for (int i = 0; i < 4; i++) {
    logicHandler.updateContainer();
    EXPECT_FALSE(logicHandler.computeInteractionsPipeline(&pairFunctor, autopas::InteractionTypeOption::pairwise));
    EXPECT_FALSE(logicHandler.computeInteractionsPipeline(&triFunctor, autopas::InteractionTypeOption::triwise));
  }

  // Beginning of the second tuning phase. Both tuners start tuning again at the same time.
  EXPECT_CALL(pairFunctor, AoSFunctor).Times(4 * 3);
  EXPECT_CALL(triFunctor, AoSFunctor).Times(4 * 3);
  for (int i = 0; i < 4; i++) {
    logicHandler.updateContainer();
    EXPECT_TRUE(logicHandler.computeInteractionsPipeline(&pairFunctor, autopas::InteractionTypeOption::pairwise));
    EXPECT_TRUE(logicHandler.computeInteractionsPipeline(&triFunctor, autopas::InteractionTypeOption::triwise));
  }
}

/**
 * Test the behavior when the tuning phase takes more iterations than the specified tuning interval.
 */
TEST_F(TuningManagerTest, testTuningPhaseLongerThanTuningInterval) {
  constexpr unsigned int verletRebuildFrequency = 15;
  constexpr unsigned int tuningInterval = 10;
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = tuningInterval,
      .maxSamples = 6,
  };
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const autopas::AutoTuner::SearchSpaceType searchSpace{
      _confLc_c18_noN3,
      _confLc_c08_N3,
  };

  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, verletRebuildFrequency, "");

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));

  constexpr size_t iterationsToDo = 35;
  size_t iterationsDone = 0;

  while (iterationsDone < iterationsToDo) {
    [[maybe_unused]] auto emigrants = logicHandler.updateContainer();
    const bool stillTuning =
        logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);

    if (iterationsDone < searchSpace.size() * autoTunerInfo.maxSamples) {
      EXPECT_TRUE(stillTuning) << "AutoTuner should tune.";
    }
    if (iterationsDone >= searchSpace.size() * autoTunerInfo.maxSamples and iterationsDone < 2 * tuningInterval) {
      EXPECT_FALSE(stillTuning) << "AutoTuner should not tune.";
    }
    if (iterationsDone >= 2 * tuningInterval and
        iterationsDone < 2 * tuningInterval + searchSpace.size() * autoTunerInfo.maxSamples) {
      EXPECT_TRUE(stillTuning) << "AutoTuner should tune.";
    }
    if (iterationsDone >= 2 * tuningInterval + searchSpace.size() * autoTunerInfo.maxSamples) {
      EXPECT_FALSE(stillTuning) << "AutoTuner should not tune.";
    }
    ++iterationsDone;
  }
}

/**
 * Initialize a tuner, do a tuning phase, then instead of waiting for the full tuning interval restart the tuning
 * earlier via forceRetune.
 */
TEST_F(TuningManagerTest, testForceRetuneBetweenPhases) {
  constexpr unsigned int verletRebuildFrequency = 20;
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{2, 4, 2},
  };
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };

  autopas::AutoTuner::SearchSpaceType searchSpace{_confLc_c01_noN3, _confLc_c18_noN3, _confLc_c08_noN3};
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, verletRebuildFrequency, "");

  const size_t numExpectedTuningIterations = searchSpace.size() * autoTunerInfo.maxSamples;
  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));

  size_t iteration = 0;
  // 1. Expect a full tuning phase
  for (size_t i = 0; i < numExpectedTuningIterations; ++i, ++iteration) {
    logicHandler.updateContainer();
    EXPECT_TRUE((logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise)))
        << "Tuner should still be tuning in iteration " << iteration;
  }

  // 2. First iteration after tuning phase
  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuilding after tuning finished.";
  EXPECT_FALSE((logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise)))
      << "Tuner should be done tuning in the first iteration after the tuning phase.";
  ++iteration;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "No rebuilding expected here.";

  // Advance an extra iteration so we are clearly out of the tuning start boundary
  logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
  ++iteration;

  // 3. Instead of waiting the full tuning interval, restart tuning immediately
  tuningManager->forceRetune();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration))
      << "TuningManager does not recognize forceRetune() required rebuild";

  // 4. Expect a second full tuning phase
  for (size_t i = 0; i < numExpectedTuningIterations; ++i, ++iteration) {
    logicHandler.updateContainer();
    EXPECT_TRUE((logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise)))
        << "Tuner should still be tuning after a forced retune in iteration " << iteration;
  }

  // 5. First iteration after the forced tuning phase
  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration))
      << "Expect rebuilding after the forced tuning phase finished.";
  EXPECT_FALSE((logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise)))
      << "Tuner should be done tuning after the forced phase.";
}

/**
 * Test the correct AutoTuner behavior if a retuning is forced during a tuning phase.
 */
TEST_F(TuningManagerTest, testForceRetuneInPhase) {
  constexpr unsigned int verletRebuildFrequency = 20;
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{2, 4, 2},
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const auto searchSpace = {_confLc_c01_noN3, _confLc_c18_noN3, _confLc_c08_noN3};

  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, verletRebuildFrequency, "");

  const size_t numExpectedTuningIterations = searchSpace.size() * autoTunerInfo.maxSamples;
  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));

  // Do part of the tuning phase. After the loop we should be in the middle of sampling the second configuration.
  ASSERT_GT(autoTunerInfo.maxSamples, 1);
  ASSERT_GT(searchSpace.size(), 1);
  size_t iteration = 0;
  for (; iteration < autoTunerInfo.maxSamples + 1; ++iteration) {
    logicHandler.updateContainer();
    EXPECT_TRUE((logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise)))
        << "Tuner should still be tuning.\nPhase 1, Iteration " << iteration;
  }

  // restart the full tuning phase
  tuningManager->forceRetune();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "TuningManager does not recognize forceRetune()";

  // expect a full tuning phase
  for (size_t i = 0; i < numExpectedTuningIterations; ++i, ++iteration) {
    logicHandler.updateContainer();
    EXPECT_TRUE((logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise)))
        << "Tuner should still be tuning.\nPhase 2, Iteration " << iteration;
  }

  // first iteration after tuning phase
  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration))
      << "Expect rebuilding after the forced tuning phase finished.";
  EXPECT_FALSE((logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise)))
      << "Tuner should be done be tuning.\nIteration " << iteration;
}

/**
 * Test the AutoTuner with all valid pairwise configurations in the search space.
 */
TEST_F(TuningManagerTest, testAllConfigurations) {
  const autopas::NumberSetFinite<double> cellSizeFactors({1});
  constexpr unsigned int rebuildFrequency = 20;
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{8., 8., 8.},
      .cutoff = 1.8,
      .verletSkin = .2,
  };
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = std::numeric_limits<unsigned int>::max(),
      .maxSamples = 2,
  };

  // the NiceMock wrapper suppresses warnings from uninteresting function calls
  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));
  // Need to resize cells during loading, otherwise we get exceptions in SoAFunctors
  EXPECT_CALL(functor, SoALoader(::testing::Matcher<autopas::ReferenceParticleCell<Molecule> &>(_), _, _, _))
      .Times(testing::AtLeast(1))
      .WillRepeatedly(
          testing::WithArgs<0, 1>(testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.size()); })));
  EXPECT_CALL(functor, SoALoader(::testing::Matcher<FMCell &>(_), _, _, _))
      .Times(testing::AtLeast(1))
      .WillRepeatedly(
          testing::WithArgs<0, 1>(testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.size()); })));

  const auto searchSpace = generateAllValidConfigurations(autopas::InteractionTypeOption::pairwise);
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);
  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, rebuildFrequency, "");
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::off);

  std::map<autopas::ContainerOption, size_t> configsPerContainer;

  // Number of configs manually counted:
  // Configurations considered valid by cartesian product but are typically invalid due to
  // Traversal::isApplicableToDomain are still included in this count.
  //
  // Direct Sum:            ds_sequential               (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //   Subtotal:                                                                                              =   4
  // Direct Sum only supports CSF 1 => Multiply by 1
  // Direct Sum supports all 4 vectorization patterns => Multiply by 4

  configsPerContainer[autopas::ContainerOption::directSum] = 16;
  // LinkedCells:           lc_c08                      (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        lc_sliced                   (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        lc_sliced_balanced          (AoS <=> SoA, newton3 <=> noNewton3, 2 heuristics)    =   8
  //                        lc_sliced_c02               (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        lc_c18                      (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        lc_c01                      (AoS <=> SoA, noNewton3)                              =   2
  //                        lc_c01_combined_SoA         (SoA, noNewton3)                                      =   1
  //                        lc_c04                      (AoS <=> SoA, newton3 <=> noNewton3) *                =   4
  //                        lc_c04_combined_SoA         (SoA, newton3 <=> noNewton3)         *                =   2
  //                        lc_c04_HCP                  (AoS <=> SoA, newton3 <=> noNewton3) *                =   4
  //   Subtotal:                                                                                              =  37
  // Linked Cells supports all CSFs (i.e. 0.5, 1.0, 1.5) => Multiply by 3
  // * C04 traversals will almost always not support CSF 0.5, but these will be excluded later.
  // Linked Cells supports all 4 vectorization patterns => Multiply by 4
  configsPerContainer[autopas::ContainerOption::linkedCells] = 444;

  // Linked Cells References:
  // same as linked Cells but without the squaredParticlesPerCell load estimator for sliced_balanced (-4 * 3 * 4 = -48)
  configsPerContainer[autopas::ContainerOption::linkedCellsReferences] =
      configsPerContainer[autopas::ContainerOption::linkedCells] - 48;

  // VerletLists:           vl_list_iteration           (AoS <=> SoA, noNewton3)                              =   2
  //   Subtotal:                                                                                              =   2
  // Verlet Lists only supports CSF >= 1 (i.e. 1.0, 1.5) => Multiply by 2
  // Verlet Lists only supports 1xVec pattern => Multiply by 1
  configsPerContainer[autopas::ContainerOption::verletLists] = 4;

  // VerletListsCells:      vlc_sliced                  (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        vlc_sliced_balanced         (AoS <=> SoA, newton3 <=> noNewton3, 3 LB heuristics) =  12
  //                        vlc_sliced_colored          (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        vlc_c18                     (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        vlc_c01                     (AoS <=> SoA, noNewton3)                              =   2
  //                        vlc_c08                     (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //   Subtotal:                                                                                              =  30
  // Verlet Lists Cells only supports CSF >= 1 (i.e. 1.0, 1.5) => Multiply by 2
  // Verlet Lists Cells only supports 1xVec pattern => Multiply by 1
  configsPerContainer[autopas::ContainerOption::verletListsCells] = 60;

  // VerletClusterLists:    vcl_cluster_iteration       (AoS <=> SoA, noNewton3)                              =   2
  //                        vcl_c06                     (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        vcl_c01_balanced            (AoS <=> SoA, noNewton3)                              =   2
  //                        vcl_sliced                  (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        vcl_sliced_c02              (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        vcl_sliced_balanced         (AoS <=> SoA, newton3 <=> noNewton3, 2 LB heuristics) =   8
  //   Subtotal:                                                                                              =  24
  // Verlet Cluster Lists only support CSF 1 => Multiply by 1
  // Verlet Cluster Lists supports all 4 vectorization patterns => Multiply by 4
  configsPerContainer[autopas::ContainerOption::verletClusterLists] = 96;

  // VarVerletListsAsBuild: vvl_as_built                (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //   Subtotal:                                                                                              =   4
  // Var Verlet Lists only supports CSF >= 1 (i.e. 1.0, 1.5) => Multiply by 2
  // Var Verlet Lists only supports 1xVec pattern => Multiply by 1
  configsPerContainer[autopas::ContainerOption::varVerletListsAsBuild] = 8;

  // PairwiseVerletLists:   vlp_sliced                  (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        vlp_sliced_balanced         (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        vlp_sliced_colored          (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        vlp_c18                     (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //                        vlp_c01                     (AoS <=> SoA, noNewton3)                              =   2
  //                        vlp_c08                     (AoS <=> SoA, newton3 <=> noNewton3)                  =   4
  //   Subtotal:                                                                                              =  22
  // Pairwise Verlet Lists only supports CSF >= 1 (i.e. 1.0, 1.5) => Multiply by 2
  // Pairwise Verlet Lists only supports 1xVec pattern => Multiply by 1
  configsPerContainer[autopas::ContainerOption::pairwiseVerletLists] = 44;

  // PseudoVerletLists:     psvl_c01                    (AoS only, noNewton3)                                 =   1
  //                        psvl_c08                    (AoS only, newton3 <=> noNewton3)                     =   2
  //                        psvl_c18                    (AoS only, newton3 <=> noNewton3)                     =   2
  //   Subtotal:                                                                                              =   5
  // Pseudo Verlet Lists only supports CSF >= 1 (i.e. 1.0, 1.5) => Multiply 2
  // Pseudo Verlet Lists only supports 1xVec pattern => Multiply by 1
  configsPerContainer[autopas::ContainerOption::pseudoVerletLists] = 10;

  // Octree:                ot_c01                      (AoS <=> SoA, noNewton3)                              =   2
  //                        ot_c18                      (AoS <=> SoA, newton3)                                =   2
  //   Subtotal:                                                                                              =   4
  // Octree only supports CSF 1 => Multiply by 1
  // Octree supports all 4 vectorization patterns => Multiply by 4
  configsPerContainer[autopas::ContainerOption::octree] = 16;

  // --------------- Check that the manually determined values above match that automatically generated ---------------

  // check that there is an entry for every container.
  ASSERT_EQ(configsPerContainer.size(), autopas::ContainerOption::getAllOptions().size());

  // Check that there are the correct number of configurations per container, and print detailed breakdowns if not.
  for (const auto &container : autopas::ContainerOption::getAllOptions()) {
    const auto expectedNumConfigs = configsPerContainer.at(container);

    const auto actualNumConfigs =
        std::ranges::count_if(searchSpace, [container](const auto &config) { return config.container == container; });

    auto detailedConfigsPerTraversal = [&]() {
      // Filter configs for this container and group by traversal
      auto containerConfigs =
          searchSpace | std::views::filter([container](const auto &config) { return config.container == container; });

      // Count configurations per traversal
      std::map<autopas::TraversalOption, size_t> traversalCounts;
      std::ranges::for_each(containerConfigs,
                            [&traversalCounts](const auto &config) { ++traversalCounts[config.traversal]; });

      // Build error message
      std::ostringstream oss;
      oss << "\n Number of configurations per traversal for this container:";
      std::ranges::for_each(traversalCounts,
                            [&](const auto &pair) { oss << "\n  " << pair.first.to_string() << ": " << pair.second; });
      return oss.str();
    };

    EXPECT_EQ(actualNumConfigs, expectedNumConfigs)
        << "Unexpected number of configurations for container " << container.to_string()
        << ".\n Check that the manually determined values in AutoTunerTest::testAllConfigurations are correct!"
        << detailedConfigsPerTraversal();
  }

  // Check that the total number of configurations is correct and stop if not
  const size_t numberOfConfigs = std::accumulate(configsPerContainer.begin(), configsPerContainer.end(), 0ul,
                                                 [](auto acc, auto &pair) { return acc + pair.second; });
  ASSERT_EQ(numberOfConfigs, searchSpace.size())
      << "The calculated number of configurations is not equal to the cross product search space!";
  // total number of possible configurations * number of samples + last iteration after tuning

  // Manually determine how many configurations will not be applicable to the domain.
  // As the domain is rather ordinary, this probably will only result from configurations with traversals not supporting
  // CSF < 1.0, but where the container in general could support this.

  // All LC C04 traversals with CSF 0.5 are expected to not be applicable
  // => lc_c04 (4), lc_c04_combined_SoA (2), lc_c04_HCP (4) x4 vecPatterns => 40 configurations
  // Similarly for Linked Cells References => 40 more

  constexpr size_t numConfigsExpectedNotApplicable{80};

  const size_t expectedNumberOfIterations =
      (numberOfConfigs - numConfigsExpectedNotApplicable) * autoTunerInfo.maxSamples + 1;

  ASSERT_GT(autoTunerInfo.tuningInterval, expectedNumberOfIterations)
      << "The tuning interval must be greater than the expected number of iterations! Fix this test.";

  size_t collectedSamples = 0;
  size_t iterations = 0;
  bool stillTuning = true;
  auto prevConfig = autopas::Configuration();

  // We should have the loop end at iterations == expectedNumberOfIterations, but for ease of debugging we extend this
  // to double. We don't use a while(stillTuning) to avoid infinite looping.
  for (; iterations < 2 * expectedNumberOfIterations; ++iterations) {
    if (collectedSamples == autoTunerInfo.maxSamples) {
      collectedSamples = 0;
      logicHandler.getContainer().deleteAllParticles();
    }

    // Should not have any leaving particles in this test
    auto leavingParticlesVec = logicHandler.updateContainer();
    EXPECT_EQ(leavingParticlesVec.size(), 0)
        << "A non-zero number of particles have left the container when none should!";
    stillTuning = logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
    ++collectedSamples;

    const auto currentConfig = tuningManager->getCurrentConfig(autopas::InteractionTypeOption::pairwise);
    if (stillTuning) {
      if (collectedSamples == 1) {
        EXPECT_NE(currentConfig, prevConfig)
            << "Iteration: " << iterations << " collectedSamples: " << collectedSamples;
      } else {
        EXPECT_EQ(currentConfig, prevConfig)
            << "Iteration: " << iterations << " collectedSamples: " << collectedSamples;
      }
    } else {
      // Increment iteration counter to "complete" this iteration
      ++iterations;
      break;
    }
    prevConfig = currentConfig;
  }

  EXPECT_EQ(expectedNumberOfIterations, iterations);
}

/**
 * This tests the rebuild logic of the AutoTuner for two configurations. (Direct Sum and Linked Cells)
 */
TEST_F(TuningManagerTest, testWillRebuildDL) {
  // also check if rebuild is detected if next config is invalid
  constexpr unsigned int rebuildFrequency = 20;
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 2,
  };
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  const autopas::AutoTuner::SearchSpaceType searchSpace{
      _confDs_seq_noN3,
      _confLc_c08_noN3,
  };

  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, rebuildFrequency, "");

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));

  size_t iteration = 0;
  // Intended false positive
  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuild in first iteration.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // DS NoN3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // DS NoN3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuild because we change config.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // LC NoN3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // LC NoN3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuild because reached end of tuning phase.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // optimum
  iteration++;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "Expect no rebuild because not tuning.";
}

/**
 * Tests if AutoPas will correctly rebuild for new configurations (DS NoN3, DS N3, LC NoN3)
 */
TEST_F(TuningManagerTest, testWillRebuildDDL) {
  // also check if rebuild is detected if next config is invalid
  constexpr unsigned int rebuildFrequency = 20;
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 2,
  };
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const autopas::AutoTuner::SearchSpaceType searchSpace{
      _confDs_seq_noN3,
      _confDs_seq_N3,
      _confLc_c08_noN3,
  };

  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, rebuildFrequency, "");

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));

  // updateContainer increments the logic handler's iteration counters for this time step and checks the rebuild
  // conditions afterwards.

  size_t iteration = 0;
  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuild for first iteration.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // DS NoN3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // DS NoN3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuild because we change config.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // DS N3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // DS N3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuild because we change config.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // LC NoN3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // LC NoN3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuild because reached end of tuning phase.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // optimum
  iteration++;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "Expect no rebuild because not tuning.";
}

/**
 * This test simulates that the next config (which is checked by willRebuildNeighborLists) is invalid.
 */
TEST_F(TuningManagerTest, testWillRebuildDDLOneConfigKicked) {
  // also check if rebuild is detected if next config is invalid

  constexpr unsigned int rebuildFrequency = 20;
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 2,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const autopas::AutoTuner::SearchSpaceType searchSpace{
      _confDs_seq_noN3,
      _confDs_seq_N3,
      _confLc_c08_N3,
  };

  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, rebuildFrequency, "");

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(false));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));

  size_t iteration = 0;
  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuild for first iteration.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // DS N3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // DS N3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuild because we change config.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // LC N3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "Expect no rebuild because more samples needed.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // LC N3
  iteration++;

  logicHandler.updateContainer();
  EXPECT_TRUE(tuningManager->requiresRebuilding(iteration)) << "Expect rebuild because reached end of tuning phase.";
  logicHandler.computeInteractionsPipeline(&functor,
                                           autopas::InteractionTypeOption::pairwise);  // optimum
  iteration++;

  logicHandler.updateContainer();
  EXPECT_FALSE(tuningManager->requiresRebuilding(iteration)) << "Expect no rebuild because not tuning.";
}

/**
 * Generates exactly one valid configuration.
 */
TEST_F(TuningManagerTest, testOneConfig) {
  constexpr unsigned int rebuildFrequency = 20;
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  const auto searchSpace = {_confLc_c08_noN3};
  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, rebuildFrequency, "");

  EXPECT_EQ(_confLc_c08_noN3, tuningManager->getCurrentConfig(autopas::InteractionTypeOption::pairwise));

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));

  size_t numSamples = 0;
  for (int i = 0; i < 5; ++i) {
    if (numSamples == autoTunerInfo.maxSamples) {
      numSamples = 0;
    }
    logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
    ++numSamples;
    EXPECT_EQ(_confLc_c08_noN3, tuningManager->getCurrentConfig(autopas::InteractionTypeOption::pairwise));
  }
}

/**
 * Generates exactly one valid and one invalid configuration.
 */
TEST_F(TuningManagerTest, testConfigSecondInvalid) {
  constexpr unsigned int rebuildFrequency = 20;
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  const auto searchSpace = {_confLc_c08_noN3, _confLc_c08_N3};
  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, rebuildFrequency, "");

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;

  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(false));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));

  logicHandler.updateContainer();
  logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
  EXPECT_EQ(_confLc_c08_N3, tuningManager->getCurrentConfig(autopas::InteractionTypeOption::pairwise));
  logicHandler.updateContainer();
  logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
  EXPECT_EQ(_confLc_c08_N3, tuningManager->getCurrentConfig(autopas::InteractionTypeOption::pairwise));
  logicHandler.updateContainer();
  logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
  EXPECT_EQ(_confLc_c08_N3, tuningManager->getCurrentConfig(autopas::InteractionTypeOption::pairwise));
}

/**
 * All generated configurations are thrown out at runtime.
 */
TEST_F(TuningManagerTest, testLastConfigThrownOut) {
  constexpr unsigned int rebuildFrequency = 20;
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
  };
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const auto searchSpace = {_confLc_c08_noN3, _confLc_c18_noN3};
  const auto tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, rebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);

  autopas::LogicHandler<Molecule> logicHandler(tuningManager, logicHandlerInfo, rebuildFrequency, "");

  testing::NiceMock<MockPairwiseFunctor<Molecule>> functor;
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(false));

  EXPECT_THROW((logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise)),
               autopas::utils::ExceptionHandler::AutoPasException);
}

/**
 * Test that the TuningManager correctly calculates the cross-tuner optimum
 * and forces the common container at the end of a tuning phase.
 */
TEST_F(TuningManagerTest, testSetOptimalConfigurationsCommonContainer) {
  constexpr unsigned int rebuildFrequency = 10;
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 100,
      .maxSamples = 1,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies1{};
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies2{};

  // Pairwise Search Space
  const autopas::AutoTuner::SearchSpaceType pairwiseSpace{_confLc_c08_noN3, _confDs_seq_noN3};
  // Triwise Search Space
  const autopas::AutoTuner::SearchSpaceType triwiseSpace{_confLc_c01_3b_noN3, _confDs_3b_N3};

  auto tunerManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tunerManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies1, pairwiseSpace, autoTunerInfo, rebuildFrequency, "2B"),
      autopas::InteractionTypeOption::pairwise);
  tunerManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies2, triwiseSpace, autoTunerInfo, rebuildFrequency, "3B"),
      autopas::InteractionTypeOption::triwise);

  const autopas::LiveInfo info{};
  size_t iteration = 0;

  // --- Phase: Sample the first configurations twice(DS) ---
  // Give DS a good rebuild time, but bad traversal time for triwise
  tunerManager->tune(iteration, info);
  tunerManager->addMeasurement(100, 10, true, iteration, autopas::InteractionTypeOption::pairwise);
  tunerManager->addMeasurement(100, 1000, true, iteration, autopas::InteractionTypeOption::triwise);
  iteration++;

  // --- Phase: Sample the second configurations (LC) ---
  // Give LC a worse rebuild, but equally good traversal time
  tunerManager->tune(iteration, info);
  tunerManager->addMeasurement(1000, 100, true, iteration, autopas::InteractionTypeOption::pairwise);
  tunerManager->addMeasurement(1000, 100, true, iteration, autopas::InteractionTypeOption::triwise);
  iteration++;

  // --- Trigger End of Tuning Phase ---
  // The TuningManager evaluates the evidence and finds that the best combination is to stay with Linked Cells
  // Pairwise would prefer DS, but we stay with LC to avoid excessive rebuilding
  tunerManager->tune(iteration, info);

  // Assert both tuners were forcefully locked into the LC container
  EXPECT_EQ(tunerManager->getCurrentConfig(autopas::InteractionTypeOption::pairwise), _confLc_c08_noN3);
  EXPECT_EQ(tunerManager->getCurrentConfig(autopas::InteractionTypeOption::triwise), _confLc_c01_3b_noN3);

  // Assert both tuners correctly exited the tuning phase
  EXPECT_FALSE(tunerManager->getAutoTuners().at(autopas::InteractionTypeOption::pairwise)->inTuningPhase());
  EXPECT_FALSE(tunerManager->getAutoTuners().at(autopas::InteractionTypeOption::triwise)->inTuningPhase());
}

/**
 * Test that the TuningManager correctly calculates the cross-tuner optimum
 * and forces the common container at the end of a tuning phase.
 */
TEST_F(TuningManagerTest, testSetOptimalConfigurationsDifferentContainers) {
  constexpr unsigned int rebuildFrequency = 10;
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 100,
      .maxSamples = 1,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies1{};
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies2{};

  // Pairwise Search Space
  const autopas::AutoTuner::SearchSpaceType pairwiseSpace{_confLc_c08_noN3, _confDs_seq_noN3};
  // Triwise Search Space
  const autopas::AutoTuner::SearchSpaceType triwiseSpace{_confLc_c01_3b_noN3, _confDs_3b_N3};

  const auto tunerManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tunerManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies1, pairwiseSpace, autoTunerInfo, rebuildFrequency, "2B"),
      autopas::InteractionTypeOption::pairwise);
  tunerManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies2, triwiseSpace, autoTunerInfo, rebuildFrequency, "3B"),
      autopas::InteractionTypeOption::triwise);

  const autopas::LiveInfo info{};
  size_t iteration = 0;

  // --- Phase: Sample the first configurations (DS) ---
  tunerManager->tune(iteration, info);  // Tuners transition to DS
  // Let Direct Sum perform good for pairwise, but bad for triwise
  tunerManager->addMeasurement(10, 10, true, iteration, autopas::InteractionTypeOption::pairwise);
  tunerManager->addMeasurement(10, 1000, true, iteration, autopas::InteractionTypeOption::triwise);
  iteration++;

  // --- Phase: Sample the second configurations (LC) ---
  tunerManager->tune(iteration, info);  // Tuners transition to LC
  // Let Linked Cells perform good for triwise, but bad for pairwise
  tunerManager->addMeasurement(10, 1000, true, iteration, autopas::InteractionTypeOption::pairwise);
  tunerManager->addMeasurement(10, 10, true, iteration, autopas::InteractionTypeOption::triwise);
  iteration++;

  // --- Trigger End of Tuning Phase ---
  // The TuningManager should evaluate the evidence and realize that switching containers between every pairwise and
  // triwise functor call is still better, even if we have to rebuild every time.
  tunerManager->tune(iteration, info);

  EXPECT_EQ(tunerManager->getCurrentConfig(autopas::InteractionTypeOption::pairwise), _confDs_seq_noN3);
  EXPECT_EQ(tunerManager->getCurrentConfig(autopas::InteractionTypeOption::triwise), _confLc_c01_3b_noN3);

  // Assert both tuners correctly exited the tuning phase
  EXPECT_FALSE(tunerManager->getAutoTuners().at(autopas::InteractionTypeOption::pairwise)->inTuningPhase());
  EXPECT_FALSE(tunerManager->getAutoTuners().at(autopas::InteractionTypeOption::triwise)->inTuningPhase());
}

/**
 * Test that TuningManager only reports needing LiveInfo exactly when a tuning phase
 * is about to begin or has just started, and blocks it otherwise.
 */
TEST_F(TuningManagerTest, testLiveInfoRouting) {
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 100,
      .maxSamples = 1,
  };

  // We need a tuner that actually requests LiveInfo to test the gatekeeper
  const autopas::AutoTuner::SearchSpaceType searchSpace{_confLc_c08_noN3};

  // Dummy strategy 1: Only needs standard LiveInfo (at the start of a phase)
  class DummyStrategyLiveInfo : public autopas::TuningStrategyInterface {
   public:
    bool needsLiveInfo() const override { return true; }
    bool needsDomainSimilarityStatistics() const override { return false; }
    autopas::TuningStrategyOption getOptionType() const override { return autopas::TuningStrategyOption::fullSearch; }
    bool optimizeSuggestions(std::vector<autopas::Configuration> &configQueue,
                             const autopas::EvidenceCollection &evidenceCollection) override {
      return true;
    }
    bool reset(size_t iteration, size_t tuningPhase, std::vector<autopas::Configuration> &configQueue,
               const autopas::EvidenceCollection &evidenceCollection) override {
      return true;
    }
  };

  // Dummy strategy 2: Only needs Domain Similarity Statistics (before a phase begins)
  class DummyStrategyDomainStats : public autopas::TuningStrategyInterface {
   public:
    bool needsLiveInfo() const override { return false; }
    bool needsDomainSimilarityStatistics() const override { return true; }
    autopas::TuningStrategyOption getOptionType() const override { return autopas::TuningStrategyOption::fullSearch; }
    bool optimizeSuggestions(std::vector<autopas::Configuration> &configQueue,
                             const autopas::EvidenceCollection &evidenceCollection) override {
      return true;
    }
    bool reset(size_t iteration, size_t tuningPhase, std::vector<autopas::Configuration> &configQueue,
               const autopas::EvidenceCollection &evidenceCollection) override {
      return true;
    }
  };

  // Loop over all 4 permutations:
  // [None, Only LiveInfoStrat, Only DomainStatsStrat, Both]
  for (const bool includeLiveInfoStrat : {false, true}) {
    for (const bool includeDomainStatsStrat : {false, true}) {
      autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

      if (includeLiveInfoStrat) {
        tuningStrategies.emplace_back(std::make_unique<DummyStrategyLiveInfo>());
      }
      if (includeDomainStatsStrat) {
        tuningStrategies.emplace_back(std::make_unique<DummyStrategyDomainStats>());
      }

      auto autoTuner = std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, 20, "");

      const auto tunerManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
      tunerManager->addAutoTuner(std::move(autoTuner), autopas::InteractionTypeOption::pairwise);

      // Iteration 0: Start of Phase -> Should need info ONLY if a strategy requests LiveInfo
      EXPECT_EQ(tunerManager->needsLiveInfo(0), includeLiveInfoStrat)
          << "Failed at Iteration 0 with includeLiveInfoStrat=" << includeLiveInfoStrat
          << " and includeDomainStatsStrat=" << includeDomainStatsStrat;

      // Iteration 50: Middle of Phase -> Should ALWAYS be false, regardless of strategies
      EXPECT_FALSE(tunerManager->needsLiveInfo(50))
          << "Failed at Iteration 50 with includeLiveInfoStrat=" << includeLiveInfoStrat
          << " and includeDomainStatsStrat=" << includeDomainStatsStrat;

      // Iteration 95: Phase about to begin -> Should need info ONLY if a strategy requests DomainStats
      EXPECT_EQ(tunerManager->needsLiveInfo(95), includeDomainStatsStrat)
          << "Failed at Iteration 95 with includeLiveInfoStrat=" << includeLiveInfoStrat
          << " and includeDomainStatsStrat=" << includeDomainStatsStrat;

      // Iteration 100: Start of new Phase -> Should need info ONLY if a strategy requests LiveInfo
      EXPECT_EQ(tunerManager->needsLiveInfo(100), includeLiveInfoStrat)
          << "Failed at Iteration 100 with includeLiveInfoStrat=" << includeLiveInfoStrat
          << " and includeDomainStatsStrat=" << includeDomainStatsStrat;
    }
  }
}