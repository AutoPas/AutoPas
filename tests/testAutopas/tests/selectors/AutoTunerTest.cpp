/**
 * @file AutoTunerTest.cpp
 * @author F. Gratl
 * @date 8/10/18
 */

#include "AutoTunerTest.h"
#include "autopas/selectors/AutoTuner.h"

TEST_F(AutoTunerTest, testAllConfigurations) {
  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 42};
  // adaptive domain size so sliced is always applicable.
  bBoxMax[2] = autopas::autopas_get_max_threads() * 2;
  const double cutoff = 1;
  const double verletSkin = 0;
  const unsigned int verletRebuildFrequency = 1;
  const unsigned int maxSamples = 2;

  autopas::LJFunctor<Particle, FPCell> functor(cutoff, 1., 1., 0.);
  autopas::AutoTuner<Particle, FPCell> autoTuner(bBoxMin, bBoxMax, cutoff, verletSkin, verletRebuildFrequency,
                                                 autopas::allContainerOptions, autopas::allTraversalOptions,
                                                 autopas::allDataLayoutOptions, autopas::allNewton3Options,
                                                 autopas::SelectorStrategy::fastestAbs, 100, maxSamples);

  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::off);
  bool stillTuning = true;
  auto prevConfig = autopas::Configuration(autopas::ContainerOption(-1), autopas::TraversalOption(-1),
                                           autopas::DataLayoutOption(-1), autopas::Newton3Option(-1));

  int collectedSamples = 0;
  while (stillTuning) {
    if (collectedSamples == maxSamples) {
      collectedSamples = 0;
    }

    stillTuning = autoTuner.iteratePairwise(&functor);
    ++collectedSamples;

    auto currentConfig = autoTuner.getCurrentConfig();
    if (collectedSamples == 1) {
      EXPECT_NE(currentConfig, prevConfig);
    } else {
      EXPECT_EQ(currentConfig, prevConfig);
    }
    prevConfig = currentConfig;
  }
}

// @TODO:
// TEST_F(AutoTunerTest, testSelectOptimalTraversalFastestAbs) {
//  auto strategy = autopas::SelectorStrategy::fastestAbs;
//
//  mapConfigTime measurements;
//
//  measurements[] = {22, 14};
//  measurements[] = {30, 10};
//
//  mapConfigTime ignoredMeasurements;
//  ignoredMeasurements[autopas::TraversalOption::c08] = {1};
//
//  testFastest(strategy, measurements, autopas::TraversalOption::sliced, ignoredMeasurements);
//}

// TEST_F(TraversalSelectorTest, testSelectOptimalTraversalFastestMean) {
//  auto strategy = autopas::SelectorStrategy::fastestMean;
//
//  mapOptionsTime measurements;
//
//  measurements[autopas::TraversalOption::c08] = {2, 20};
//  measurements[autopas::TraversalOption::sliced] = {5, 7};
//
//  testFastest(strategy, measurements, autopas::TraversalOption::sliced);
//}
//
// TEST_F(TraversalSelectorTest, testSelectOptimalTraversalFastestMedian) {
//  auto strategy = autopas::SelectorStrategy::fastestMedian;
//
//  mapOptionsTime measurements;
//
//  measurements[autopas::TraversalOption::c08] = {4, 1, 5};
//  measurements[autopas::TraversalOption::sliced] = {2, 3, 3, 100};
//
//  testFastest(strategy, measurements, autopas::TraversalOption::sliced);
//}

//void AutoTunerTest::testFastest(autopas::SelectorStrategy strategy, mapConfigTime configAndTimes,
//                                autopas::Configuration expectedBest, mapConfigTime ignoredConfigAndTimes) {
//  MFunctor functor;
//
//  unsigned int maxSamples = 4;
//
//  std::set<autopas::Configuration> uniqueConfigs;
//
//  std::set<autopas::ContainerOption> containerOptions;
//  std::set<autopas::TraversalOption> traversalOptions;
//  std::set<autopas::DataLayoutOption> dataLayoutOptions;
//  std::set<autopas::Newton3Option> newton3Options;
//
//  for (auto &&m : configAndTimes) {
//    uniqueConfigs.insert(m.first);
//  }
//  for (auto &&m : ignoredConfigAndTimes) {
//    uniqueConfigs.insert(m.first);
//  }
//  for (auto &&c : uniqueConfigs) {
//    containerOptions.insert(c._container);
//    traversalOptions.insert(c._traversal);
//    dataLayoutOptions.insert(c._dataLayout);
//    newton3Options.insert(c._newton3);
//  }
//
//  constexpr size_t domainSize = 1000;
//
//  autopas::AutoTuner<Particle, FPCell> tuner(
//      {0, 0, 0}, {domainSize, domainSize, domainSize}, 1, 0, 100,
//      std::vector<autopas::ContainerOption>(containerOptions.begin(), containerOptions.end()),
//      std::vector<autopas::TraversalOption>(traversalOptions.begin(), traversalOptions.end()),
//      std::vector<autopas::DataLayoutOption>(dataLayoutOptions.begin(), dataLayoutOptions.end()),
//      std::vector<autopas::Newton3Option>(newton3Options.begin(), newton3Options.end()), strategy, 1000, maxSamples);
//
//  // check that all expected configurations are created
//  ASSERT_EQ(uniqueConfigs.size(), tuner.getAllowedConfigurations().size());

  //  autopas::TraversalSelector<FPCell> traversalSelector({domainSize, domainSize, domainSize}, optionVector);
  //
  //  EXPECT_THROW((traversalSelector.selectOptimalTraversal<MFunctor, true, true>(strategy, functor)), std::exception);
  //
  //  for (auto &&m : configAndTimes) {
  //    for (auto &&t : m.second) {
  //      EXPECT_CALL(functor, isRelevantForTuning()).WillOnce(Return(true));
  //      traversalSelector.addTimeMeasurement(functor, m.first, t);
  //    }
  //  }
  //
  //  for (auto &&m : ignoredConfigAndTimes) {
  //    for (auto &&t : m.second) {
  //      EXPECT_CALL(functor, isRelevantForTuning()).WillOnce(Return(false));
  //      traversalSelector.addTimeMeasurement(functor, m.first, t);
  //    }
  //  }
  //
  //  auto traversal = traversalSelector.selectOptimalTraversal<MFunctor, true, true>(strategy, functor);
  //  EXPECT_EQ(expectedBest, traversal->getTraversalType());
  //
  //  // select optimal traversal should delete all configAndTimes
  //  EXPECT_THROW((traversalSelector.selectOptimalTraversal<MFunctor, true, true>(strategy, functor)), std::exception);
//}