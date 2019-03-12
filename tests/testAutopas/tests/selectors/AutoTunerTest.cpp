/**
 * @file AutoTunerTest.cpp
 * @author F. Gratl
 * @date 8/10/18
 */

#include "AutoTunerTest.h"
#include "autopas/selectors/AutoTuner.h"

using ::testing::_;

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

  // number of possible configurations * number of samples
  size_t expectedNumberOfIterations = 40 * maxSamples;

  int collectedSamples = 0;
  int iterations = 0;
  while (stillTuning) {
    if (collectedSamples == maxSamples) {
      collectedSamples = 0;
    }

    stillTuning = autoTuner.iteratePairwise(&functor);
    ++iterations;
    ++collectedSamples;

    auto currentConfig = autoTuner.getCurrentConfig();
    if (collectedSamples == 1) {
      EXPECT_NE(currentConfig, prevConfig);
    } else {
      EXPECT_EQ(currentConfig, prevConfig);
    }
    prevConfig = currentConfig;
  }

  EXPECT_EQ(expectedNumberOfIterations, iterations);
}

TEST_F(AutoTunerTest, testSelectOptimalTraversalFastestAbs) {
  mapConfigTime configTimes;

  autopas::Configuration bestConfig(autopas::ContainerOption::directSum, autopas::TraversalOption::directSumTraversal,
                                    autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  configTimes[bestConfig] = {30, 10};
  configTimes[autopas::Configuration(autopas::ContainerOption::directSum, autopas::TraversalOption::directSumTraversal,
                                     autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled)] = {22, 14};

  mapConfigTime ignoredConfigTimes;
  ignoredConfigTimes[autopas::Configuration(autopas::ContainerOption::directSum,
                                            autopas::TraversalOption::directSumTraversal,
                                            autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled)] = {1};

  testFastest(autopas::SelectorStrategy::fastestAbs, configTimes, bestConfig, ignoredConfigTimes);
}

TEST_F(AutoTunerTest, testSelectOptimalTraversalFastestMean) {
  mapConfigTime configTimes;

  autopas::Configuration bestConfig(autopas::ContainerOption::directSum, autopas::TraversalOption::directSumTraversal,
                                    autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  configTimes[bestConfig] = {5, 7};
  configTimes[autopas::Configuration(autopas::ContainerOption::directSum, autopas::TraversalOption::directSumTraversal,
                                     autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled)] = {2, 20};

  mapConfigTime ignoredConfigTimes;
  ignoredConfigTimes[autopas::Configuration(autopas::ContainerOption::directSum,
                                            autopas::TraversalOption::directSumTraversal,
                                            autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled)] = {1, 5};

  testFastest(autopas::SelectorStrategy::fastestMean, configTimes, bestConfig, ignoredConfigTimes);
}

TEST_F(AutoTunerTest, testSelectOptimalTraversalFastestMedian) {
  mapConfigTime configTimes;

  autopas::Configuration bestConfig(autopas::ContainerOption::directSum, autopas::TraversalOption::directSumTraversal,
                                    autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled);

  configTimes[bestConfig] = {2, 3, 3, 100};
  configTimes[autopas::Configuration(autopas::ContainerOption::directSum, autopas::TraversalOption::directSumTraversal,
                                     autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled)] = {4, 1, 5, 2};

  mapConfigTime ignoredConfigTimes;
  ignoredConfigTimes[autopas::Configuration(autopas::ContainerOption::directSum,
                                            autopas::TraversalOption::directSumTraversal,
                                            autopas::DataLayoutOption::soa, autopas::Newton3Option::disabled)] = {1, 5};

  testFastest(autopas::SelectorStrategy::fastestMedian, configTimes, bestConfig, ignoredConfigTimes);
}

void AutoTunerTest::testFastest(autopas::SelectorStrategy strategy, mapConfigTime configAndTimes,
                                autopas::Configuration expectedBest, mapConfigTime ignoredConfigAndTimes) {
  unsigned int maxSamples = 4;

  std::set<autopas::Configuration> uniqueConfigs;

  std::set<autopas::ContainerOption> containerOptions;
  std::set<autopas::TraversalOption> traversalOptions;
  std::set<autopas::DataLayoutOption> dataLayoutOptions;
  std::set<autopas::Newton3Option> newton3Options;

  for (auto &&m : configAndTimes) {
    uniqueConfigs.insert(m.first);
  }
  for (auto &&m : ignoredConfigAndTimes) {
    uniqueConfigs.insert(m.first);
  }
  for (auto &&c : uniqueConfigs) {
    containerOptions.insert(c._container);
    traversalOptions.insert(c._traversal);
    dataLayoutOptions.insert(c._dataLayout);
    newton3Options.insert(c._newton3);
  }

  constexpr size_t domainSize = 1000;

  autopas::AutoTuner<Particle, FPCell> tuner(
      {0, 0, 0}, {domainSize, domainSize, domainSize}, 1, 0, 100,
      std::vector<autopas::ContainerOption>(containerOptions.begin(), containerOptions.end()),
      std::vector<autopas::TraversalOption>(traversalOptions.begin(), traversalOptions.end()),
      std::vector<autopas::DataLayoutOption>(dataLayoutOptions.begin(), dataLayoutOptions.end()),
      std::vector<autopas::Newton3Option>(newton3Options.begin(), newton3Options.end()), strategy, 1000, maxSamples);

  // check that all expected configurations are created
  ASSERT_EQ(uniqueConfigs.size(), tuner.getAllowedConfigurations().size());

  MFunctor functor;

  for (auto &&m : configAndTimes) {
    for (auto &&t : m.second) {
      EXPECT_CALL(functor, isRelevantForTuning()).WillOnce(::testing::Return(true));
      tuner.addTimeMeasurement(functor, t);
    }
  }

  for (auto &&m : ignoredConfigAndTimes) {
    for (auto &&t : m.second) {
      EXPECT_CALL(functor, isRelevantForTuning()).WillOnce(::testing::Return(false));
      tuner.addTimeMeasurement(functor, t);
    }
  }

  // some necessary calls...
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(::testing::Return(true));
  EXPECT_CALL(functor, SoALoader(_, _)).Times(::testing::AtLeast(1));
  EXPECT_CALL(functor, SoAExtractor(_, _)).Times(::testing::AtLeast(1));
  // Only needed if verlet Lists are present
  //  EXPECT_CALL(functor, SoALoader(_,_,_)).Times(::testing::AtLeast(1));
  //  EXPECT_CALL(functor, SoAExtractor(_,_,_)).Times(::testing::AtLeast(1));

  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(::testing::Return(true));
  // rush through the tuning procedure
  for (size_t i = 0; i < tuner.getAllowedConfigurations().size() * maxSamples; ++i) {
    tuner.iteratePairwise(&functor);
  }

  EXPECT_EQ(expectedBest, tuner.getCurrentConfig());
}