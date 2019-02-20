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

TEST_F(AutoTunerTest, testTuneAoS) { testTune(autopas::DataLayoutOption::aos); }

TEST_F(AutoTunerTest, testTuneSoA) { testTune(autopas::DataLayoutOption::soa); }