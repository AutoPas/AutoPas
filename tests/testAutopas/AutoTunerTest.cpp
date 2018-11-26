/**
 * @file AutoTunerTest.cpp
 * @author F. Gratl
 * @date 8/10/18
 */

#include "AutoTunerTest.h"
#include <autopas/selectors/AutoTuner.h>

void AutoTunerTest::testTune(autopas::DataLayoutOption dataLayoutOption) {
  autopas::LJFunctor<Particle, FPCell> functor(1., 1., 1., 0.);
  std::vector<autopas::ContainerOptions> containers = {autopas::ContainerOptions::verletLists,
                                                       autopas::ContainerOptions::linkedCells,
                                                       autopas::ContainerOptions::directSumContainer};
  std::vector<autopas::TraversalOptions> traversals = {
      autopas::TraversalOptions::sliced, autopas::TraversalOptions::c08, autopas::TraversalOptions::directSum};

  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 42};
  // adaptive domain size so sliced is always applicable.
  bBoxMax[2] = autopas::autopas_get_max_threads() * 2;
  const double cutoff = 1;
  const double verletSkin = 0;
  const unsigned int verletRebuildFrequency = 1;
  const unsigned int numSamples = 2;
  autopas::AutoTuner<Particle, FPCell> autoTuner(bBoxMin, bBoxMax, cutoff, verletSkin, verletRebuildFrequency,
                                                 containers, traversals, 100, numSamples);

  std::shared_ptr<autopas::ParticleContainer<Particle, FPCell>> fastestContainer;
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  bool stillTuning = true;
  int i = 0;
  for (; stillTuning; ++i) {
    std::cout << "ITERATION " << i << std::endl;
    stillTuning = autoTuner.iteratePairwise(&functor, dataLayoutOption);

    auto container = autoTuner.getContainer();

    // tuning phases:
    // 0 -> test verlet
    // 1 -> test linked with sliced
    // 2 -> test linked with c08
    // 3 -> test directSumContainer
    // 4 -> choose best combination -> tuning finished and  normal iteration using optimal combination
    switch (i / numSamples) {
      case 0: {
        ASSERT_EQ(containers[0], container->getContainerType());
        break;
      }
      // only here both traversals are checked
      case 1:
      case 2: {
        ASSERT_EQ(containers[1], container->getContainerType());
        break;
      }
      case 3: {
        ASSERT_EQ(containers[2], container->getContainerType());
        break;
      }
      case 4: {
        // the fastest container might be nondeterministic here due to hardware constrains so just remember it
        // and check if the selector returns the same later
        fastestContainer = container;
        ASSERT_FALSE(stillTuning) << "tune() returns true(=still tuning) after checking all options!";
        break;
      }
      default:
        FAIL() << "Tuning took more iterations than expected!";
    }
  }

  EXPECT_EQ(i, 4 * numSamples + 1) << "Unexpected number of tuning iterations!";

  auto container = autoTuner.getContainer();
  EXPECT_EQ(fastestContainer->getContainerType(), container->getContainerType())
      << "tune() returned the wrong container after tuning phase";
}

TEST_F(AutoTunerTest, testTuneAoS) { testTune(autopas::DataLayoutOption::aos); }

TEST_F(AutoTunerTest, testTuneSoA) { testTune(autopas::DataLayoutOption::soa); }