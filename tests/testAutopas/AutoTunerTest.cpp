/**
 * @file AutoTunerTest.cpp
 * @author F. Gratl
 * @date 8/10/18
 */

#include "AutoTunerTest.h"
#include <autopas/selectors/AutoTuner.h>

TEST_F(AutoTunerTest, testTune) {
  autopas::LJFunctor<Particle, FPCell> functor;
  std::vector<autopas::ContainerOptions> containers = {autopas::ContainerOptions::verletLists,
                                                       autopas::ContainerOptions::directSum,
                                                       autopas::ContainerOptions::linkedCells};
  std::vector<autopas::TraversalOptions> traversals = {autopas::TraversalOptions::sliced,
                                                       autopas::TraversalOptions::c08};

  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 100};
  const double cutoff = 1;
  const double verletSkin = 0;
  const unsigned int verletRebuildFrequency = 1;
  autopas::AutoTuner<Particle, FPCell> autoTuner(bBoxMin, bBoxMax, cutoff, verletSkin, verletRebuildFrequency,
                                                 containers, traversals, 100);

  std::shared_ptr<autopas::ParticleContainer<Particle, FPCell>> fastestContainer;
  autopas::Logger::get()->set_level(autopas::Logger::LogLevel::debug);
  bool stillTuning = true;
  int i = 0;
  for (; stillTuning; ++i) {
    // AutoPasLog(debug, "Iteration {}", i);
    stillTuning = autoTuner.iteratePairwise(&functor, autopas::DataLayoutOption::aos);

    auto container = autoTuner.getContainer();

    // tuning phases:
    // 0 -> test verlet
    // 1 -> test directSum
    // 2 -> test linked with sliced
    // 3 -> test linked with c08
    // 4 -> choose best lc traversal -> traversal tuning finished
    // 5 -> choose best container -> tuning finished
    // 6 -> normal iteration using optimal combination
    switch (i) {
      case 0: {
        EXPECT_TRUE((dynamic_cast<autopas::VerletLists<Particle>*>(container.get())));
        break;
      }
      case 1: {
        EXPECT_TRUE((dynamic_cast<autopas::DirectSum<Particle, FPCell>*>(container.get())));
        break;
      }
      // only here both traversals are checked
      case 2:
      case 3:
      case 4: {
        EXPECT_TRUE((dynamic_cast<autopas::LinkedCells<Particle, FPCell>*>(container.get())));
        break;
      }
      case 5: {
        // the fastest container might be nondeterministic here due to hardware constrains so just remember it
        // and check if the selector returns the same later
        fastestContainer = container;
        EXPECT_FALSE(stillTuning) << "tune() returns true(=still tuning) after checking all options!";
        break;
      }
      default:
        FAIL() << "Tuning took more iterations than expected!";
    }
  }

  EXPECT_EQ(i, 6) << "Unexpected number of tuning iterations!";

  auto container = autoTuner.getContainer();
  EXPECT_EQ(fastestContainer->getContainerType(), container->getContainerType())
      << "tune() returned the wrong container after tuning phase";
}