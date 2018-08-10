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

  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double verletSkin = 0;
  const unsigned int verletRebuildFrequency = 1;
  autopas::AutoTuner<Particle, FPCell> autoTuner(bBoxMin, bBoxMax, cutoff, verletSkin, verletRebuildFrequency,
                                                 containers, traversals, 100);

  bool stillTuning = true;
  int i = 0;
  for (; stillTuning; ++i) {
    stillTuning = autoTuner.iteratePairwise(&functor, autopas::DataLayoutOption::aos);

    auto container = autoTuner.getContainer();

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
      case 3: {
        EXPECT_TRUE((dynamic_cast<autopas::LinkedCells<Particle, FPCell>*>(container.get())));
        break;
      }
      case 4: {
        // direct sum should be the fastest since it has the least overhead
        EXPECT_TRUE((dynamic_cast<autopas::DirectSum<Particle, FPCell>*>(container.get())))
            << "tune() selected the wrong container after collecting all timings";
        EXPECT_FALSE(stillTuning) << "tune() returns true(=still tuning) after checking all options!";
        break;
      }
      default:
        FAIL() << "Tuning took more iterations than expected!";
    }
  }

  EXPECT_EQ(i, 5) << "Too unexpected number of tuning iterations!";

  auto container = autoTuner.getContainer();
  EXPECT_TRUE((dynamic_cast<autopas::DirectSum<Particle, FPCell>*>(container.get())))
            << "tune() returned the wrong container after tuning phase";
}