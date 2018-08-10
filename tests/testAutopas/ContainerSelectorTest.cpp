/**
 * @file ContainerSelectorTest.cpp
 * @author F. Gratl
 * @date 22.06.18
 */

#include "ContainerSelectorTest.h"

// must be TEST_F because the logger which is called in the LC constructor is part of the fixture
TEST_F(ContainerSelectorTest, testGetOptimalContainerOneOption) {
  std::vector<autopas::ContainerOptions> optionVectorDir = {autopas::ContainerOptions::directSum};
  std::vector<autopas::ContainerOptions> optionVectorLC = {autopas::ContainerOptions::linkedCells};
  std::vector<autopas::ContainerOptions> optionVectorVerlet = {autopas::ContainerOptions::verletLists};

  std::vector<autopas::TraversalOptions> traversals = {autopas::TraversalOptions::c08};

  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double verletSkin = 0;
  const unsigned int verletRebuildFrequency = 1;
  autopas::ContainerSelector<Particle, FPCell> containerSelectorDir(
      bBoxMin, bBoxMax, cutoff, verletSkin, verletRebuildFrequency, optionVectorDir, traversals);
  autopas::ContainerSelector<Particle, FPCell> containerSelectorLC(bBoxMin, bBoxMax, cutoff, verletSkin,
                                                                   verletRebuildFrequency, optionVectorLC, traversals);
  autopas::ContainerSelector<Particle, FPCell> containerSelectorVerlet(
      bBoxMin, bBoxMax, cutoff, verletSkin, verletRebuildFrequency, optionVectorVerlet, traversals);

  auto containerDir = containerSelectorDir.getOptimalContainer();
  auto containerLC = containerSelectorLC.getOptimalContainer();
  auto containerVerlet = containerSelectorVerlet.getOptimalContainer();

  EXPECT_TRUE((dynamic_cast<autopas::DirectSum<Particle, FPCell>*>(containerDir.get())));
  EXPECT_TRUE((dynamic_cast<autopas::LinkedCells<Particle, FPCell>*>(containerLC.get())));
  EXPECT_TRUE((dynamic_cast<autopas::VerletLists<Particle>*>(containerVerlet.get())));
}

TEST_F(ContainerSelectorTest, testTune) {
  std::vector<autopas::ContainerOptions> containers = {autopas::ContainerOptions::directSum,
                                                       autopas::ContainerOptions::verletLists,
                                                       autopas::ContainerOptions::linkedCells};
  std::vector<autopas::TraversalOptions> traversals = {autopas::TraversalOptions::c08,
                                                       autopas::TraversalOptions::sliced};

  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double verletSkin = 0;
  const unsigned int verletRebuildFrequency = 1;
  autopas::ContainerSelector<Particle, FPCell> containerSelector(bBoxMin, bBoxMax, cutoff, verletSkin,
                                                                 verletRebuildFrequency, containers, traversals);

  bool stillTuning = true;
  int i = 0;
  for (; stillTuning; ++i) {
    stillTuning = containerSelector.tune();

    auto container = containerSelector.getOptimalContainer();

    switch (i) {
      case 0: {
        EXPECT_TRUE((dynamic_cast<autopas::DirectSum<Particle, FPCell>*>(container.get())));
        containerSelector.addTimeMeasurement(autopas::ContainerOptions::directSum, 20);
        break;
      }
      case 1: {
        EXPECT_TRUE((dynamic_cast<autopas::VerletLists<Particle>*>(container.get())));
        containerSelector.addTimeMeasurement(autopas::ContainerOptions::verletLists, 30);
        break;
      }
      case 2: {
        EXPECT_TRUE((dynamic_cast<autopas::LinkedCells<Particle, FPCell>*>(container.get())));
        containerSelector.addTimeMeasurement(autopas::ContainerOptions::linkedCells, 10);
        break;
      }
      case 3: {
        EXPECT_TRUE((dynamic_cast<autopas::LinkedCells<Particle, FPCell>*>(container.get())))
            << "tune() selected the wrong container after collecting all timings";
        EXPECT_FALSE(stillTuning) << "tune() returns true(=still tuning) after checking all options!";
        break;
      }
      default:
        FAIL() << "Tuning took more turns than expected!";
    }
  }

  EXPECT_EQ(i, 4) << "Too unexpected number of tuning iterations!";

  auto container = containerSelector.getOptimalContainer();
  EXPECT_TRUE((dynamic_cast<autopas::LinkedCells<Particle, FPCell>*>(container.get())))
      << "tune() returned the wrong container after tuning phase";
}