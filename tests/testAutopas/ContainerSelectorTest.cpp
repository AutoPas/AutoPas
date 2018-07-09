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