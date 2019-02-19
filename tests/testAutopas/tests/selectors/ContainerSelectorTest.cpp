/**
 * @file ContainerSelectorTest.cpp
 * @author F. Gratl
 * @date 22.06.18
 */

#include "ContainerSelectorTest.h"

// must be TEST_F because the logger which is called in the LC constructor is part of the fixture
TEST_F(ContainerSelectorTest, testGetOptimalContainerOneOption) {
  std::vector<autopas::ContainerOption> optionVectorDir = {autopas::ContainerOption::directSum};
  std::vector<autopas::ContainerOption> optionVectorLC = {autopas::ContainerOption::linkedCells};
  std::vector<autopas::ContainerOption> optionVectorVerlet = {autopas::ContainerOption::verletLists};

  std::vector<autopas::TraversalOption> traversals = {autopas::TraversalOption::c08};

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

  auto containerDir = containerSelectorDir.getCurrentContainer();
  auto containerLC = containerSelectorLC.getCurrentContainer();
  auto containerVerlet = containerSelectorVerlet.getCurrentContainer();

  EXPECT_TRUE((dynamic_cast<autopas::DirectSum<Particle, FPCell>*>(containerDir.get())));
  EXPECT_TRUE((dynamic_cast<autopas::LinkedCells<Particle, FPCell>*>(containerLC.get())));
  EXPECT_TRUE((dynamic_cast<autopas::VerletLists<Particle>*>(containerVerlet.get())));
}

TEST_F(ContainerSelectorTest, testNextContainer) {
  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double verletSkin = 0;
  const unsigned int verletRebuildFrequency = 1;

  std::vector<autopas::ContainerOption> containerOptions = {autopas::ContainerOption::verletLists,
                                                             autopas::ContainerOption::directSum,
                                                             autopas::ContainerOption::linkedCells};
  std::vector<autopas::TraversalOption> traversalOptions = {autopas::TraversalOption::c08};
  autopas::ContainerSelector<Particle, FPCell> containerSelector(
      bBoxMin, bBoxMax, cutoff, verletSkin, verletRebuildFrequency, containerOptions, traversalOptions);

  auto container = containerSelector.selectNextContainer();
  EXPECT_EQ(autopas::ContainerOption::verletLists, container->getContainerType());

  container = containerSelector.selectNextContainer();
  EXPECT_EQ(autopas::ContainerOption::directSum, container->getContainerType());

  container = containerSelector.selectNextContainer();
  EXPECT_EQ(autopas::ContainerOption::linkedCells, container->getContainerType());

  container = containerSelector.selectNextContainer();
  EXPECT_EQ(nullptr, container);

  container = containerSelector.selectNextContainer();
  EXPECT_EQ(nullptr, container);
}

TEST_F(ContainerSelectorTest, testSelectOptimalContainer) {
  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double verletSkin = 0;
  const unsigned int verletRebuildFrequency = 1;

  std::vector<autopas::ContainerOption> containerOptions = {autopas::ContainerOption::verletLists,
                                                             autopas::ContainerOption::directSum,
                                                             autopas::ContainerOption::linkedCells};
  std::vector<autopas::TraversalOption> traversalOptions = {autopas::TraversalOption::c08};
  autopas::ContainerSelector<Particle, FPCell> containerSelector(
      bBoxMin, bBoxMax, cutoff, verletSkin, verletRebuildFrequency, containerOptions, traversalOptions);

  EXPECT_THROW(containerSelector.selectOptimalContainer(), std::exception);

  containerSelector.addTimeMeasurement(containerOptions[2], 20);
  containerSelector.addTimeMeasurement(containerOptions[1], 10);
  containerSelector.addTimeMeasurement(containerOptions[0], 30);
  containerSelector.addTimeMeasurement(containerOptions[1], 22);

  auto container = containerSelector.selectOptimalContainer();
  EXPECT_EQ(containerOptions[1], container->getContainerType());

  // select optimal container should delete all measurements
  EXPECT_THROW(containerSelector.selectOptimalContainer(), std::exception);
}