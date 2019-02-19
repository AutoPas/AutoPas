/**
 * @file ContainerSelectorTest.cpp
 * @author F. Gratl
 * @date 22.06.18
 */

#include "ContainerSelectorTest.h"

// must be TEST_F because the logger which is called in the LC constructor is part of the fixture
TEST_F(ContainerSelectorTest, testSelectAndGetCurrentContainer) {
  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double verletSkin = 0;
  const unsigned int verletRebuildFrequency = 1;

  autopas::ContainerSelector<Particle, FPCell> containerSelector(bBoxMin, bBoxMax, cutoff, verletSkin,
                                                                 verletRebuildFrequency, autopas::allContainerOptions,
                                                                 autopas::allTraversalOptions);
  for (auto containerOp : autopas::allContainerOptions) {
    containerSelector.selectContainer(containerOp);

    EXPECT_EQ(containerOp, containerSelector.getCurrentContainer()->getContainerType());
  }
}
