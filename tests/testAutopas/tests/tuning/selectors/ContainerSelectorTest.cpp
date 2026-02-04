/**
 * @file ContainerSelectorTest.cpp
 * @author F. Gratl
 * @date 22.06.18
 */

#include "ContainerSelectorTest.h"

#include "autopas/tuning/selectors/ContainerSelector.h"
#include "testingHelpers/commonTypedefs.h"

TEST_F(ContainerSelectorTest, testSelectAndGetCurrentContainer) {
  autopas::ContainerSelectorInfo containerInfo(bBoxMin, bBoxMax, cutoff, cellSizeFactor, verletSkin, 64, 8,
  autopas::LoadEstimatorOption::none, orderCellsByMortonIndex, preloadLJMixingPtr, useSoAIndex, reserveVLSizes,
  bucketSortParticles, sortVerletLists, sortingFrequency);

  // test all individual options
  for (auto containerOp : autopas::ContainerOption::getAllOptions()) {
    const auto container = autopas::ContainerSelector<ParticleFP64>::generateContainer(containerOp, containerInfo);

    EXPECT_EQ(containerOp, container->getContainerType());
  }
}