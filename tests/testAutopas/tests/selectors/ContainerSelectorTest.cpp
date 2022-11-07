/**
 * @file ContainerSelectorTest.cpp
 * @author F. Gratl
 * @date 22.06.18
 */

#include "ContainerSelectorTest.h"

TEST_F(ContainerSelectorTest, testSelectAndGetCurrentContainer) {
  autopas::ContainerSelector<Particle> containerSelector(bBoxMin, bBoxMax, cutoff);
  autopas::ContainerSelectorInfo containerInfo(cellSizeFactor, verletSkinPerTimestep, verletRebuildFrequency, 64,
                                               autopas::LoadEstimatorOption::none);

  // expect an exception if nothing is selected yet
  EXPECT_THROW((containerSelector.getCurrentContainer()), autopas::utils::ExceptionHandler::AutoPasException);

  // test all individual options
  for (auto containerOp : autopas::ContainerOption::getAllOptions()) {
    containerSelector.selectContainer(containerOp, containerInfo);

    EXPECT_EQ(containerOp, containerSelector.getCurrentContainer()->getContainerType());
  }
}