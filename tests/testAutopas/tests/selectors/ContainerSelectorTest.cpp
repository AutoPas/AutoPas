/**
 * @file ContainerSelectorTest.cpp
 * @author F. Gratl
 * @date 22.06.18
 */

#include "ContainerSelectorTest.h"

namespace ContainerSelectorTest {

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::ValuesIn;

// must be TEST_F because the logger which is called in the LC constructor is part of the fixture
TEST_F(ContainerSelectorTest, testSelectAndGetCurrentContainer) {
  autopas::ContainerSelector<Particle> containerSelector(bBoxMin, bBoxMax, cutoff);
  autopas::ContainerSelectorInfo containerInfo(cellSizeFactor, verletSkin, 64, autopas::LoadEstimatorOption::none);

  // expect an exception if nothing is selected yet
  EXPECT_THROW((containerSelector.getCurrentContainer()), autopas::utils::ExceptionHandler::AutoPasException);

  // test all individual options
  for (auto containerOp : autopas::ContainerOption::getAllOptions()) {
    containerSelector.selectContainer(containerOp, containerInfo);

    EXPECT_EQ(containerOp, containerSelector.getCurrentContainer()->getContainerType());
  }
}

INSTANTIATE_TEST_SUITE_P(Generated, ContainerSelectorTest,
                         ::testing::ValuesIn(autopas::ContainerOption::getAllOptions()),
                         ContainerSelectorTest::oneParamToString());

}  // end namespace ContainerSelectorTest
