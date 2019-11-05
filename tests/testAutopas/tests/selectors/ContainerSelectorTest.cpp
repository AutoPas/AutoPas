/**
 * @file ContainerSelectorTest.cpp
 * @author F. Gratl
 * @date 22.06.18
 */

#include "ContainerSelectorTest.h"

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::ValuesIn;

// must be TEST_F because the logger which is called in the LC constructor is part of the fixture
TEST_F(ContainerSelectorTest, testSelectAndGetCurrentContainer) {
  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double cellSizeFactor = 1;
  const double verletSkin = 0;

  autopas::ContainerSelector<Particle, FPCell> containerSelector(bBoxMin, bBoxMax, cutoff);
  autopas::ContainerSelectorInfo containerInfo(cellSizeFactor, verletSkin);

  // expect an exception if nothing is selected yet
  EXPECT_THROW((containerSelector.getCurrentContainer()), autopas::utils::ExceptionHandler::AutoPasException);

  // test all individual options
  for (auto containerOp : autopas::ContainerOption::getAllOptions()) {
    containerSelector.selectContainer(containerOp, containerInfo);

    EXPECT_EQ(containerOp, containerSelector.getCurrentContainer()->getContainerType());
  }
}

/**
 * This function stores a copy of each particle depending on the position in ListInner, ListHalo or ListHaloVerletOnly.
 * @param bBoxMin Bounding box min.
 * @param bBoxMax Bounding box max.
 * @param cutoff Cutoff radius.
 * @param containerSelector Container selector used to retrieve the current container.
 * @param ListInner All particles inside the bounding box.
 * @param ListHalo All particles in the halo.
 * @param ListHaloVerletOnly All particles in the halo.
 */
void getStatus(const std::array<double, 3> &bBoxMin, const std::array<double, 3> &bBoxMax, const double cutoff,
               autopas::ContainerSelector<Particle, FPCell> &containerSelector, std::vector<Particle> &ListInner,
               std::vector<Particle> &ListHalo, std::vector<Particle> &ListHaloVerletOnly) {
  for (auto iter = containerSelector.getCurrentContainer()->begin(autopas::IteratorBehavior::ownedOnly); iter.isValid();
       ++iter) {
    ListInner.push_back(*iter);
  }
  const auto haloBoxMin = autopas::ArrayMath::subScalar(bBoxMin, cutoff);
  const auto haloBoxMax = autopas::ArrayMath::addScalar(bBoxMax, cutoff);
  for (auto iter = containerSelector.getCurrentContainer()->begin(autopas::IteratorBehavior::haloOnly); iter.isValid();
       ++iter) {
    if (autopas::utils::inBox(iter->getR(), haloBoxMin, haloBoxMax)) {
      ListHalo.push_back(*iter);
    } else {
      ListHaloVerletOnly.push_back(*iter);
    }
  }
}

TEST_P(ContainerSelectorTest, testContainerConversion) {
  auto from = std::get<0>(GetParam());
  auto to = std::get<1>(GetParam());

  const std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double cellSizeFactor = 1;
  const double verletSkin = 0.1;

  autopas::ContainerSelector<Particle, FPCell> containerSelector(bBoxMin, bBoxMax, cutoff);
  autopas::ContainerSelectorInfo containerInfo(cellSizeFactor, verletSkin);

  // select container from which we want to convert from
  containerSelector.selectContainer(from, containerInfo);

  // fill witch problematic particles
  {
    auto container = containerSelector.getCurrentContainer();
    auto getPossible1DPositions = [&](double min, double max) -> auto {
      return std::array<double, 6>{min - cutoff - verletSkin,       min - cutoff, min, max, max + cutoff - 1e-3,
                                   min - cutoff - verletSkin - 1e-3};
    };
    size_t id = 0;
    // const auto haloBoxMin = autopas::ArrayMath::subScalar(bBoxMin, cutoff);
    // const auto haloBoxMax = autopas::ArrayMath::addScalar(bBoxMax, cutoff);

    for (auto x : getPossible1DPositions(bBoxMin[0], bBoxMax[0])) {
      for (auto y : getPossible1DPositions(bBoxMin[1], bBoxMax[1])) {
        for (auto z : getPossible1DPositions(bBoxMin[2], bBoxMax[2])) {
          const std::array<double, 3> pos{x, y, z};
          Particle p(pos, {0., 0., 0.}, id);
          if (autopas::utils::inBox(pos, bBoxMin, bBoxMax)) {
            container->addParticle(p);
          } else {
            container->addHaloParticle(p);
          }
          ++id;
        }
      }
    }
  }

  std::vector<Particle> beforeListInner, beforeListHalo,
      beforeListHaloVerletOnly /*for particles only in verlet containers*/;

  getStatus(bBoxMin, bBoxMax, cutoff, containerSelector, beforeListInner, beforeListHalo, beforeListHaloVerletOnly);

  // select container to which we want to convert to
  containerSelector.selectContainer(to, containerInfo);

  std::vector<Particle> afterListInner, afterListHalo, afterListHaloVerletOnly;

  getStatus(bBoxMin, bBoxMax, cutoff, containerSelector, afterListInner, afterListHalo, afterListHaloVerletOnly);

  EXPECT_THAT(afterListInner, UnorderedElementsAreArray(beforeListInner));
  EXPECT_THAT(afterListHalo, UnorderedElementsAreArray(beforeListHalo));
  if (to.to_string().find("Verlet") != std::string::npos and from.to_string().find("Verlet") != std::string::npos) {
    EXPECT_THAT(afterListHaloVerletOnly, UnorderedElementsAreArray(beforeListHaloVerletOnly));
  }
}

/// @todo: use this instead of below to enable testing of VerletClusterLists.
// INSTANTIATE_TEST_SUITE_P(Generated, ContainerSelectorTest,
//                         Combine(ValuesIn(autopas::ContainerOption::getAllOptions()),
//                         ValuesIn(autopas::ContainerOption::getAllOptions())),
//                         ContainerSelectorTest::PrintToStringParamName());

INSTANTIATE_TEST_SUITE_P(Generated, ContainerSelectorTest,
                         Combine(ValuesIn([]() -> std::set<autopas::ContainerOption> {
                                   auto all = autopas::ContainerOption::getAllOptions();
                                   all.erase(all.find(autopas::ContainerOption::verletClusterLists));
                                   return all;
                                 }()),
                                 ValuesIn([]() -> std::set<autopas::ContainerOption> {
                                   auto all = autopas::ContainerOption::getAllOptions();
                                   all.erase(all.find(autopas::ContainerOption::verletClusterLists));
                                   return all;
                                 }())),
                         ContainerSelectorTest::PrintToStringParamName());
