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
  const unsigned int verletRebuildFrequency = 1;

  autopas::ContainerSelector<Particle, FPCell> containerSelector(bBoxMin, bBoxMax, cutoff, cellSizeFactor, verletSkin,
                                                                 verletRebuildFrequency);

  // expect an exception if nothing is selected yet
  EXPECT_THROW((containerSelector.getCurrentContainer()), autopas::utils::ExceptionHandler::AutoPasException);

  // test all individual options
  for (auto containerOp : autopas::allContainerOptions) {
    containerSelector.selectContainer(containerOp);

    EXPECT_EQ(containerOp, containerSelector.getCurrentContainer()->getContainerType());
  }
}

TEST_P(ContainerSelectorTest, testContainerConversion) {
  auto from = std::get<0>(GetParam());
  auto to = std::get<1>(GetParam());

  std::array<double, 3> bBoxMin = {0, 0, 0}, bBoxMax = {10, 10, 10};
  const double cutoff = 1;
  const double cellSizeFactor = 1;
  const double verletSkin = 0.1;
  const unsigned int verletRebuildFrequency = 1;

  autopas::ContainerSelector<Particle, FPCell> containerSelector(bBoxMin, bBoxMax, cutoff, cellSizeFactor, verletSkin,
                                                                 verletRebuildFrequency);
  // select container from which we want to convert from
  containerSelector.selectContainer(from);

  // fill witch problematic particles
  {
    auto container = containerSelector.getCurrentContainer();
    auto getPossible1D = [&](double min, double max) -> auto {
      return std::initializer_list<double>{min - cutoff - verletSkin,       min - cutoff, min, max, max + cutoff - 1e-3,
                                           min - cutoff - verletSkin - 1e-3};
    };
    size_t id = 0;
    for (auto x : getPossible1D(bBoxMin[0], bBoxMax[0])) {
      for (auto y : getPossible1D(bBoxMin[1], bBoxMax[1])) {
        for (auto z : getPossible1D(bBoxMin[2], bBoxMax[2])) {
          std::array<double, 3> pos{x, y, z};
          Particle p(pos, {0., 0., 0.}, id);
          if (autopas::utils::inBox(pos, bBoxMin, bBoxMax)) {
            container->addParticle(p);
          } else {
            if (autopas::utils::inBox(pos, autopas::ArrayMath::sub(bBoxMin, std::array<double, 3>{cutoff}),
                                      autopas::ArrayMath::add(bBoxMin, std::array<double, 3>{cutoff})) or
                autopas::utils::StringUtils::to_string(container->getContainerType()).find("Verlet") !=
                    std::string::npos) {
              container->addHaloParticle(p);
            }
          }
          ++id;
        }
      }
    }
  }
  std::vector<Particle> beforeListInner, beforeListHalo,
      beforeListHaloVerletOnly /*for particles only in verlet containers*/;
  for (auto iter = containerSelector.getCurrentContainer()->begin(autopas::IteratorBehavior::ownedOnly); iter.isValid();
       ++iter) {
    beforeListInner.push_back(*iter);
  }
  for (auto iter = containerSelector.getCurrentContainer()->begin(autopas::IteratorBehavior::haloOnly); iter.isValid();
       ++iter) {
    if (autopas::utils::inBox(iter->getR(), autopas::ArrayMath::sub(bBoxMin, std::array<double, 3>{cutoff}),
                              autopas::ArrayMath::add(bBoxMin, std::array<double, 3>{cutoff}))) {
      beforeListHalo.push_back(*iter);
    } else {
      beforeListHaloVerletOnly.push_back(*iter);
    }
  }

  // select container to which we want to convert to
  containerSelector.selectContainer(to);

  std::vector<Particle> afterListInner, afterListHalo, afterListHaloVerletOnly;
  for (auto iter = containerSelector.getCurrentContainer()->begin(autopas::IteratorBehavior::ownedOnly); iter.isValid();
       ++iter) {
    afterListInner.push_back(*iter);
  }
  for (auto iter = containerSelector.getCurrentContainer()->begin(autopas::IteratorBehavior::haloOnly); iter.isValid();
       ++iter) {
    if (autopas::utils::inBox(iter->getR(), autopas::ArrayMath::sub(bBoxMin, std::array<double, 3>{cutoff}),
                              autopas::ArrayMath::add(bBoxMin, std::array<double, 3>{cutoff}))) {
      afterListHalo.push_back(*iter);
    } else {
      afterListHaloVerletOnly.push_back(*iter);
    }
  }

  EXPECT_THAT(afterListInner, UnorderedElementsAreArray(beforeListInner));
  EXPECT_THAT(afterListHalo, UnorderedElementsAreArray(beforeListHalo));
  if (autopas::utils::StringUtils::to_string(to).find("Verlet") != std::string::npos and
      autopas::utils::StringUtils::to_string(from).find("Verlet") != std::string::npos) {
    EXPECT_THAT(afterListHaloVerletOnly, UnorderedElementsAreArray(beforeListHaloVerletOnly));
  }
}

INSTANTIATE_TEST_SUITE_P(
    Generated, ContainerSelectorTest,
    Combine(ValuesIn([]() -> std::vector<autopas::ContainerOption> {
              auto all = autopas::allContainerOptions;
              all.erase(std::remove(all.begin(), all.end(), autopas::ContainerOption::verletClusterLists), all.end());
              return all;
            }()),
            ValuesIn([]() -> std::vector<autopas::ContainerOption> {
              auto all = autopas::allContainerOptions;
              all.erase(std::remove(all.begin(), all.end(), autopas::ContainerOption::verletClusterLists), all.end());
              return all;
            }())),
    ContainerSelectorTest::PrintToStringParamName());
