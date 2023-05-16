/**
 * @file ContainerSelectorTestFromTo.cpp
 * @author F. Gratl
 * @date 14.12.2020
 */

#include "ContainerSelectorTestFromTo.h"

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::ValuesIn;

/**
 * This function stores a copy of each particle depending on the position in ListInner, ListHaloWithinCutoff or
 * ListHaloOutsideCutoff.
 * @param bBoxMin Bounding box min.
 * @param bBoxMax Bounding box max.
 * @param cutoff Cutoff radius.
 * @param containerSelector Container selector used to retrieve the current container.
 * @param ListInner All particles inside the bounding box.
 * @param ListHaloWithinCutoff All particles in the halo.
 * @param ListHaloOutsideCutoff All particles in the halo.
 */
void getStatus(const std::array<double, 3> &bBoxMin, const std::array<double, 3> &bBoxMax, const double cutoff,
               autopas::ContainerSelector<Particle> &containerSelector, std::vector<Particle> &ListInner,
               std::vector<Particle> &ListHaloWithinCutoff, std::vector<Particle> &ListHaloOutsideCutoff) {
  using namespace autopas::utils::ArrayMath::literals;

  for (auto iter = containerSelector.getCurrentContainer()->begin(autopas::IteratorBehavior::owned); iter.isValid();
       ++iter) {
    ListInner.push_back(*iter);
  }
  const auto cutoffBoxMin = bBoxMin - cutoff;
  const auto cutoffBoxMax = bBoxMax + cutoff;
  for (auto iter = containerSelector.getCurrentContainer()->begin(autopas::IteratorBehavior::halo); iter.isValid();
       ++iter) {
    if (autopas::utils::inBox(iter->getR(), cutoffBoxMin, cutoffBoxMax)) {
      ListHaloWithinCutoff.push_back(*iter);
    } else {
      ListHaloOutsideCutoff.push_back(*iter);
    }
  }
}

TEST_P(ContainerSelectorTestFromTo, testContainerConversion) {
  const auto &[from, to] = GetParam();

  autopas::ContainerSelector<Particle> containerSelector(bBoxMin, bBoxMax, cutoff);
  autopas::ContainerSelectorInfo containerInfo(cellSizeFactor, verletSkinPerTimestep, verletRebuildFrequency, 64,
                                               autopas::LoadEstimatorOption::none);

  // select container from which we want to convert from
  containerSelector.selectContainer(from, containerInfo);

  // fill with problematic particles
  {
    auto container = containerSelector.getCurrentContainer();
    auto getPossible1DPositions = [&](double min, double max) -> auto{
      return std::array<double, 6>{min - cutoff - verletSkinPerTimestep * verletRebuildFrequency,
                                   min - cutoff,
                                   min,
                                   max,
                                   max + cutoff - 1e-3,
                                   max + cutoff + verletSkinPerTimestep * verletRebuildFrequency - 1e-3};
    };
    size_t id = 0;

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

  std::vector<Particle> beforeListInner, beforeListHaloWithinCutoff,
      beforeListHaloOutsideCutoff /*for particles only in verlet containers*/;

  getStatus(bBoxMin, bBoxMax, cutoff, containerSelector, beforeListInner, beforeListHaloWithinCutoff,
            beforeListHaloOutsideCutoff);

  // select container to which we want to convert to
  containerSelector.selectContainer(to, containerInfo);

  std::vector<Particle> afterListInner, afterListHaloWithinCutoff, afterListHaloOutsideCutoff;

  getStatus(bBoxMin, bBoxMax, cutoff, containerSelector, afterListInner, afterListHaloWithinCutoff,
            afterListHaloOutsideCutoff);

  EXPECT_EQ(afterListInner.size(), beforeListInner.size());
  EXPECT_EQ(afterListHaloWithinCutoff.size(), beforeListHaloWithinCutoff.size());
  EXPECT_EQ(afterListHaloOutsideCutoff.size(), beforeListHaloOutsideCutoff.size());

  // abort the test if any list size does not match as element comparisons do not work then.
  if (::testing::Test::HasFailure()) {
    FAIL();
  }

  EXPECT_THAT(afterListInner, UnorderedElementsAreArray(beforeListInner));
  EXPECT_THAT(afterListHaloWithinCutoff, UnorderedElementsAreArray(beforeListHaloWithinCutoff));
  EXPECT_THAT(afterListHaloOutsideCutoff, UnorderedElementsAreArray(beforeListHaloOutsideCutoff));
}

INSTANTIATE_TEST_SUITE_P(Generated, ContainerSelectorTestFromTo,
                         Combine(ValuesIn(autopas::ContainerOption::getAllOptions()),
                                 ValuesIn(autopas::ContainerOption::getAllOptions())),
                         ContainerSelectorTestFromTo::twoParamToString());
