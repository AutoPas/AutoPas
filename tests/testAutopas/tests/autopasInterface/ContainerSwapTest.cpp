/**
 * @file ContainerSwap.cpp
 * @author F. Gratl
 * @date 14.12.2020
 */

#include "ContainerSwapTest.h"

#include "autopas/LogicHandler.h"
#include "autopas/particles/OwnershipState.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::ValuesIn;

namespace autopas {
/**
 * This function stores a copy of each particle depending on the position in ListInner, ListHaloWithinCutoff or
 * ListHaloOutsideCutoff.
 * @tparam Container
 * @param bBoxMin Bounding box min.
 * @param bBoxMax Bounding box max.
 * @param cutoff Cutoff radius.
 * @param container Container selector used to retrieve the current container.
 * @param ListInner All particles inside the bounding box.
 * @param ListHaloWithinCutoff All particles in the halo.
 * @param ListHaloOutsideCutoff All particles in the halo.
 */
template <class Container>
void getStatus(const std::array<double, 3> &bBoxMin, const std::array<double, 3> &bBoxMax, const double cutoff,
               Container &container, std::vector<ParticleFP64> &ListInner,
               std::vector<ParticleFP64> &ListHaloWithinCutoff, std::vector<ParticleFP64> &ListHaloOutsideCutoff) {
  using namespace autopas::utils::ArrayMath::literals;

  for (auto iter = container.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    ListInner.push_back(*iter);
  }
  const auto cutoffBoxMin = bBoxMin - cutoff;
  const auto cutoffBoxMax = bBoxMax + cutoff;
  for (auto iter = container.begin(autopas::IteratorBehavior::halo); iter.isValid(); ++iter) {
    if (autopas::utils::inBox(iter->getR(), cutoffBoxMin, cutoffBoxMax)) {
      ListHaloWithinCutoff.push_back(*iter);
    } else {
      ListHaloOutsideCutoff.push_back(*iter);
    }
  }
}

TEST_P(ContainerSwapTest, testContainerConversion) {
  const auto &[from, to] = GetParam();

  constexpr autopas::DataLayoutOption dataLayout = autopas::DataLayoutOption::aos;
  constexpr autopas::Newton3Option newton3 = autopas::Newton3Option::disabled;

  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{bBoxMin},
      .boxMax{bBoxMax},
      .cutoff = cutoff,
      .verletSkin = verletSkin,
  };
  constexpr autopas::AutoTunerInfo autoTunerInfo{};
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const std::set<autopas::Configuration> searchSpace(
      {{from, cellSizeFactor, autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none, dataLayout, newton3,
        autopas::InteractionTypeOption::pairwise}});

  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<ParticleFP64> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");

  // fill with problematic particles
  auto &container = logicHandler.getContainer();
  auto getPossible1DPositions = [&](double min, double max) -> auto{
    return std::array<double, 6>{min - cutoff - verletSkin,       min - cutoff, min, max, max + cutoff - 1e-3,
                                 max + cutoff + verletSkin - 1e-3};
  };
  size_t id = 0;

  for (auto x : getPossible1DPositions(bBoxMin[0], bBoxMax[0])) {
    for (auto y : getPossible1DPositions(bBoxMin[1], bBoxMax[1])) {
      for (auto z : getPossible1DPositions(bBoxMin[2], bBoxMax[2])) {
        const std::array<double, 3> pos{x, y, z};
        ParticleFP64 p(pos, {0., 0., 0.}, id);
        if (autopas::utils::inBox(pos, bBoxMin, bBoxMax)) {
          container.addParticle(p);
        } else {
          p.setOwnershipState(autopas::OwnershipState::halo);
          container.addHaloParticle(p);
        }
        ++id;
      }
    }
  }

  std::vector<ParticleFP64> beforeListInner, beforeListHaloWithinCutoff,
      beforeListHaloOutsideCutoff /*for particles only in verlet containers*/;

  getStatus(bBoxMin, bBoxMax, cutoff, container, beforeListInner, beforeListHaloWithinCutoff,
            beforeListHaloOutsideCutoff);

  // select container to which we want to convert to
  autopas::ContainerSelectorInfo containerInfo(bBoxMin, bBoxMax, cutoff, cellSizeFactor, verletSkin, 64, 64,
                                               autopas::LoadEstimatorOption::none);
  auto newContainer = autopas::ContainerSelector<ParticleFP64>::generateContainer(to, containerInfo);
  logicHandler.setCurrentContainer(std::move(newContainer));

  std::vector<ParticleFP64> afterListInner, afterListHaloWithinCutoff, afterListHaloOutsideCutoff;

  getStatus(bBoxMin, bBoxMax, cutoff, logicHandler.getContainer(), afterListInner, afterListHaloWithinCutoff,
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

INSTANTIATE_TEST_SUITE_P(Generated, ContainerSwapTest,
                         Combine(ValuesIn(autopas::ContainerOption::getAllOptions()),
                                 ValuesIn(autopas::ContainerOption::getAllOptions())),
                         ContainerSwapTest::twoParamToString());
}  // namespace autopas
