/**
 * @file ContainerSwapTest.cpp
 * @author F. Gratl
 * @date 14.12.2020
 */

#include "ContainerSwapTest.h"

#include "autopas/LogicHandler.h"
#include "autopas/particles/OwnershipState.h"
#include "testingHelpers/commonTypedefs.h"
#include "autopas/containers/CompatibleTraversals.h"

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::ValuesIn;
using ::testing::Return;

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
  auto fromTraversalOptions = autopas::compatibleTraversals::allCompatibleTraversals(from.container, autopas::InteractionTypeOption::pairwise);
  auto toTraversalOptions = autopas::compatibleTraversals::allCompatibleTraversals(to.container, autopas::InteractionTypeOption::pairwise);

  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{bBoxMin},
      .boxMax{bBoxMax},
      .cutoff = cutoff,
      .verletSkin = verletSkin,
  };
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .maxSamples = 1,
      };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const std::set searchSpace({from, to});

  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<ParticleFP64> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");

  // The first computeInteractions run to initialize the 'from' container.
  MockPairwiseFunctor<ParticleFP64> functor{};
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(Return(true));
  logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
  ASSERT_EQ(logicHandler.getContainer().getContainerType(), from.container);

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


  // Second iteration should go for the second configuration
  static_cast<void>(logicHandler.updateContainer());
  logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
  ASSERT_EQ(logicHandler.getContainer().getContainerType(), to.container);

  std::vector<ParticleFP64> afterListInner, afterListHaloWithinCutoff, afterListHaloOutsideCutoff;
  getStatus(bBoxMin, bBoxMax, cutoff, logicHandler.getContainer(), afterListInner, afterListHaloWithinCutoff,
            afterListHaloOutsideCutoff);

  EXPECT_EQ(afterListInner.size(), beforeListInner.size());
  EXPECT_EQ(afterListHaloWithinCutoff.size(), beforeListHaloWithinCutoff.size());
  EXPECT_EQ(afterListHaloOutsideCutoff.size(), beforeListHaloOutsideCutoff.size());

  // abort the test if any list size does not match as element comparisons do not work then.
  if (HasFailure()) {
    FAIL();
  }

  EXPECT_THAT(afterListInner, UnorderedElementsAreArray(beforeListInner));
  EXPECT_THAT(afterListHaloWithinCutoff, UnorderedElementsAreArray(beforeListHaloWithinCutoff));
  EXPECT_THAT(afterListHaloOutsideCutoff, UnorderedElementsAreArray(beforeListHaloOutsideCutoff));

  // Third iteration should go for the first configuration again
  tunerMap[autopas::InteractionTypeOption::pairwise]->forceRetune();
  static_cast<void>(logicHandler.updateContainer());
  logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
  ASSERT_EQ(logicHandler.getContainer().getContainerType(), from.container);

  std::vector<ParticleFP64> after2ListInner, after2ListHaloWithinCutoff, after2ListHaloOutsideCutoff;
  getStatus(bBoxMin, bBoxMax, cutoff, logicHandler.getContainer(), after2ListInner, after2ListHaloWithinCutoff,
            after2ListHaloOutsideCutoff);

  EXPECT_EQ(after2ListInner.size(), afterListInner.size());
  EXPECT_EQ(after2ListHaloWithinCutoff.size(), afterListHaloWithinCutoff.size());
  EXPECT_EQ(after2ListHaloOutsideCutoff.size(), afterListHaloOutsideCutoff.size());

  // abort the test if any list size does not match as element comparisons do not work then.
  if (HasFailure()) {
    FAIL();
  }

  EXPECT_THAT(after2ListInner, UnorderedElementsAreArray(afterListInner));
  EXPECT_THAT(after2ListHaloWithinCutoff, UnorderedElementsAreArray(afterListHaloWithinCutoff));
  EXPECT_THAT(after2ListHaloOutsideCutoff, UnorderedElementsAreArray(afterListHaloOutsideCutoff));
}

std::vector<autopas::Configuration> containerConfigs = {
  {autopas::ContainerOption::directSum, 1, autopas::TraversalOption::ds_sequential, autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise},
{autopas::ContainerOption::linkedCells, 1, autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise},
{autopas::ContainerOption::linkedCellsReferences, 1, autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise},
{autopas::ContainerOption::verletLists, 1, autopas::TraversalOption::vl_list_iteration, autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise},
{autopas::ContainerOption::varVerletListsAsBuild, 1, autopas::TraversalOption::vvl_as_built, autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise},
  {autopas::ContainerOption::verletClusterLists, 1, autopas::TraversalOption::vcl_cluster_iteration, autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise},
  {autopas::ContainerOption::verletListsCells, 1, autopas::TraversalOption::vlc_c01, autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise},
  {autopas::ContainerOption::pairwiseVerletLists, 1, autopas::TraversalOption::vlp_c01, autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise},
{autopas::ContainerOption::octree, 1, autopas::TraversalOption::ot_c01, autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise}};

std::vector<std::pair<autopas::Configuration, autopas::Configuration>> GenerateUniquePairs(
    const std::vector<autopas::Configuration>& configs
) {
  std::vector<std::pair<autopas::Configuration, autopas::Configuration>> pairs;
  for (size_t i = 0; i < configs.size(); ++i) {
    for (size_t j = i + 1; j < configs.size(); ++j) {
      pairs.emplace_back(configs[i], configs[j]);
    }
  }
  return pairs;
}

INSTANTIATE_TEST_SUITE_P(Generated, ContainerSwapTest,
                         ValuesIn(GenerateUniquePairs(containerConfigs)),
                         ContainerSwapTest::twoParamToString());
