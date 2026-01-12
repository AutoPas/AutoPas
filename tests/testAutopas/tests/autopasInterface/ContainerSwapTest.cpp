/**
 * @file ContainerSwapTest.cpp
 * @author muehlhaeusser
 * @date 30.07.2025
 */

#include "ContainerSwapTest.h"

#include "autopas/LogicHandler.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/particles/OwnershipState.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::Combine;
using ::testing::Return;
using ::testing::UnorderedElementsAreArray;
using ::testing::ValuesIn;

/**
 * This function stores a copy of each particle depending on the position in ListInner, ListHaloWithinCutoff or
 * ListHaloOutsideCutoff.
 * @tparam Container
 * @param bBoxMin Bounding box min.
 * @param bBoxMax Bounding box max.
 * @param cutoff Cutoff radius.
 * @param container Container selector used to retrieve the current container.
 * @param ListInner All particles inside the bounding box.
 * @param ListHaloWithinCutoff All particles in the halo but within cutoff of the bounding box.
 * @param ListHaloOutsideCutoff All particles in the halo and outside the cutoff of the bounding box.
 */
template <class Container>
void gatherContainerParticles(const std::array<double, 3> &bBoxMin, const std::array<double, 3> &bBoxMax,
                              const double cutoff, Container &container, std::vector<ParticleFP64> &ListInner,
                              std::vector<ParticleFP64> &ListHaloWithinCutoff,
                              std::vector<ParticleFP64> &ListHaloOutsideCutoff) {
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

/**
 * This tests whether the logic handler swaps between different kind of containers without losing particles.
 * It initializes a searchspace of two configs and swaps between the first and second config back and forth.
 */
TEST_P(ContainerSwapTest, testContainerConversion) {
  const auto &[config1, config2] = GetParam();

  constexpr autopas::DataLayoutOption dataLayout = autopas::DataLayoutOption::aos;
  constexpr autopas::Newton3Option newton3 = autopas::Newton3Option::disabled;
  auto config1TraversalOptions = autopas::compatibleTraversals::allCompatibleTraversals(
      config1.container, autopas::InteractionTypeOption::pairwise);
  auto config2TraversalOptions = autopas::compatibleTraversals::allCompatibleTraversals(
      config2.container, autopas::InteractionTypeOption::pairwise);

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

  const std::set searchSpace({config1, config2});

  std::unordered_map<autopas::InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> tunerMap;
  tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  autopas::LogicHandler<ParticleFP64> logicHandler(tunerMap, logicHandlerInfo, verletRebuildFrequency, "");

  // Helper to add particles to the container.
  auto addParticlesToContainer = [&](auto &containerToFill) {
    auto getPossible1DPositions = [&](double min, double max) -> auto{
      return std::array<double, 6>{min - cutoff - verletSkin,       min - cutoff, min, max, max + cutoff - 1e-3,
                                   max + cutoff + verletSkin - 1e-3};
    };
    containerToFill.deleteAllParticles();
    size_t id = 0;
    for (auto x : getPossible1DPositions(bBoxMin[0], bBoxMax[0])) {
      for (auto y : getPossible1DPositions(bBoxMin[1], bBoxMax[1])) {
        for (auto z : getPossible1DPositions(bBoxMin[2], bBoxMax[2])) {
          const std::array<double, 3> pos{x, y, z};
          ParticleFP64 p(pos, {0., 0., 0.}, id);
          if (autopas::utils::inBox(pos, bBoxMin, bBoxMax)) {
            containerToFill.addParticle(p);
          } else {
            p.setOwnershipState(autopas::OwnershipState::halo);
            containerToFill.addHaloParticle(p);
          }
          ++id;
        }
      }
    }
  };

  // The first computeInteractions run to initialize the 'from' container.
  MockPairwiseFunctor<ParticleFP64> functor{};
  EXPECT_CALL(functor, isRelevantForTuning()).WillRepeatedly(Return(true));
  EXPECT_CALL(functor, allowsNewton3()).WillRepeatedly(Return(true));
  EXPECT_CALL(functor, allowsNonNewton3()).WillRepeatedly(Return(true));
  EXPECT_CALL(functor, isVecPatternAllowed(::testing::_)).WillRepeatedly(::testing::Return(true));
  logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
  const auto firstContainerType = logicHandler.getContainer().getContainerType();
  const auto secondContainerType =
      logicHandler.getContainer().getContainerType() == config1.container ? config2.container : config1.container;

  // Start the second iteration which should swap the container to the secondContainerType configuration.
  auto emigrants = logicHandler.updateContainer();
  ASSERT_TRUE(emigrants.empty()) << "There should be no emigrating particles in this test.";
  addParticlesToContainer(logicHandler.getContainer());

  std::vector<ParticleFP64> beforeListInner, beforeListHaloWithinCutoff,
      beforeListHaloOutsideCutoff /*for particles only in verlet containers*/;
  gatherContainerParticles(bBoxMin, bBoxMax, cutoff, logicHandler.getContainer(), beforeListInner,
                           beforeListHaloWithinCutoff, beforeListHaloOutsideCutoff);

  // Container swap should happen during computeInteractions
  logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
  ASSERT_EQ(logicHandler.getContainer().getContainerType(), secondContainerType);

  std::vector<ParticleFP64> afterListInner, afterListHaloWithinCutoff, afterListHaloOutsideCutoff;
  gatherContainerParticles(bBoxMin, bBoxMax, cutoff, logicHandler.getContainer(), afterListInner,
                           afterListHaloWithinCutoff, afterListHaloOutsideCutoff);

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

  // Reset tuning, so third iteration should swap back to firstContainerType
  tunerMap[autopas::InteractionTypeOption::pairwise]->forceRetune();
  emigrants = logicHandler.updateContainer();
  ASSERT_TRUE(emigrants.empty()) << "There should be no emigrating particles in this test.";

  addParticlesToContainer(logicHandler.getContainer());

  logicHandler.computeInteractionsPipeline(&functor, autopas::InteractionTypeOption::pairwise);
  ASSERT_EQ(logicHandler.getContainer().getContainerType(), firstContainerType);

  std::vector<ParticleFP64> after2ListInner, after2ListHaloWithinCutoff, after2ListHaloOutsideCutoff;
  gatherContainerParticles(bBoxMin, bBoxMax, cutoff, logicHandler.getContainer(), after2ListInner,
                           after2ListHaloWithinCutoff, after2ListHaloOutsideCutoff);

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
    {autopas::ContainerOption::directSum, 1, autopas::TraversalOption::ds_sequential,
     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled,
     autopas::InteractionTypeOption::pairwise},
    {autopas::ContainerOption::linkedCells, 1, autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
     autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise},
    {autopas::ContainerOption::linkedCellsReferences, 1, autopas::TraversalOption::lc_c01,
     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled,
     autopas::InteractionTypeOption::pairwise},
    {autopas::ContainerOption::verletLists, 1, autopas::TraversalOption::vl_list_iteration,
     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled,
     autopas::InteractionTypeOption::pairwise},
    {autopas::ContainerOption::varVerletListsAsBuild, 1, autopas::TraversalOption::vvl_as_built,
     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled,
     autopas::InteractionTypeOption::pairwise},
    {autopas::ContainerOption::verletClusterLists, 1, autopas::TraversalOption::vcl_cluster_iteration,
     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled,
     autopas::InteractionTypeOption::pairwise},
    {autopas::ContainerOption::verletListsCells, 1, autopas::TraversalOption::vlc_c01,
     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled,
     autopas::InteractionTypeOption::pairwise},
    {autopas::ContainerOption::pairwiseVerletLists, 1, autopas::TraversalOption::vlp_c01,
     autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled,
     autopas::InteractionTypeOption::pairwise},
    {autopas::ContainerOption::octree, 1, autopas::TraversalOption::ot_c01, autopas::LoadEstimatorOption::none,
     autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled, autopas::InteractionTypeOption::pairwise}};

// Generates all unique pairs of configurations, order does not matter and no pairs of the same configuration.
std::vector<std::pair<autopas::Configuration, autopas::Configuration>> GenerateUniquePairs(
    const std::vector<autopas::Configuration> &configs) {
  // Check that all container options are covered.
  std::set<autopas::ContainerOption> givenConfigs;
  for (const auto &config : configs) {
    givenConfigs.insert(config.container);
  }
  if (givenConfigs != autopas::ContainerOption::getAllOptions()) {
    throw std::runtime_error("ContainerSwapTest: Given configurations do not cover all container options!");
  }

  // Generate all unique pairs.
  std::vector<std::pair<autopas::Configuration, autopas::Configuration>> pairs;
  for (size_t i = 0; i < configs.size(); ++i) {
    for (size_t j = i + 1; j < configs.size(); ++j) {
      pairs.emplace_back(configs[i], configs[j]);
    }
  }
  return pairs;
}

INSTANTIATE_TEST_SUITE_P(Generated, ContainerSwapTest, ValuesIn(GenerateUniquePairs(containerConfigs)),
                         ContainerSwapTest::twoParamToString());
