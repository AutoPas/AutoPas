/**
 * @file AllContainersTests.h
 * @author humig
 * @date 08.07.2019
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/tuning/selectors/ContainerSelector.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "testingHelpers/EmptyPairwiseFunctor.h"
#include "testingHelpers/GenerateValidConfigurations.h"
#include "testingHelpers/commonTypedefs.h"

class AllContainersTestsBase : public AutoPasTestBase {
 protected:
  std::array<double, 3> boxMin = {0, 0, 0};
  std::array<double, 3> boxMax = {10, 10, 10};
  double cutoff = 1;
  const double skin = 0.2;
  const unsigned int rebuildFrequency = 20;

  template <class Particle_T>
  auto getInitializedContainer(const ContainerConfiguration &containerConfig) {
    const autopas::ContainerSelectorInfo selectorInfo{boxMin, boxMax, cutoff, containerConfig.cellSizeFactor,    skin,
                                                      32,     8,      8,      autopas::LoadEstimatorOption::none};
    auto container = autopas::ContainerSelector<Particle_T>::generateContainer(containerConfig.container, selectorInfo);
    return std::move(container);
  }

  /**
   * Builds the internal structure of the container as the LogicHandler does before the first traversal after a rebuild.
   * E.g. VerletClusterLists sorts the added particles into its towers.
   * @tparam Particle_T
   * @param container
   * @param containerConfig The configuration the container was initialized with.
   */
  template <class Particle_T>
  void rebuildNeighborLists(autopas::ParticleContainerInterface<Particle_T> &container,
                            const ContainerConfiguration &containerConfig) {
    const auto config = containerConfig.generateFullConfig(autopas::InteractionTypeOption::pairwise).value();
    EmptyPairwiseFunctor<Particle_T> functor;
    const auto traversal = autopas::TraversalSelector::generateTraversalFromConfig<Particle_T>(
        config, functor, container.getTraversalSelectorInfo());
    ASSERT_NE(traversal, nullptr) << "No traversal could be generated for " << config.toString();
    container.rebuildNeighborLists(traversal.get());
  }
};

using ParamType = std::tuple<ContainerConfiguration>;

class AllContainersTests : public AllContainersTestsBase, public ::testing::WithParamInterface<ParamType> {
 public:
  static auto getParamToStringFunction() {
    static const auto paramToString = [](const testing::TestParamInfo<ParamType> &info) {
      auto [containerConfig] = info.param;
      return containerConfig.toShortString();
    };
    return paramToString;
  }
};

using ParamTypeBothUpdates = std::tuple<ContainerConfiguration, bool /*keep Lists Valid*/>;

class AllContainersTestsBothUpdates : public AllContainersTestsBase,
                                      public ::testing::WithParamInterface<ParamTypeBothUpdates> {
 public:
  static auto getParamToStringFunction() {
    static const auto paramToString = [](const testing::TestParamInfo<ParamType> &info) {
      auto [containerConfig, keepListValid] = info.param;
      return containerConfig.toShortString() + "_" + (keepListValid ? "keepListsValid" : "allowListInvalidation");
    };
    return paramToString;
  }

 protected:
  void testUpdateContainerDeletesDummy(bool previouslyOwned);
};