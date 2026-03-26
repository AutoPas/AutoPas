/**
 * @file Newton3OnOffTest.cpp
 * @author seckler
 * @date 18.04.18
 */

#include "Newton3OnOffTest.h"

#include "autopas/tuning/selectors/ContainerSelector.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "autopas/utils/logging/Logger.h"
#include "autopasTools/generators/UniformGenerator.h"

using ::testing::_;
using ::testing::Combine;
using ::testing::Return;
using ::testing::ValuesIn;

// Parse combination strings and call actual test function
TEST_P(Newton3OnOffTest, countFunctorCallsTest) {
  auto config = GetParam();

  if (config.interactionType == autopas::InteractionTypeOption::pairwise) {
    countFunctorCalls<MPairwiseFunctor>(config);
  } else if (config.interactionType == autopas::InteractionTypeOption::triwise) {
    countFunctorCalls<MTriwiseFunctor>(config);
  }
}

// Generate Unittests for all Configurations that support both Newton3 modes
INSTANTIATE_TEST_SUITE_P(Generated, Newton3OnOffTest, ValuesIn([]() -> std::set<autopas::Configuration> {
                           // Start with all valid configurations
                           auto configs = generateAllValidConfigurations(autopas::InteractionTypeOption::all);

                           // Remove all configs that do not support both Newton3 modes
                           std::erase_if(configs, [&](const autopas::Configuration &config) {
                             const auto configsSupportingN3EnabledOnly =
                                 autopas::compatibleTraversals::allTraversalsSupportingOnlyNewton3Enabled();
                             const auto configsSupportingN3DisabledOnly =
                                 autopas::compatibleTraversals::allTraversalsSupportingOnlyNewton3Disabled();
                             return configsSupportingN3EnabledOnly.contains(config.traversal) or
                                    configsSupportingN3DisabledOnly.contains(config.traversal);
                           });

                           // Remove all newton3 enabled configurations so that for each config-N3 enabled/disabled
                           // pair, only one config is provided, from which the other can be generated.
                           std::erase_if(configs, [&](auto &config) {
                             return config.newton3 == autopas::Newton3Option::enabled;
                           });

                           return configs;
                         }()),
                         Newton3OnOffTest::PrintToStringParamName());

template <typename Functor_T>
void Newton3OnOffTest::countFunctorCalls(autopas::Configuration config) {
  // TODO: Make test possible for direct sum SoA
  if (config.container == autopas::ContainerOption::directSum and config.dataLayout == autopas::DataLayoutOption::soa) {
    return;
  }
  const autopas::ContainerSelectorInfo containerInfo(getBoxMin(), getBoxMax(), getCutoff(), config.cellSizeFactor,
                                                     getVerletSkin(), getClusterSize(), getSortingThreshold(),
                                                     config.loadEstimator);
  auto container = autopas::ContainerSelector<ParticleFP64>::generateContainer(config.container, containerInfo);

  const ParticleFP64 defaultParticle{};
  autopasTools::generators::UniformGenerator::fillWithParticles(*container, defaultParticle, container->getBoxMin(),
                                                                container->getBoxMax(), 200);
  // Do not add any halo particles to this test!
  // Given an owned particle p1 and a halo particle p2 the following interactions are necessary:
  // Newton3   : p1 <-> p2
  // No Newton3: p1 <-  p2 (but NOT p1 -> p2)
  // Depending if the container is able to avoid all halo <- owned interactions, the total number of functor calls
  // is not trivially to calculate.

  // Create the functor
  Functor_T mockFunctor{};

  EXPECT_CALL(mockFunctor, isRelevantForTuning()).WillRepeatedly(Return(true));

  if (config.dataLayout == autopas::DataLayoutOption::soa) {
    if (config.container == autopas::ContainerOption::linkedCellsReferences) {
      // loader and extractor will be called, we don't care how often.
      auto *expectation =
          &EXPECT_CALL(mockFunctor,
                       SoALoader(testing::Matcher<autopas::ReferenceParticleCell<ParticleFP64> &>(_), _, _, _))
               .Times(testing::AtLeast(1));
      expectation->WillRepeatedly(
          testing::WithArgs<0, 1>(testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.size()); })));
      EXPECT_CALL(mockFunctor, SoAExtractor(testing::Matcher<autopas::ReferenceParticleCell<ParticleFP64> &>(_), _, _))
          .Times(testing::AtLeast(1));
    } else {
      // loader and extractor will be called, we don't care how often.
      auto *expectation =
          &EXPECT_CALL(mockFunctor, SoALoader(testing::Matcher<autopas::FullParticleCell<ParticleFP64> &>(_), _, _, _))
               .Times(testing::AtLeast(1));
      // Verlet based containers resize the SoA before they call SoALoader, so no need for the testing::Invoke here
      if (std::set{autopas::ContainerOption::varVerletListsAsBuild, autopas::ContainerOption::pairwiseVerletLists,
                   autopas::ContainerOption::verletListsCells}
              .count(container->getContainerType()) == 0) {
        expectation->WillRepeatedly(
            testing::WithArgs<0, 1>(testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.size()); })));
      }

      EXPECT_CALL(mockFunctor, SoAExtractor(testing::Matcher<autopas::FullParticleCell<ParticleFP64> &>(_), _, _))
          .Times(testing::AtLeast(1));
    }
  }

  // "SC" = single cell
  const auto [callsNewton3SC, callsNewton3Pair, callsNewton3Triple] =
      eval(config, /*useNewton3*/ true, *container, mockFunctor);
  const auto [callsNonNewton3SC, callsNonNewton3Pair, callsNonNewton3Triple] =
      eval(config, /*useNewton3*/ false, *container, mockFunctor);

  EXPECT_GT(callsNewton3Triple + callsNewton3Pair + callsNewton3SC, 0)
      << "No interactions were found WITH Newton3. Are the generated particles too far apart?";
  EXPECT_GT(callsNonNewton3Triple + callsNonNewton3Pair + callsNonNewton3SC, 0)
      << "No interactions were found WITHOUT Newton3. Are the generated particles too far apart?";

  if (config.dataLayout == autopas::DataLayoutOption::soa) {
    // within one cell no N3 optimization
    EXPECT_EQ(callsNewton3SC, callsNonNewton3SC)
        << "Mismatch in number of interactions within cells for container option: " << config.container;
    // two times the calls for no N3, once for each "owning" cell
    EXPECT_EQ(callsNewton3Pair * 2, callsNonNewton3Pair)
        << "Mismatch in number of interactions between cell pairs for container option: " << config.container;
    // three times the calls for no N3, once for each "owning" cell
    EXPECT_EQ(callsNewton3Triple * 3, callsNonNewton3Triple)
        << "Mismatch in number of interactions between cell triplets for container option: " << config.container;

  } else {  // aos
    // should be called exactly two times
    EXPECT_EQ(callsNewton3Pair * 2, callsNonNewton3Pair)
        << "Mismatch in number of interactions between cell pairs for container option: " << config.container;
    // should be called exactly three times
    EXPECT_EQ(callsNewton3Triple * 3, callsNonNewton3Triple) << "for container option: " << config.container;
  }

  if (HasFailure()) {
    std::cerr << "Failures for Container: " << config.container.to_string()
              << ", Traversal: " << config.traversal.to_string() << ", Data Layout: " << config.dataLayout.to_string()
              << std::endl;
  }
}

template <class Container_T, class Traversal_T>
void Newton3OnOffTest::iterate(Container_T &container, Traversal_T traversal) {
  container.rebuildNeighborLists(traversal.get());
  container.computeInteractions(traversal.get());
}

template <class Functor_T, class Container_T>
std::tuple<size_t, size_t, size_t> Newton3OnOffTest::eval(autopas::Configuration config, bool useNewton3,
                                                          Container_T &container, Functor_T &mockFunctor) {
  // set up counters for particle interactions
  std::atomic<unsigned int> callsSC(0ul);
  std::atomic<unsigned int> callsPair(0ul);
  std::atomic<unsigned int> callsTriple(0ul);

  auto getInteractionTypeFromFunctor = [&]() {
    if constexpr (autopas::utils::isPairwiseFunctor<Functor_T>()) {
      return autopas::InteractionTypeOption::pairwise;
    } else if constexpr (autopas::utils::isTriwiseFunctor<Functor_T>()) {
      return autopas::InteractionTypeOption::triwise;
    }
  };
  constexpr auto interactionType = getInteractionTypeFromFunctor();

  EXPECT_CALL(mockFunctor, allowsNewton3()).WillRepeatedly(Return(useNewton3));
  EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillRepeatedly(Return(not useNewton3));

  // depending on the layout, set up some expectations on what Functors are called
  switch (config.dataLayout) {
    case autopas::DataLayoutOption::soa: {
      // some containers actually use different SoA functor calls so expect them instead of the regular ones
      if (container.getContainerType() == autopas::ContainerOption::varVerletListsAsBuild ||
          container.getContainerType() == autopas::ContainerOption::pairwiseVerletLists ||
          container.getContainerType() == autopas::ContainerOption::verletListsCells) {
        EXPECT_CALL(mockFunctor, SoAFunctorVerlet(_, _, _, useNewton3))
            .Times(testing::AtLeast(1))
            // signature: SoAFunctorVerlet(soa, particleIndex, neighbors, useNewton3)
            .WillRepeatedly(testing::Invoke([&](auto, auto, auto neighbors, auto) { callsPair += neighbors.size(); }));

        // non useNewton3 variant should not happen
        EXPECT_CALL(mockFunctor, SoAFunctorVerlet(_, _, _, not useNewton3)).Times(0);

      } else {
        // single cell
        EXPECT_CALL(mockFunctor, SoAFunctorSingle(_, useNewton3))
            .Times(testing::AtLeast(1))
            // signature: SoAFunctorSingle(soa, newton3)
            .WillRepeatedly(testing::Invoke([&](auto soa, auto) { callsSC += soa.size(); }));

        // pair of cells
        EXPECT_CALL(mockFunctor, SoAFunctorPair(_, _, useNewton3))
            .Times(testing::AtLeast(1))
            // signature: SoAFunctorPair(soa1, soa2, newton3)
            .WillRepeatedly(
                testing::Invoke([&](auto soa1, auto soa2, auto) { callsPair += soa1.size() * soa2.size(); }));

        // non useNewton3 variant should not happen
        EXPECT_CALL(mockFunctor, SoAFunctorPair(_, _, not useNewton3)).Times(0);

        if constexpr (interactionType == autopas::InteractionTypeOption::triwise) {
          // triplet of cells
          EXPECT_CALL(mockFunctor, SoAFunctorTriple(_, _, _, useNewton3))
              .Times(testing::AtLeast(0))  // Direct sum has only 2 cells
              .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsTriple++; }));

          // non useNewton3 variant should not happen
          EXPECT_CALL(mockFunctor, SoAFunctorTriple(_, _, _, not useNewton3)).Times(0);
        }
      }
      break;
    }
    case autopas::DataLayoutOption::aos: {
      if constexpr (interactionType == autopas::InteractionTypeOption::pairwise) {
        EXPECT_CALL(mockFunctor, AoSFunctor(_, _, useNewton3))
            .Times(testing::AtLeast(1))
            .WillRepeatedly(testing::InvokeWithoutArgs([&]() { ++callsPair; }));

        // non useNewton3 variant should not happen
        EXPECT_CALL(mockFunctor, AoSFunctor(_, _, not useNewton3)).Times(0);
      } else if (interactionType == autopas::InteractionTypeOption::triwise) {
        EXPECT_CALL(mockFunctor, AoSFunctor(_, _, _, useNewton3))
            .Times(testing::AtLeast(1))
            .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsTriple++; }));

        // non useNewton3 variant should not happen
        EXPECT_CALL(mockFunctor, AoSFunctor(_, _, _, not useNewton3)).Times(0);
      }

      break;
    }
    default: {
      ADD_FAILURE() << "This test does not support data layout : " << config.dataLayout.to_string();
    }
  }

  // We expect a newton3 disabled configuration, from which we need to generate a newton3 enabled one if useNewton3 is
  // true
  if (useNewton3) {
    config.newton3 = autopas::Newton3Option::enabled;
  }

  auto traversalSelectorInfo = container.getTraversalSelectorInfo();

  // simulate iteration
  iterate(container, autopas::TraversalSelector::generateTraversalFromConfig<ParticleFP64, Functor_T>(
                         config, mockFunctor, traversalSelectorInfo));

  return std::make_tuple(callsSC.load(), callsPair.load(), callsTriple.load());
}
