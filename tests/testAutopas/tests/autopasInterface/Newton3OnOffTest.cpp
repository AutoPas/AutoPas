/**
 * @file Newton3OnOffTest.cpp
 * @author seckler
 * @date 18.04.18
 */

#include "Newton3OnOffTest.h"

#include "autopas/tuning/selectors/ContainerSelector.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "autopas/utils/StaticCellSelector.h"
#include "autopas/utils/logging/Logger.h"
#include "autopasTools/generators/UniformGenerator.h"

using ::testing::_;
using ::testing::Combine;
using ::testing::Return;
using ::testing::ValuesIn;

// Parse combination strings and call actual test function
TEST_P(Newton3OnOffTest, countFunctorCallsTest) {
  auto [containerOption, traversalOption, dataLayoutOption, interactionTypeOption] = GetParam();

  if (interactionTypeOption == autopas::InteractionTypeOption::pairwise) {
    countFunctorCalls<MPairwiseFunctor>(containerOption, traversalOption, dataLayoutOption, interactionTypeOption);
  } else if (interactionTypeOption == autopas::InteractionTypeOption::triwise) {
    countFunctorCalls<MTriwiseFunctor>(containerOption, traversalOption, dataLayoutOption, interactionTypeOption);
  }
}

// Generate Unittests for all Configurations that support both Newton3 modes
INSTANTIATE_TEST_SUITE_P(
    Generated, Newton3OnOffTest,
    ValuesIn([]() -> std::vector<std::tuple<autopas::ContainerOption, autopas::TraversalOption,
                                            autopas::DataLayoutOption, autopas::InteractionTypeOption>> {
      // needed because CellBlock3D (called when building containers) logs always
      autopas::Logger::create();

      std::vector<std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::DataLayoutOption,
                             autopas::InteractionTypeOption>>
          applicableCombinations;

      // container factory
      autopas::ContainerSelector<autopas::ParticleBaseFP64> containerSelector(
          Newton3OnOffTest::getBoxMin(), Newton3OnOffTest::getBoxMax(), Newton3OnOffTest::getCutoff());
      autopas::ContainerSelectorInfo containerInfo(
          Newton3OnOffTest::getCellSizeFactor(), Newton3OnOffTest::getVerletSkin(),
          Newton3OnOffTest::getRebuildFrequency(), Newton3OnOffTest::getClusterSize(),
          autopas::LoadEstimatorOption::none);

      // generate for all containers
      for (auto containerOption : autopas::ContainerOption::getAllOptions()) {
        containerSelector.selectContainer(containerOption, containerInfo);
        autopas::ParticleContainerInterface<ParticleFP64> &container = containerSelector.getCurrentContainer();

        for (auto dataLayoutOption : autopas::DataLayoutOption::getAllOptions()) {
          for (auto interactionType : autopas::InteractionTypeOption::getMostOptions()) {
            for (auto traversalOption : container.getAllTraversals(interactionType)) {
              // this is the functor that will be used in the test.
              bool configOk = false;

              if (interactionType == autopas::InteractionTypeOption::pairwise) {
                MockPairwiseFunctor<ParticleFP64> functor;
                // generate both newton3 versions of the same traversal and check that both are applicable
                configOk = autopas::utils::withStaticCellType<ParticleFP64>(
                    container.getParticleCellTypeEnum(), [&](auto particleCellDummy) {
                      auto traversalSelector = autopas::TraversalSelector<decltype(particleCellDummy)>();
                      auto traversalWithN3 =
                          traversalSelector.template generateTraversal<MockPairwiseFunctor<ParticleFP64>>(
                              traversalOption, functor, container.getTraversalSelectorInfo(), dataLayoutOption,
                              autopas::Newton3Option::enabled);
                      auto traversalWithoutN3 =
                          traversalSelector.template generateTraversal<MockPairwiseFunctor<ParticleFP64>>(
                              traversalOption, functor, container.getTraversalSelectorInfo(), dataLayoutOption,
                              autopas::Newton3Option::disabled);

                      return traversalWithN3->isApplicable() and traversalWithoutN3->isApplicable();
                    });
              } else if (interactionType == autopas::InteractionTypeOption::triwise) {
                MockTriwiseFunctor<ParticleFP64> functor;
                // generate both newton3 versions of the same traversal and check that both are applicable
                configOk = autopas::utils::withStaticCellType<ParticleFP64>(
                    container.getParticleCellTypeEnum(), [&](auto particleCellDummy) {
                      auto traversalSelector = autopas::TraversalSelector<decltype(particleCellDummy)>();
                      auto traversalWithN3 =
                          traversalSelector.template generateTraversal<MockTriwiseFunctor<ParticleFP64>>(
                              traversalOption, functor, container.getTraversalSelectorInfo(), dataLayoutOption,
                              autopas::Newton3Option::enabled);
                      auto traversalWithoutN3 =
                          traversalSelector.template generateTraversal<MockTriwiseFunctor<ParticleFP64>>(
                              traversalOption, functor, container.getTraversalSelectorInfo(), dataLayoutOption,
                              autopas::Newton3Option::disabled);

                      return traversalWithN3->isApplicable() and traversalWithoutN3->isApplicable();
                    });
              }
              if (configOk) {
                applicableCombinations.emplace_back(containerOption, traversalOption, dataLayoutOption,
                                                    interactionType);
              }
            }
          }
        }
      }

      autopas::Logger::unregister();

      return applicableCombinations;
    }()),
    Newton3OnOffTest::PrintToStringParamName());

// Count number of Functor calls with and without newton 3 and compare
template <typename FunctorType>
void Newton3OnOffTest::countFunctorCalls(autopas::ContainerOption containerOption,
                                         autopas::TraversalOption traversalOption, autopas::DataLayoutOption dataLayout,
                                         autopas::InteractionTypeOption interactionType) {
  // TODO: Make test possible for direct sum SoA
  if (containerOption == autopas::ContainerOption::directSum and dataLayout == autopas::DataLayoutOption::soa) {
    return;
  }
  autopas::ContainerSelector<ParticleFP64> containerSelector(getBoxMin(), getBoxMax(), getCutoff());
  const autopas::ContainerSelectorInfo containerInfo(getCellSizeFactor(), getVerletSkin(), getRebuildFrequency(),
                                                     getClusterSize(), autopas::LoadEstimatorOption::none);

  containerSelector.selectContainer(containerOption, containerInfo);

  autopas::ParticleContainerInterface<ParticleFP64> &container = containerSelector.getCurrentContainer();

  const Molecule defaultParticle{};
  autopasTools::generators::UniformGenerator::fillWithParticles(container, defaultParticle, container.getBoxMin(),
                                                                container.getBoxMax(), 100);
  // Do not add any halo particles to this test!
  // Given an owned particle p1 and a halo particle p2 the following interactions are necessary:
  // Newton3   : p1 <-> p2
  // No Newton3: p1 <-  p2 (but NOT p1 -> p2)
  // Depending if the container is able to avoid all halo <- owned interactions, the total number of functor calls
  // is not trivially to calculate.

  // Create the functor
  FunctorType mockFunctor{};

  EXPECT_CALL(mockFunctor, isRelevantForTuning()).WillRepeatedly(Return(true));

  if (dataLayout == autopas::DataLayoutOption::soa) {
    // loader and extractor will be called, we don't care how often.
    autopas::utils::withStaticCellType<ParticleFP64>(container.getParticleCellTypeEnum(), [&](auto particleCellDummy) {
      auto *expectation =
          &EXPECT_CALL(mockFunctor, SoALoader(::testing::Matcher<decltype(particleCellDummy) &>(_), _, _, _))
               .Times(testing::AtLeast(1));
      // Verlet based containers resize the SoA before they call SoALoader, so no need for the testing::Invoke here
      if (std::set{autopas::ContainerOption::varVerletListsAsBuild, autopas::ContainerOption::pairwiseVerletLists,
                   autopas::ContainerOption::verletListsCells}
              .count(container.getContainerType()) == 0) {
        expectation->WillRepeatedly(
            testing::WithArgs<0, 1>(testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.size()); })));
      }

      EXPECT_CALL(mockFunctor, SoAExtractor(::testing::Matcher<decltype(particleCellDummy) &>(_), _, _))
          .Times(testing::AtLeast(1));
    });
  }

  // "SC" = single cell
  const auto [callsNewton3SC, callsNewton3Pair, callsNewton3Triple] =
      eval(dataLayout, /*useNewton3*/ true, container, traversalOption, mockFunctor);
  const auto [callsNonNewton3SC, callsNonNewton3Pair, callsNonNewton3Triple] =
      eval(dataLayout, /*useNewton3*/ false, container, traversalOption, mockFunctor);

  EXPECT_GT(callsNewton3Triple + callsNewton3Pair + callsNewton3SC, 0)
      << "No interactions were found WITH Newton3. Are the generated particles too far apart?";
  EXPECT_GT(callsNonNewton3Triple + callsNonNewton3Pair + callsNonNewton3SC, 0)
      << "No interactions were found WITHOUT Newton3. Are the generated particles too far apart?";

  if (dataLayout == autopas::DataLayoutOption::soa) {
    // within one cell no N3 optimization
    EXPECT_EQ(callsNewton3SC, callsNonNewton3SC)
        << "Mismatch in number of interactions within cells for container option: " << containerOption;
    // two times the calls for no N3, once for each "owning" cell
    EXPECT_EQ(callsNewton3Pair * 2, callsNonNewton3Pair)
        << "Mismatch in number of interactions between cell pairs for container option: " << containerOption;
    // three times the calls for no N3, once for each "owning" cell
    EXPECT_EQ(callsNewton3Triple * 3, callsNonNewton3Triple)
        << "Mismatch in number of interactions between cell triplets for container option: " << containerOption;

  } else {  // aos
    // should be called exactly two times
    EXPECT_EQ(callsNewton3Pair * 2, callsNonNewton3Pair)
        << "Mismatch in number of interactions between cell pairs for container option: " << containerOption;
    // should be called exactly three times
    EXPECT_EQ(callsNewton3Triple * 3, callsNonNewton3Triple) << "for container option: " << containerOption;
  }

  if (HasFailure()) {
    std::cerr << "Failures for Container: " << containerOption.to_string()
              << ", Traversal: " << traversalOption.to_string() << ", Data Layout: " << dataLayout.to_string()
              << std::endl;
  }
}

template <class Container, class Traversal>
void Newton3OnOffTest::iterate(Container &container, Traversal traversal,
                               autopas::InteractionTypeOption interactionType) {
  container.rebuildNeighborLists(traversal.get());
  container.computeInteractions(traversal.get());
}

template <class Functor, class Container>
std::tuple<size_t, size_t, size_t> Newton3OnOffTest::eval(autopas::DataLayoutOption dataLayout, bool useNewton3,
                                                          Container &container,
                                                          autopas::TraversalOption traversalOption,
                                                          Functor &mockFunctor) {
  // set up counters for particle interactions
  std::atomic<unsigned int> callsSC(0ul);
  std::atomic<unsigned int> callsPair(0ul);
  std::atomic<unsigned int> callsTriple(0ul);

  auto getInteractionTypeFromFunctor = [&]() {
    if constexpr (autopas::utils::isPairwiseFunctor<Functor>()) {
      return autopas::InteractionTypeOption::pairwise;
    } else if constexpr (autopas::utils::isTriwiseFunctor<Functor>()) {
      return autopas::InteractionTypeOption::triwise;
    }
  };
  constexpr auto interactionType = getInteractionTypeFromFunctor();

  EXPECT_CALL(mockFunctor, allowsNewton3()).WillRepeatedly(Return(useNewton3));
  EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillRepeatedly(Return(not useNewton3));

  // depending on the layout, set up some expectations on what Functors are called
  switch (dataLayout) {
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
      ADD_FAILURE() << "This test does not support data layout : " << dataLayout.to_string();
    }
  }

  const autopas::Newton3Option n3Option =
      useNewton3 ? autopas::Newton3Option::enabled : autopas::Newton3Option::disabled;

  auto traversalSelectorInfo = container.getTraversalSelectorInfo();
  // simulate iteration
  autopas::utils::withStaticCellType<ParticleFP64>(container.getParticleCellTypeEnum(), [&](auto particleCellDummy) {
    iterate(container,
            autopas::TraversalSelector<decltype(particleCellDummy)>::template generateTraversal<Functor>(
                traversalOption, mockFunctor, traversalSelectorInfo, dataLayout, n3Option),
            interactionType);
  });

  return std::make_tuple(callsSC.load(), callsPair.load(), callsTriple.load());
}
