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
#include "autopasTools/generators/RandomGenerator.h"

using ::testing::_;
using ::testing::Combine;
using ::testing::Return;
using ::testing::ValuesIn;

// Parse combination strings and call actual test function
TEST_P(Newton3OnOffTest, countFunctorCallsTest) {
  auto [containerOption, traversalOption, dataLayoutOption] = GetParam();

  countFunctorCalls(containerOption, traversalOption, dataLayoutOption);
}

// Generate Unittests for all Configurations that support both Newton3 modes
INSTANTIATE_TEST_SUITE_P(
    Generated, Newton3OnOffTest,
    ValuesIn(
        []() -> std::vector<std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::DataLayoutOption>> {
          // needed because CellBlock3D (called when building containers) logs always
          autopas::Logger::create();

          std::vector<std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::DataLayoutOption>>
              applicableCombinations;

          // container factory
          autopas::ContainerSelector<Particle> containerSelector({0., 0., 0.}, {10., 10., 10.}, 1.);
          autopas::ContainerSelectorInfo containerInfo(
              Newton3OnOffTest::getCellSizeFactor(), Newton3OnOffTest::getVerletSkinPerTimestep(),
              Newton3OnOffTest::getRebuildFrequency(), Newton3OnOffTest::getClusterSize(),
              autopas::LoadEstimatorOption::none);

          // generate for all containers
          for (auto containerOption : autopas::ContainerOption::getAllOptions()) {
            containerSelector.selectContainer(containerOption, containerInfo);
            autopas::ParticleContainerInterface<Particle> &container = containerSelector.getCurrentContainer();

            for (auto traversalOption : container.getAllTraversals()) {
              for (auto dataLayoutOption : autopas::DataLayoutOption::getAllOptions()) {
                // this is the functor that will be used in the test.
                MockFunctor<Particle> f;
                // generate both newton3 versions of the same traversal and check that both are applicable
                bool configOk = autopas::utils::withStaticCellType<Particle>(
                    container.getParticleCellTypeEnum(), [&](auto particleCellDummy) {
                      auto traversalSelector = autopas::TraversalSelector<decltype(particleCellDummy)>();
                      auto traversalWithN3 =
                          traversalSelector.generateTraversal(traversalOption, f, container.getTraversalSelectorInfo(),
                                                              dataLayoutOption, autopas::Newton3Option::enabled);
                      auto traversalWithoutN3 =
                          traversalSelector.generateTraversal(traversalOption, f, container.getTraversalSelectorInfo(),
                                                              dataLayoutOption, autopas::Newton3Option::disabled);

                      return traversalWithN3->isApplicable() and traversalWithoutN3->isApplicable();
                    });

                if (configOk) {
                  applicableCombinations.emplace_back(containerOption, traversalOption, dataLayoutOption);
                }
              }
            }
          }

          autopas::Logger::unregister();

          return applicableCombinations;
        }()),
    Newton3OnOffTest::PrintToStringParamName());

// Count number of Functor calls with and without newton 3 and compare
void Newton3OnOffTest::countFunctorCalls(autopas::ContainerOption containerOption,
                                         autopas::TraversalOption traversalOption,
                                         autopas::DataLayoutOption dataLayout) {
  autopas::ContainerSelector<Particle> containerSelector(getBoxMin(), getBoxMax(), getCutoff());
  autopas::ContainerSelectorInfo containerInfo(getCellSizeFactor(), getVerletSkinPerTimestep(), getRebuildFrequency(),
                                               getClusterSize(), autopas::LoadEstimatorOption::none);

  containerSelector.selectContainer(containerOption, containerInfo);

  autopas::ParticleContainerInterface<Particle> &container = containerSelector.getCurrentContainer();

  Molecule defaultParticle;
  autopasTools::generators::RandomGenerator::fillWithParticles(container, defaultParticle, container.getBoxMin(),
                                                               container.getBoxMax(), 100);
  // Do not add any halo particles to this test!
  // Given an owned particle p1 and a halo particle p2 the following interactions are necessary:
  // Newton3   : p1 <-> p2
  // No Newton3: p1 <-  p2 (but NOT p1 -> p2)
  // Depending if the container is able to avoid all halo <- owned interactions, the total number of functor calls
  // is not trivially to calculate.

  EXPECT_CALL(mockFunctor, isRelevantForTuning()).WillRepeatedly(Return(true));

  if (dataLayout == autopas::DataLayoutOption::soa) {
    // loader and extractor will be called, we don't care how often.
    autopas::utils::withStaticCellType<Particle>(container.getParticleCellTypeEnum(), [&](auto particleCellDummy) {
      auto *expectation =
          &EXPECT_CALL(mockFunctor, SoALoader(::testing::Matcher<decltype(particleCellDummy) &>(_), _, _, _))
               .Times(testing::AtLeast(1));
      // Verlet based containers resize the SoA before they call SoALoader, so no need for the testing::Invoke here
      if (std::set{autopas::ContainerOption::varVerletListsAsBuild, autopas::ContainerOption::pairwiseVerletLists,
                   autopas::ContainerOption::verletListsCells, autopas::ContainerOption::dynamicVerletListsCells,
                   autopas::ContainerOption::dynamicPairwiseVerletLists}
              .count(container.getContainerType()) == 0) {
        expectation->WillRepeatedly(
            testing::WithArgs<0, 1>(testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.size()); })));
      }

      EXPECT_CALL(mockFunctor, SoAExtractor(::testing::Matcher<decltype(particleCellDummy) &>(_), _, _))
          .Times(testing::AtLeast(1));
    });
  }
  // "SC" = single cell
  const auto [callsNewton3SC, callsNewton3Pair] = eval(dataLayout, /*useNewton3*/ true, container, traversalOption);
  const auto [callsNonNewton3SC, callsNonNewton3Pair] =
      eval(dataLayout, /*useNewton3*/ false, container, traversalOption);

  if (dataLayout == autopas::DataLayoutOption::soa) {
    // within one cell no N3 optimization.
    // Don't check if there is at least one, because verlet style
    // algorithms don't necessarily have a single cell interaction
    EXPECT_EQ(callsNewton3SC, callsNonNewton3SC)
        << "Mismatch in cell self interactions for container option: " << containerOption;
  }

  EXPECT_GT(callsNewton3Pair, 0)
      << "Test generated particles too far apart so there aren't be any functor calls to count";
  // should be called exactly two times
  EXPECT_EQ(callsNewton3Pair * 2, callsNonNewton3Pair)
      << "Mismatch in cell pair interactions for container option: " << containerOption;

  if (::testing::Test::HasFailure()) {
    std::cerr << "Failures for Container: " << containerOption.to_string()
              << ", Traversal: " << traversalOption.to_string() << ", Data Layout: " << dataLayout.to_string()
              << std::endl;
  }
}

template <class Container, class Traversal>
void Newton3OnOffTest::iterate(Container &container, Traversal traversal) {
  container.rebuildNeighborLists(traversal.get());
  container.iteratePairwise(traversal.get());
}

template <class Container, class Traversal>
std::pair<size_t, size_t> Newton3OnOffTest::eval(autopas::DataLayoutOption dataLayout, bool useNewton3,
                                                 Container &container, Traversal traversalOption) {
  std::atomic<unsigned int> callsSC(0ul);
  std::atomic<unsigned int> callsPair(0ul);
  EXPECT_CALL(mockFunctor, allowsNewton3()).WillRepeatedly(Return(useNewton3));
  EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillRepeatedly(Return(not useNewton3));

  auto traversalSelectorInfo = container.getTraversalSelectorInfo();

  // depending on the layout set up expectations on what Functors are called
  switch (dataLayout) {
    case autopas::DataLayoutOption::soa: {
      // some containers actually use different SoA functor calls so expect them instead of the regular ones
      if (container.getContainerType() == autopas::ContainerOption::varVerletListsAsBuild ||
          container.getContainerType() == autopas::ContainerOption::pairwiseVerletLists ||
          container.getContainerType() == autopas::ContainerOption::verletListsCells ||
          container.getContainerType() == autopas::ContainerOption::dynamicVerletListsCells ||
          container.getContainerType() == autopas::ContainerOption::dynamicPairwiseVerletLists) {
        EXPECT_CALL(mockFunctor, SoAFunctorVerlet(_, _, _, useNewton3))
            .Times(testing::AtLeast(1))
            .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsPair++; }));

        // non useNewton3 variant should not happen
        EXPECT_CALL(mockFunctor, SoAFunctorVerlet(_, _, _, not useNewton3)).Times(0);

      } else {
        // single cell
        EXPECT_CALL(mockFunctor, SoAFunctorSingle(_, useNewton3))
            .Times(testing::AtLeast(1))
            .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsSC++; }));

        // pair of cells
        EXPECT_CALL(mockFunctor, SoAFunctorPair(_, _, useNewton3))
            .Times(testing::AtLeast(1))
            .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsPair++; }));

        // non useNewton3 variant should not happen
        EXPECT_CALL(mockFunctor, SoAFunctorPair(_, _, not useNewton3)).Times(0);
      }
      break;
    }
    case autopas::DataLayoutOption::aos: {
      EXPECT_CALL(mockFunctor, AoSFunctor(_, _, useNewton3))
          .Times(testing::AtLeast(1))
          .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsPair++; }));

      // non useNewton3 variant should not happen
      EXPECT_CALL(mockFunctor, AoSFunctor(_, _, not useNewton3)).Times(0);

      break;
    }
    default: {
      ADD_FAILURE() << "This test does not support data layout : " << dataLayout.to_string();
    }
  }

  const autopas::Newton3Option n3Option =
      useNewton3 ? autopas::Newton3Option::enabled : autopas::Newton3Option::disabled;

  // simulate iteration
  autopas::utils::withStaticCellType<Particle>(container.getParticleCellTypeEnum(), [&](auto particleCellDummy) {
    iterate(container,
            autopas::TraversalSelector<decltype(particleCellDummy)>::template generateTraversal<MockFunctor<Particle>>(
                traversalOption, mockFunctor, traversalSelectorInfo, dataLayout, n3Option));
  });

  return std::make_pair(callsSC.load(), callsPair.load());
}
