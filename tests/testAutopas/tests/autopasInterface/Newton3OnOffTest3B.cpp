/**
 * @file Newton3OnOffTest3B.cpp
 * @author muehlhaeusser
 * @date 26.10.23
 */

#include "Newton3OnOffTest3B.h"

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
TEST_P(Newton3OnOffTest3B, countFunctorCallsTest) {
  auto [containerOption, traversalOption, dataLayoutOption] = GetParam();

  countFunctorCalls(containerOption, traversalOption, dataLayoutOption);
}

// Generate Unittests for all Configurations that support both Newton3 modes
INSTANTIATE_TEST_SUITE_P(
    Generated, Newton3OnOffTest3B,
    ValuesIn(
        []() -> std::vector<std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::DataLayoutOption>> {
          // needed because CellBlock3D (called when building containers) logs always
          autopas::Logger::create();

          std::vector<std::tuple<autopas::ContainerOption, autopas::TraversalOption, autopas::DataLayoutOption>>
              applicableCombinations;

          // container factory
          autopas::ContainerSelector<Particle> containerSelector({0., 0., 0.}, {10., 10., 10.}, 1.);
          autopas::ContainerSelectorInfo containerInfo(
              Newton3OnOffTest3B::getCellSizeFactor(), Newton3OnOffTest3B::getVerletSkinPerTimestep(),
              Newton3OnOffTest3B::getRebuildFrequency(), Newton3OnOffTest3B::getClusterSize(),
              autopas::LoadEstimatorOption::none);

          // generate for all containers
          for (auto containerOption : autopas::ContainerOption::getAllOptions()) {
            containerSelector.selectContainer(containerOption, containerInfo);
            autopas::ParticleContainerInterface<Particle> &container = containerSelector.getCurrentContainer();

            for (auto traversalOption : container.getAllTraversals(autopas::InteractionTypeOption::triwise)) {
              for (auto dataLayoutOption : autopas::DataLayoutOption::getAllOptions()) {
                // this is the functor that will be used in the test.
                MTriwiseFunctor f;
                // generate both newton3 versions of the same traversal and check that both are applicable
                bool configOk = autopas::utils::withStaticCellType<Particle>(
                    container.getParticleCellTypeEnum(), [&](auto particleCellDummy) {
                      auto traversalSelector = autopas::TraversalSelector<decltype(particleCellDummy)>();
                      auto traversalWithN3 =
                          traversalSelector.template generateTraversal<MTriwiseFunctor, autopas::InteractionTypeOption::triwise>(traversalOption, f, container.getTraversalSelectorInfo(),
                                                              dataLayoutOption, autopas::Newton3Option::enabled);
                      auto traversalWithoutN3 =
                          traversalSelector.template generateTraversal<MTriwiseFunctor, autopas::InteractionTypeOption::triwise>(traversalOption, f, container.getTraversalSelectorInfo(),
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
    Newton3OnOffTest3B::PrintToStringParamName());

// Count number of Functor calls with and without newton 3 and compare
void Newton3OnOffTest3B::countFunctorCalls(autopas::ContainerOption containerOption,
                                           autopas::TraversalOption traversalOption,
                                           autopas::DataLayoutOption dataLayout) {
  autopas::ContainerSelector<Particle> containerSelector(getBoxMin(), getBoxMax(), getCutoff());
  autopas::ContainerSelectorInfo containerInfo(getCellSizeFactor(), getVerletSkinPerTimestep(), getRebuildFrequency(),
                                               getClusterSize(), autopas::LoadEstimatorOption::none);

  containerSelector.selectContainer(containerOption, containerInfo);

  autopas::ParticleContainerInterface<Particle> &container = containerSelector.getCurrentContainer();

  Particle defaultParticle;
  autopasTools::generators::RandomGenerator::fillWithParticles(container, defaultParticle, container.getBoxMin(),
                                                               container.getBoxMax(), 80);
  autopasTools::generators::RandomGenerator::fillWithHaloParticles(container, defaultParticle, container.getCutoff(),
                                                                   8);

  EXPECT_CALL(mockFunctor, isRelevantForTuning()).WillRepeatedly(Return(true));

  if (dataLayout == autopas::DataLayoutOption::soa) {
    // loader and extractor will be called, we don't care how often.
    autopas::utils::withStaticCellType<Particle>(container.getParticleCellTypeEnum(), [&](auto particleCellDummy) {
      EXPECT_CALL(mockFunctor, SoALoader(::testing::Matcher<decltype(particleCellDummy) &>(_), _, _, _))
          .Times(testing::AtLeast(1))
          .WillRepeatedly(
              testing::WithArgs<0, 1>(testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.size()); })));
      EXPECT_CALL(mockFunctor, SoAExtractor(::testing::Matcher<decltype(particleCellDummy) &>(_), _, _))
          .Times(testing::AtLeast(1));
    });
  }

  const auto [callsNewton3SC, callsNewton3Pair, callsNewton3Triple] =
      eval<true>(dataLayout, container, traversalOption);
  const auto [callsNonNewton3SC, callsNonNewton3Pair, callsNonNewton3Triple] =
      eval<false>(dataLayout, container, traversalOption);

  if (dataLayout == autopas::DataLayoutOption::soa) {
    // within one cell no N3 optimization
    EXPECT_EQ(callsNewton3SC, callsNonNewton3SC) << "for containeroption: " << containerOption;
    // two times the calls for no N3, once for each "owning" cell
    EXPECT_EQ(callsNewton3Pair * 2, callsNonNewton3Pair) << "for containeroption: " << containerOption;
    // three times the calls for no N3, once for each "owning" cell
    EXPECT_EQ(callsNewton3Triple * 3, callsNonNewton3Triple) << "for containeroption: " << containerOption;
  } else {  // aos
    // should be called exactly three times
    EXPECT_EQ(callsNewton3Triple * 3, callsNonNewton3Triple) << "for containeroption: " << containerOption;
  }

  if (::testing::Test::HasFailure()) {
    std::cerr << "Failures for Container: " << containerOption.to_string()
              << ", Traversal: " << traversalOption.to_string() << ", Data Layout: " << dataLayout.to_string()
              << std::endl;
  }
}

template <class Container, class Traversal>
void Newton3OnOffTest3B::iterate(Container &container, Traversal traversal) {
  auto triwiseTraversal =
      dynamic_cast<autopas::TraversalInterface<autopas::InteractionTypeOption::triwise> *>(traversal.get());
  //  container.rebuildNeighborLists(triwiseTraversal);
  container.iterateTriwise(triwiseTraversal);
}

template <bool useNewton3, class Container, class Traversal>
std::tuple<size_t, size_t, size_t> Newton3OnOffTest3B::eval(autopas::DataLayoutOption dataLayout, Container &container,
                                                            Traversal traversalOption) {
  std::atomic<unsigned int> callsSC(0ul);
  std::atomic<unsigned int> callsPair(0ul);
  std::atomic<unsigned int> callsTriple(0ul);
  EXPECT_CALL(mockFunctor, allowsNewton3()).WillRepeatedly(Return(useNewton3));
  EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillRepeatedly(Return(not useNewton3));

  auto traversalSelectorInfo = container.getTraversalSelectorInfo();

  // depending on the layout set up expectations on what Functors are called
  switch (dataLayout) {
    case autopas::DataLayoutOption::soa: {
      // some containers actually use different SoA functor calls so expect them instead of the regular ones
      if (container.getContainerType() == autopas::ContainerOption::varVerletListsAsBuild ||
          container.getContainerType() == autopas::ContainerOption::pairwiseVerletLists ||
          container.getContainerType() == autopas::ContainerOption::verletListsCells) {
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

        // triplet of cells
        EXPECT_CALL(mockFunctor, SoAFunctorTriple(_, _, _, useNewton3))
            .Times(testing::AtLeast(0))  // Direct sum has only 2 cells
            .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsTriple++; }));

        // non useNewton3 variant should not happen
        EXPECT_CALL(mockFunctor, SoAFunctorPair(_, _, not useNewton3)).Times(0);

        // non useNewton3 variant should not happen
        EXPECT_CALL(mockFunctor, SoAFunctorTriple(_, _, _, not useNewton3)).Times(0);
      }
      break;
    }
    case autopas::DataLayoutOption::aos: {
      EXPECT_CALL(mockFunctor, AoSFunctor(_, _, _, useNewton3))
          .Times(testing::AtLeast(1))
          .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsTriple++; }));

      // non useNewton3 variant should not happen
      EXPECT_CALL(mockFunctor, AoSFunctor(_, _, _, not useNewton3)).Times(0);

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
            autopas::TraversalSelector<decltype(particleCellDummy)>::
                template generateTraversal<MTriwiseFunctor, autopas::InteractionTypeOption::triwise>(traversalOption, mockFunctor, traversalSelectorInfo,
                                                            dataLayout, n3Option));
  });

  return std::make_tuple(callsSC.load(), callsPair.load(), callsTriple.load());
}
