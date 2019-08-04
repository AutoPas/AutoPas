/**
 * @file Newton3OnOffTest.cpp
 * @author seckler
 * @date 18.04.18
 */

#include "Newton3OnOffTest.h"
#include "autopas/utils/Logger.h"

using ::testing::_;  // anything is ok
using ::testing::Combine;
using ::testing::Return;
using ::testing::ValuesIn;

// Parse combination strings and call actual test function
TEST_P(Newton3OnOffTest, countFunctorCallsTest) {
  auto contTravStr = std::get<0>(GetParam());
  auto dataLayoutStr = std::get<1>(GetParam());

  transform(contTravStr.begin(), contTravStr.end(), contTravStr.begin(), ::tolower);
  transform(dataLayoutStr.begin(), dataLayoutStr.end(), dataLayoutStr.begin(), ::tolower);
  auto containerOption = autopas::utils::StringUtils::parseContainerOptions(contTravStr).begin().operator*();
  auto traversalOption = autopas::utils::StringUtils::parseTraversalOptions(contTravStr).begin().operator*();
  auto dataLayoutOption = autopas::utils::StringUtils::parseDataLayout(dataLayoutStr).begin().operator*();
  countFunctorCalls(containerOption, traversalOption, dataLayoutOption);
}

// Generate Unittests for all Container / Traversal / Datalayout combinations
// Use strings such that generated test names are readable
INSTANTIATE_TEST_SUITE_P(
    Generated, Newton3OnOffTest,
    Combine(ValuesIn([]() -> std::vector<std::string> {
              // needed because CellBlock3D (called when building containers) logs always
              autopas::Logger::create();

              std::vector<std::string> ret;

              for (auto containerOption : autopas::allContainerOptions) {
                // @TODO: let verlet lists support Newton 3
                if (containerOption == autopas::ContainerOption::verletLists ||
                    containerOption == autopas::ContainerOption::verletListsCells ||
                    containerOption == autopas::ContainerOption::verletClusterLists ||
                    containerOption == autopas::ContainerOption::varVerletListsAsBuild ||
                    containerOption == autopas::ContainerOption::adaptiveLinkedCells) {
                  continue;
                }

                autopas::ContainerSelector<Particle, FPCell> containerSelector({0, 0, 0}, {10, 10, 10}, 1);
                autopas::ContainerSelectorInfo containerInfo(1, 0);

                containerSelector.selectContainer(containerOption, containerInfo);

                auto container = containerSelector.getCurrentContainer();

                for (auto traversalOption : container->getAllTraversals()) {
                  if (traversalOption == autopas::TraversalOption::c01 ||
                      traversalOption ==
                          autopas::TraversalOption::c01CombinedSoA /*and autopas::autopas_get_max_threads() > 1*/) {
                    continue;
                  }
                  if (traversalOption == autopas::TraversalOption::c01Cuda) {
                    // Traversal provides no AoS and SoA Traversal
                    continue;
                  }
                  std::stringstream name;

                  name << autopas::utils::StringUtils::to_string(containerOption);
                  name << "+";
                  name << autopas::utils::StringUtils::to_string(traversalOption);

                  ret.push_back(name.str());
                }
              }

              autopas::Logger::unregister();

              return ret;
            }()),
            ValuesIn([]() {
              std::vector<std::string> ret;
              std::transform(
                  autopas::allDataLayoutOptions.begin(), autopas::allDataLayoutOptions.end(), std::back_inserter(ret),
                  [](autopas::DataLayoutOption d) -> std::string { return autopas::utils::StringUtils::to_string(d); });
              return ret;
            }())));

// Count number of Functor calls with and without newton 3 and compare
void Newton3OnOffTest::countFunctorCalls(autopas::ContainerOption containerOption,
                                         autopas::TraversalOption traversalOption,
                                         autopas::DataLayoutOption dataLayout) {
  if (traversalOption == autopas::TraversalOption::c04SoA and dataLayout == autopas::DataLayoutOption::aos) {
    return;
  }

  autopas::ContainerSelector<Particle, FPCell> containerSelector(getBoxMin(), getBoxMax(), getCutoff());
  autopas::ContainerSelectorInfo containerInfo(getCellSizeFactor(), getVerletSkin());

  containerSelector.selectContainer(containerOption, containerInfo);

  auto container = containerSelector.getCurrentContainer();

  autopas::MoleculeLJ defaultParticle;
  RandomGenerator::fillWithParticles(*container, defaultParticle, 100);
  RandomGenerator::fillWithHaloParticles(*container, defaultParticle, container->getCutoff(), 10);

  EXPECT_CALL(mockFunctor, isRelevantForTuning()).WillRepeatedly(Return(true));

  if (dataLayout == autopas::DataLayoutOption::soa) {
    // loader and extractor will be called, we don't care how often.
    EXPECT_CALL(mockFunctor, SoALoader(_, _))
        .Times(testing::AtLeast(1))
        .WillRepeatedly(testing::WithArgs<0, 1>(
            testing::Invoke([](auto &cell, auto &buf) { buf.resizeArrays(cell.numParticles()); })));
    EXPECT_CALL(mockFunctor, SoAExtractor(_, _)).Times(testing::AtLeast(1));
  }

  const auto [callsNewton3SC, callsNewton3Pair] = eval<true>(dataLayout, container, traversalOption);
  const auto [callsNonNewton3SC, callsNonNewton3Pair] = eval<false>(dataLayout, container, traversalOption);

  if (autopas::DataLayoutOption::soa) {
    // within one cell no N3 optimization
    EXPECT_EQ(callsNewton3SC, callsNonNewton3SC) << "for containeroption: " << containerOption;
  }
  // should be called exactly two times
  EXPECT_EQ(callsNewton3Pair * 2, callsNonNewton3Pair) << "for containeroption: " << containerOption;

  if (::testing::Test::HasFailure()) {
    std::cerr << "Failures for Container: " << autopas::utils::StringUtils::to_string(containerOption)
              << ", Traversal: " << autopas::utils::StringUtils::to_string(traversalOption)
              << ", Data Layout: " << autopas::utils::StringUtils::to_string(dataLayout) << std::endl;
  }
}

template <class ParticleFunctor, class Container, class Traversal>
void Newton3OnOffTest::iterate(Container container, Traversal traversal, autopas::DataLayoutOption dataLayout,
                               autopas::Newton3Option newton3, ParticleFunctor *f) {
  container->iteratePairwise(traversal.get());
}

template <bool useNewton3, class Container, class Traversal>
std::pair<size_t, size_t> Newton3OnOffTest::eval(autopas::DataLayoutOption dataLayout, Container &container,
                                                 Traversal traversalOption) {
  std::atomic<unsigned int> callsSC(0ul);
  std::atomic<unsigned int> callsPair(0ul);
  EXPECT_CALL(mockFunctor, allowsNewton3()).WillRepeatedly(Return(useNewton3));
  EXPECT_CALL(mockFunctor, allowsNonNewton3()).WillRepeatedly(Return(not useNewton3));

  auto traversalSelectorInfo = container->getTraversalSelectorInfo();
  const autopas::Newton3Option n3option =
      (useNewton3) ? autopas::Newton3Option::enabled : autopas::Newton3Option::disabled;

  switch (dataLayout) {
    case autopas::DataLayoutOption::soa: {
      // single cell
      EXPECT_CALL(mockFunctor, SoAFunctor(_, useNewton3))
          .Times(testing::AtLeast(1))
          .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsSC++; }));

      // pair of cells
      EXPECT_CALL(mockFunctor, SoAFunctor(_, _, useNewton3))
          .Times(testing::AtLeast(1))
          .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsPair++; }));

      // non useNewton3 variant should not happen
      EXPECT_CALL(mockFunctor, SoAFunctor(_, _, not useNewton3)).Times(0);
      iterate(
          container,
          autopas::TraversalSelector<FPCell>::template generateTraversal<MockFunctor<Particle, FPCell>,
                                                                         autopas::DataLayoutOption::soa, useNewton3>(
              traversalOption, mockFunctor, std::move(traversalSelectorInfo)),
          dataLayout, n3option, &mockFunctor);
      break;
    }
    case autopas::DataLayoutOption::aos: {
      EXPECT_CALL(mockFunctor, AoSFunctor(_, _, useNewton3))
          .Times(testing::AtLeast(1))
          .WillRepeatedly(testing::InvokeWithoutArgs([&]() { callsPair++; }));

      // non useNewton3 variant should not happen
      EXPECT_CALL(mockFunctor, AoSFunctor(_, _, not useNewton3)).Times(0);
      iterate(
          container,
          autopas::TraversalSelector<FPCell>::template generateTraversal<MockFunctor<Particle, FPCell>,
                                                                         autopas::DataLayoutOption::aos, useNewton3>(
              traversalOption, mockFunctor, std::move(traversalSelectorInfo)),
          dataLayout, n3option, &mockFunctor);
      break;
    }
    case autopas::DataLayoutOption::cuda: {
      // supresses Warning
      // TODO implement suitable test
      iterate(
          container,
          autopas::TraversalSelector<FPCell>::template generateTraversal<MockFunctor<Particle, FPCell>,
                                                                         autopas::DataLayoutOption::cuda, useNewton3>(
              traversalOption, mockFunctor, std::move(traversalSelectorInfo)),
          dataLayout, n3option, &mockFunctor);
      break;
    }
    default:
      ADD_FAILURE() << "This test does not support data layout : "
                    << autopas::utils::StringUtils::to_string(dataLayout);
  }

  return std::make_pair(callsSC.load(), callsPair.load());
}