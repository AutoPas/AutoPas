/**
 * @file ContainerIteratorTest.cpp
 * @author F.Gratl
 * @date 13.01.23
 */

#include "ContainerIteratorTest.h"

#include "IteratorTestHelper.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "testingHelpers/EmptyPairwiseFunctor.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::computeInteractions(EmptyPairwiseFunctor<Molecule> *);

using ::testing::_;

template <typename AutoPasT>
auto ContainerIteratorTestBase::defaultInit(AutoPasT &autoPas, const autopas::ContainerOption &containerOption,
                                            double cellSizeFactor) {
  using namespace autopas::utils::ArrayMath::literals;
  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1);
  autoPas.setVerletSkin(0.4);
  autoPas.setVerletRebuildFrequency(2);
  autoPas.setNumSamples(2);
  autoPas.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption});
  autoPas.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(
      containerOption, autopas::InteractionTypeOption::pairwise));
  autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>({cellSizeFactor})));

  autoPas.init();

  auto haloBoxMin = autoPas.getBoxMin() - (autoPas.getVerletSkin() + autoPas.getCutoff());
  auto haloBoxMax = autoPas.getBoxMax() + (autoPas.getVerletSkin() + autoPas.getCutoff());

  return std::make_tuple(haloBoxMin, haloBoxMax);
}

template <bool constIter, class AutoPasT, class F>
auto ContainerIteratorTestBase::deleteParticles(AutoPasT &autopas, F predicate, bool useRegionIterator,
                                                const autopas::IteratorBehavior &behavior) {
  if constexpr (not constIter) {
    IteratorTestHelper::provideIterator<false>(autopas, behavior, useRegionIterator, [&](auto &autopas, auto getIter) {
      AUTOPAS_OPENMP(parallel) {
        for (auto iter = getIter(); iter.isValid(); ++iter) {
          const auto leID = iter->getID();
          if (predicate(leID)) {
            autopas.deleteParticle(iter);
          }
        }
      }
    });
  } else {
    GTEST_FAIL() << "Calling deleteParticles with a const iterator! This indicates that the test is ill defined!";
  }
}

/**
 * This Test applies an iterator on the whole domain and expects to find that it is empty.
 * Exact Behavior might vary depending on test parameters but the general flow is:
 * - Create an AutoPas object with a specified container.
 * - Apply an iterator and confirm that it finds no particles.
 */
TEST_P(ContainerIteratorTest, emptyContainer) {
  const auto [containerOption, cellSizeFactor, useRegionIterator, useConstIterator, priorForceCalc, behavior] =
      GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyPairwiseFunctor<Molecule> eFunctor;
    autoPas.computeInteractions(&eFunctor);
  }

  // actual test
  IteratorTestHelper::provideIterator(
      useConstIterator, autoPas, behavior, useRegionIterator,
      [&](const auto &autopas, auto &iter) { IteratorTestHelper::findParticles(autoPas, iter, {}); });
}

/**
 * This Test applies an iterator on the whole domain and expects to find all particles this way.
 * Exact Behavior might vary depending on test parameters but the general flow is:
 * - Create an AutoPas object with a specified container.
 * - Place particles in a grid inside the domain.
 * - Find the particles with iterators and compare their IDs with expectations.
 */
TEST_P(ContainerIteratorTest, findAllParticlesInsideDomain) {
  const auto [containerOption, cellSizeFactor, useRegionIterator, useConstIterator, priorForceCalc, behavior] =
      GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);
  auto expectedIDs = IteratorTestHelper::fillContainerWithGrid(autoPas, 3);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyPairwiseFunctor<Molecule> eFunctor;
    autoPas.computeInteractions(&eFunctor);
  }

  // set up expectations
  switch (behavior) {
    case autopas::IteratorBehavior::ownedOrHalo:
      [[fallthrough]];
    case autopas::IteratorBehavior::owned: {
      // expectations already correct
      break;
    }
    case autopas::IteratorBehavior::halo: {
      // no particles in the halo -> expect nothing
      expectedIDs = {};
      break;
    }
    default: {
      GTEST_FAIL() << "IteratorBehavior::" << behavior
                   << "  should not be tested through this test!\n"
                      "Container behavior with dummy particles is not uniform.\n"
                      "forceSequential alone makes no sense.";
    }
  }

  // actual test
  IteratorTestHelper::provideIterator(
      useConstIterator, autoPas, behavior, useRegionIterator,
      [&](const auto &autopas, auto &iter) { IteratorTestHelper::findParticles(autoPas, iter, expectedIDs); });
}

/**
 * This Test applies an iterator on the whole domain and expects to find all particles this way.
 * Exact Behavior might vary depending on test parameters but the general flow is:
 * - Create an AutoPas object with a specified container.
 * - Strategically place particles around the boundaries.
 * - Find the particles with iterators and compare their IDs with expectations.
 */
TEST_P(ContainerIteratorTest, findAllParticlesAroundBoundaries) {
  const auto [containerOption, cellSizeFactor, useRegionIterator, useConstIterator, priorForceCalc, behavior] =
      GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);
  auto [particleIDsOwned, particleIDsHalo, _, __] = IteratorTestHelper::fillContainerAroundBoundary(autoPas);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyPairwiseFunctor<Molecule> eFunctor;
    autoPas.computeInteractions(&eFunctor);
  }

  // set up expectations
  std::vector<size_t> expectedIDs;
  switch (behavior) {
    case autopas::IteratorBehavior::owned: {
      expectedIDs = particleIDsOwned;
      break;
    }
    case autopas::IteratorBehavior::halo: {
      expectedIDs = particleIDsHalo;
      break;
    }
    case autopas::IteratorBehavior::ownedOrHalo: {
      expectedIDs = particleIDsOwned;
      expectedIDs.insert(expectedIDs.end(), particleIDsHalo.begin(), particleIDsHalo.end());
      break;
    }
    default: {
      GTEST_FAIL() << "IteratorBehavior::" << behavior
                   << "  should not be tested through this test!\n"
                      "Container behavior with dummy particles is not uniform.\n"
                      "Using forceSequential is not supported.";
    }
  }

  // actual test
  IteratorTestHelper::provideIterator(
      useConstIterator, autoPas, behavior, useRegionIterator,
      [&](const auto &autopas, auto &iter) { IteratorTestHelper::findParticles(autoPas, iter, expectedIDs); });
}

/**
 * This test uses an iterator to delete every particle with an odd ID.
 * Since deletion does not work through const iterators this test is skipped when instantiated with
 * useConstIterator==true.
 */
TEST_P(ContainerIteratorTestNonConst, deleteParticles) {
  const auto [containerOption, cellSizeFactor, useRegionIterator, useConstIterator, priorForceCalc, behavior] =
      GetParam();

  if (useConstIterator) {
    GTEST_FAIL() << "This test is not applicable for const iterators and should not be instantiated!";
  }

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);
  auto expectedIDs = IteratorTestHelper::fillContainerWithGrid(autoPas, 3);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyPairwiseFunctor<Molecule> eFunctor;
    autoPas.computeInteractions(&eFunctor);
  }

  auto isOdd = [](auto id) -> bool { return id % 2 != 0; };

  // set up expectations
  switch (behavior) {
    case autopas::IteratorBehavior::ownedOrHalo:
      [[fallthrough]];
    case autopas::IteratorBehavior::owned: {
      // remove all odd numbers from expectations
      expectedIDs.erase(std::remove_if(expectedIDs.begin(), expectedIDs.end(), [&](auto id) { return isOdd(id); }),
                        expectedIDs.end());
      break;
    }
    case autopas::IteratorBehavior::halo: {
      // nothing should be deleted so expect everything.
      break;
    }
    default: {
      GTEST_FAIL() << "IteratorBehavior::" << behavior
                   << "  should not be tested through this test!\n"
                      "Container behavior with dummy particles is not uniform.\n"
                      "Using forceSequential is not supported.";
    }
  }

  // delete all particles with odd ids
  if (useConstIterator) {
    GTEST_FAIL() << "Calling deleteParticles with a const iterator! This indicates that the test is ill defined!";
  } else {
    deleteParticles<false>(autoPas, isOdd, useRegionIterator, behavior);
  }

  // now use again an iterator to confirm only the expected ones are still there
  IteratorTestHelper::provideIterator(
      useConstIterator, autoPas, autopas::IteratorBehavior::ownedOrHalo, useRegionIterator,
      [&](const auto &autopas, auto &iter) { IteratorTestHelper::findParticles(autoPas, iter, expectedIDs); });
}

TEST_P(ContainerIteratorTestNonConstOwned, addParticlesWhileIterating) {
  const auto [containerOption, cellSizeFactor, useRegionIterator, useConstIterator, priorForceCalc, behavior] =
      GetParam();

  if (useConstIterator or behavior != autopas::IteratorBehavior::owned) {
    GTEST_FAIL() << "This test is not applicable for const iterators or iterator behaviors other than owned and should "
                    "not be instantiated!";
  }

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);
  const auto expectedInitialIDs = IteratorTestHelper::fillContainerWithGrid(autoPas, 3);
  const auto highestInitialId = *std::max_element(expectedInitialIDs.begin(), expectedInitialIDs.end());
  // offset between initial and new IDs
  const auto idOffset = static_cast<size_t>(std::pow(10, std::ceil(std::log10(highestInitialId)) + 1));
  const auto expectedNewIDs = [&]() {
    std::vector<size_t> newVec = expectedInitialIDs;
    for (auto &v : newVec) {
      v += idOffset;
    }
    return newVec;
  }();

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyPairwiseFunctor<Molecule> eFunctor;
    autoPas.computeInteractions(&eFunctor);
  }

  const auto cellSize = (autoPas.getCutoff() + autoPas.getVerletSkin()) * cellSizeFactor;
  const auto boxLength = autoPas.getBoxMax();

  std::vector<size_t> foundInitialParticles{};
  std::vector<size_t> foundNewParticles{};
  foundInitialParticles.reserve(expectedInitialIDs.size());
  foundNewParticles.reserve(expectedNewIDs.size());
  // now use again an iterator to confirm only the expected ones are still there
  IteratorTestHelper::provideIterator</* const iter */ false>(
      autoPas, behavior, useRegionIterator, [&](auto &autopas, auto &getIter) {
        for (auto iterator = getIter(); iterator.isValid(); ++iterator) {
          if (iterator->getID() > highestInitialId) {
            foundNewParticles.push_back(iterator->getID());
            continue;
          }
          foundInitialParticles.push_back(iterator->getID());
          auto newPosition = iterator->getR();
          newPosition[0] += cellSize;
          if (newPosition[0] > boxLength[0]) {
            newPosition[0] -= boxLength[0];
          }
          const Molecule newParticle{newPosition, {}, iterator->getID() + idOffset};
          autopas.addParticle(newParticle);
        }
      });
  // The octree might visit some particles multiple times
  const std::set<size_t> foundInitialParticlesUnique{foundInitialParticles.begin(), foundInitialParticles.end()};
  EXPECT_THAT(foundInitialParticlesUnique, ::testing::UnorderedElementsAreArray(expectedInitialIDs));
  // We don't support finding all new particles while iterating.
  //  EXPECT_THAT(foundNewParticles, ::testing::UnorderedElementsAreArray(expectedNewIDs));
}

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::Values;
using ::testing::ValuesIn;

static inline auto getTestableContainerOptions() { return autopas::ContainerOption::getAllOptions(); }

INSTANTIATE_TEST_SUITE_P(Generated, ContainerIteratorTest,
                         Combine(ValuesIn(getTestableContainerOptions()), /*cell size factor*/ Values(0.5, 1., 1.5),
                                 /*use region iter*/ Values(true, false),
                                 /*use const*/ Values(true, false), /*prior force calc*/ Values(true, false),
                                 ValuesIn(autopas::IteratorBehavior::getMostOptions())),
                         ContainerIteratorTestBase::PrintToStringParamName());

INSTANTIATE_TEST_SUITE_P(Generated, ContainerIteratorTestNonConst,
                         Combine(ValuesIn(getTestableContainerOptions()), /*cell size factor*/ Values(0.5, 1., 1.5),
                                 /*use region iter*/ Values(true, false),
                                 /*use const*/ Values(false), /*prior force calc*/ Values(true, false),
                                 ValuesIn(autopas::IteratorBehavior::getMostOptions())),
                         ContainerIteratorTestBase::PrintToStringParamName());

INSTANTIATE_TEST_SUITE_P(Generated, ContainerIteratorTestNonConstOwned,
                         Combine(ValuesIn(getTestableContainerOptions()), /*cell size factor*/ Values(0.5, 1., 1.5),
                                 /*use region iter*/ Values(true, false),
                                 /*use const*/ Values(false), /*prior force calc*/ Values(true, false),
                                 Values(autopas::IteratorBehavior::owned)),
                         ContainerIteratorTestBase::PrintToStringParamName());
