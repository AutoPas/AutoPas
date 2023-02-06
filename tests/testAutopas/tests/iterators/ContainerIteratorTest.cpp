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
#include "autopasTools/generators/RandomGenerator.h"
#include "testingHelpers/EmptyFunctor.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::iteratePairwise(EmptyFunctor<Molecule> *);

using ::testing::_;

template <typename AutoPasT>
auto ContainerIteratorTest::defaultInit(AutoPasT &autoPas, autopas::ContainerOption &containerOption,
                                        double cellSizeFactor) {
  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1);
  autoPas.setVerletSkinPerTimestep(0.2);
  autoPas.setVerletRebuildFrequency(2);
  autoPas.setNumSamples(2);
  autoPas.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption});
  autoPas.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(containerOption));
  autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>({cellSizeFactor})));

  autoPas.init();

  auto haloBoxMin =
      autopas::utils::ArrayMath::subScalar(autoPas.getBoxMin(), autoPas.getVerletSkin() + autoPas.getCutoff());
  auto haloBoxMax =
      autopas::utils::ArrayMath::addScalar(autoPas.getBoxMax(), autoPas.getVerletSkin() + autoPas.getCutoff());

  return std::make_tuple(haloBoxMin, haloBoxMax);
}

template <bool constIter, class AutoPasT, class F>
auto ContainerIteratorTest::deleteParticles(AutoPasT &autopas, F predicate, bool useRegionIterator,
                                            const autopas::IteratorBehavior &behavior) {
  if constexpr (not constIter) {
    IteratorTestHelper::provideIterator<false>(autopas, behavior, useRegionIterator, [&](auto &autopas, auto getIter) {
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
      {
        for (auto iter = getIter(); iter.isValid(); ++iter) {
          if (predicate(iter->getID())) {
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
  auto [containerOption, cellSizeFactor, useRegionIterator, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
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
  auto [containerOption, cellSizeFactor, useRegionIterator, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);
  auto expectedIDs = IteratorTestHelper::fillContainerWithGrid(autoPas, 3);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
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
  auto [containerOption, cellSizeFactor, useRegionIterator, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);
  auto [particleIDsOwned, particleIDsHalo, _, __] = IteratorTestHelper::fillContainerAroundBoundary(autoPas);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
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
TEST_P(ContainerIteratorTest, deleteParticles) {
  auto [containerOption, cellSizeFactor, useRegionIterator, useConstIterator, priorForceCalc, behavior] = GetParam();

  if (useConstIterator) {
    GTEST_SKIP_("Not applicable since deleting with a const iterator is not possible");
  }

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);
  auto expectedIDs = IteratorTestHelper::fillContainerWithGrid(autoPas, 3);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
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
                         ContainerIteratorTest::PrintToStringParamName());
