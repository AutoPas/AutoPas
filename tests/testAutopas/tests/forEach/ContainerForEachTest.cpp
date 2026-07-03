/**
 * @file ContainerForEachTest.cpp
 * @author lgaertner
 * @date 25.08.2021
 */
#include "ContainerForEachTest.h"

#include "ForEachTestHelper.h"
#include "autopas/AutoPasDecl.h"
#include "testingHelpers/EmptyPairwiseFunctor.h"

extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::computeInteractions(EmptyPairwiseFunctor<Molecule> *);

template <typename AutoPas_T>
auto ContainerForEachTest::defaultInit(AutoPas_T &autoPas, const ContainerConfiguration &containerConfig) {
  using namespace autopas::utils::ArrayMath::literals;

  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1);
  autoPas.setVerletSkin(0.2);
  autoPas.setVerletRebuildFrequency(2);
  autoPas.setNumSamples(2);
  autoPas.setAllowedContainers(std::set<autopas::ContainerOption>{containerConfig.container});
  // Arbitrarily allow all valid traversals, so AutoPas can be initialized.
  autoPas.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(
      containerConfig.container, autopas::InteractionTypeOption::pairwise));
  autoPas.setCellSizeFactor(containerConfig.cellSizeFactor);

  autoPas.init();

  auto haloBoxMin = autoPas.getBoxMin() - (autoPas.getVerletSkin() + autoPas.getCutoff());
  auto haloBoxMax = autoPas.getBoxMax() + (autoPas.getVerletSkin() + autoPas.getCutoff());

  return std::make_tuple(haloBoxMin, haloBoxMax);
}

/**
 * Returns the expected particle IDs for a given iterator behavior.
 *
 * @param behavior Iterator behavior that determines which IDs are selected. For `owned`, only the owned IDs are
 * returned. For `halo`, only the halo IDs are returned. For `ownedOrHalo`, the owned IDs are returned first, followed
 * by the halo IDs. Any other behavior yields an empty result.
 * @param owned IDs of particles owned by the current rank / container.
 * @param halo IDs of halo particles.
 * @return IDs expected to be visited for the given behavior.
 */
std::vector<size_t> getExpectedIds(autopas::IteratorBehavior behavior, const std::vector<size_t> &owned,
                                   const std::vector<size_t> &halo) {
  std::vector<size_t> expectedIDs;
  switch (behavior) {
    case autopas::IteratorBehavior::owned: {
      expectedIDs = owned;
      break;
    }
    case autopas::IteratorBehavior::halo: {
      expectedIDs = halo;
      break;
    }
    case autopas::IteratorBehavior::ownedOrHalo: {
      expectedIDs = owned;
      expectedIDs.insert(expectedIDs.end(), halo.begin(), halo.end());
      break;
    }
    default: {
      break;
    }
  }
  return expectedIDs;
}

/**
 * Test forEachInRegionSequential for every possible combination of [containerConfig, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices.
 */
TEST_P(ContainerForEachTest, testForEachInRegionSequential) {
  using namespace autopas::utils::ArrayMath::literals;

  auto [containerConfig, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerConfig);

  auto domainLength = autoPas.getBoxMax() - autoPas.getBoxMin();
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = domainLength * 0.3;
  std::array<double, 3> searchBoxMin = autoPas.getBoxMin() - searchBoxLengthHalf;
  std::array<double, 3> searchBoxMax = autoPas.getBoxMin() + searchBoxLengthHalf;

  auto [particleIDsOwned, particleIDsHalo, particleIDsInBoxOwned, particleIDsInBoxHalo] =
      ForEachTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyPairwiseFunctor<Molecule> eFunctor;
    autoPas.computeInteractions(&eFunctor);
  }

  if (behavior & autopas::IteratorBehavior::dummy) {
    GTEST_FAIL() << "IteratorBehavior::" << behavior
                 << "  should not be tested through this test!\n"
                    "Container behavior with dummy particles is not uniform.\n"
                    "Using forceSequential is not supported.";
  }

  std::vector<size_t> expectedIDs = getExpectedIds(behavior, particleIDsInBoxOwned, particleIDsInBoxHalo);

  // sanity check: there should be particles in the expected region
  ASSERT_THAT(expectedIDs, ::testing::Not(::testing::IsEmpty()));

  // actual test
  auto bh = behavior;  // necessary for compiler, behavior not detected as variable
  auto forEachInRegionLambda = [&, bh](auto lambda) {
    autoPas.forEachInRegion(lambda, searchBoxMin, searchBoxMax, bh);
  };
  ForEachTestHelper::forEachParticleTest(forEachInRegionLambda, expectedIDs);
}

/**
 * Test forEachSequential for every possible combination of [containerConfig, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices.
 */
TEST_P(ContainerForEachTest, testForEachSequential) {
  using namespace autopas::utils::ArrayMath::literals;

  auto [containerConfig, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerConfig);

  auto domainLength = autoPas.getBoxMax() - autoPas.getBoxMin();
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = domainLength * 0.3;
  std::array<double, 3> searchBoxMin = autoPas.getBoxMin() - searchBoxLengthHalf;
  std::array<double, 3> searchBoxMax = autoPas.getBoxMin() + searchBoxLengthHalf;

  auto [particleIDsOwned, particleIDsHalo, particleIDsInBoxOwned, particleIDsInBoxHalo] =
      ForEachTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyPairwiseFunctor<Molecule> eFunctor;
    autoPas.computeInteractions(&eFunctor);
  }

  if (behavior & autopas::IteratorBehavior::dummy) {
    GTEST_FAIL() << "IteratorBehavior::" << behavior
                 << "  should not be tested through this test!\n"
                    "Container behavior with dummy particles is not uniform.\n"
                    "Using forceSequential is not supported.";
  }

  std::vector<size_t> expectedIDs = getExpectedIds(behavior, particleIDsOwned, particleIDsHalo);

  // sanity check: there should be particles in the expected region
  ASSERT_THAT(expectedIDs, ::testing::Not(::testing::IsEmpty()));

  // actual test
  auto bh = behavior;  // necessary for compiler, behavior not detected as variable
  auto forEachLambda = [&, bh](auto lambda) { autoPas.forEach(lambda, bh); };
  ForEachTestHelper::forEachParticleTest(forEachLambda, expectedIDs);
}

/**
 * Test forEachInRegionParallel for every possible combination of [containerConfig, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices.
 */
TEST_P(ContainerForEachTest, testForEachInRegionParallel) {
  using namespace autopas::utils::ArrayMath::literals;

  auto [containerConfig, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerConfig);

  auto domainLength = autoPas.getBoxMax() - autoPas.getBoxMin();
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = domainLength * 0.3;
  std::array<double, 3> searchBoxMin = autoPas.getBoxMin() - searchBoxLengthHalf;
  std::array<double, 3> searchBoxMax = autoPas.getBoxMin() + searchBoxLengthHalf;

  auto [particleIDsOwned, particleIDsHalo, particleIDsInBoxOwned, particleIDsInBoxHalo] =
      ForEachTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyPairwiseFunctor<Molecule> eFunctor;
    autoPas.computeInteractions(&eFunctor);
  }

  if (behavior & autopas::IteratorBehavior::dummy) {
    GTEST_FAIL() << "IteratorBehavior::" << behavior
                 << "  should not be tested through this test!\n"
                    "Container behavior with dummy particles is not uniform.\n"
                    "Using forceSequential is not supported.";
  }

  std::vector<size_t> expectedIDs = getExpectedIds(behavior, particleIDsInBoxOwned, particleIDsInBoxHalo);

  // sanity check: there should be particles in the expected region
  ASSERT_THAT(expectedIDs, ::testing::Not(::testing::IsEmpty()));

  // actual test
  auto bh = behavior;  // necessary for compiler, behavior not detected as variable
  auto forEachInRegionLambda = [&, bh](auto lambda) {
    autoPas.forEachInRegionParallel(lambda, searchBoxMin, searchBoxMax, bh);
  };
  ForEachTestHelper::forEachParticleTest(forEachInRegionLambda, expectedIDs);
}

/**
 * Test forEachParallel for every possible combination of [containerConfig, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices.
 */
TEST_P(ContainerForEachTest, testForEachParallel) {
  using namespace autopas::utils::ArrayMath::literals;

  auto [containerConfig, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerConfig);

  auto domainLength = autoPas.getBoxMax() - autoPas.getBoxMin();
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = domainLength * 0.3;
  std::array<double, 3> searchBoxMin = autoPas.getBoxMin() - searchBoxLengthHalf;
  std::array<double, 3> searchBoxMax = autoPas.getBoxMin() + searchBoxLengthHalf;

  auto [particleIDsOwned, particleIDsHalo, particleIDsInBoxOwned, particleIDsInBoxHalo] =
      ForEachTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyPairwiseFunctor<Molecule> eFunctor;
    autoPas.computeInteractions(&eFunctor);
  }

  if (behavior & autopas::IteratorBehavior::dummy) {
    GTEST_FAIL() << "IteratorBehavior::" << behavior
                 << "  should not be tested through this test!\n"
                    "Container behavior with dummy particles is not uniform.\n"
                    "Using forceSequential is not supported.";
  }

  std::vector<size_t> expectedIDs = getExpectedIds(behavior, particleIDsOwned, particleIDsHalo);

  // sanity check: there should be particles in the expected region
  ASSERT_THAT(expectedIDs, ::testing::Not(::testing::IsEmpty()));

  // actual test
  auto bh = behavior;  // necessary for compiler, behavior not detected as variable
  auto forEachLambda = [&, bh](auto lambda) { autoPas.forEachParallel(lambda, bh); };
  ForEachTestHelper::forEachParticleTest(forEachLambda, expectedIDs);
}

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::Values;
using ::testing::ValuesIn;

INSTANTIATE_TEST_SUITE_P(Generated, ContainerForEachTest,
                         Combine(ValuesIn(generateAllValidContainerConfigurations()),
                                 /*use const*/ Values(true, false), /*prior force calc*/ Values(true, false),
                                 ValuesIn(autopas::IteratorBehavior::getMostOptions())),
                         ContainerForEachTest::PrintToStringParamName());
