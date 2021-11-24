/**
 * @file ContainerForEachTest.cpp
 * @author lgaertner
 * @date 25.08.2021
 */
#include "ContainerForEachTest.h"

#include "ForEachTestHelper.h"
#include "autopas/AutoPas.h"
#include "testingHelpers/EmptyFunctor.h"

template <typename AutoPasT>
auto ContainerForEachTest::defaultInit(AutoPasT &autoPas, autopas::ContainerOption &containerOption,
                                       double cellSizeFactor) {
  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1);
  autoPas.setVerletSkin(0.2);
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

template <typename Ids>
std::vector<size_t> getExpectedIds(autopas::IteratorBehavior behavior, Ids owned, Ids halo) {
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
 * Test forEachInRegionSequential for every possible combination of [containerOption, cellSizeFactor, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices.
 */
TEST_P(ContainerForEachTest, testForEachInRegionSequential) {
  auto [containerOption, cellSizeFactor, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  using ::autopas::utils::ArrayMath::add;
  using ::autopas::utils::ArrayMath::mulScalar;
  using ::autopas::utils::ArrayMath::sub;
  auto domainLength = sub(autoPas.getBoxMax(), autoPas.getBoxMin());
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = mulScalar(domainLength, 0.3);
  std::array<double, 3> searchBoxMin = sub(autoPas.getBoxMin(), searchBoxLengthHalf);
  std::array<double, 3> searchBoxMax = add(autoPas.getBoxMin(), searchBoxLengthHalf);

  auto [particleIDsOwned, particleIDsHalo, particleIDsInBoxOwned, particleIDsInBoxHalo] =
      ForEachTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
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
  ForEachTestHelper::forEachParticleTest(forEachInRegionLambda, expectedIDs,
                                         autoPas.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo));
}

/**
 * Test forEachSequential for every possible combination of [containerOption, cellSizeFactor, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices.
 */
TEST_P(ContainerForEachTest, testForEachSequential) {
  auto [containerOption, cellSizeFactor, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  using ::autopas::utils::ArrayMath::add;
  using ::autopas::utils::ArrayMath::mulScalar;
  using ::autopas::utils::ArrayMath::sub;
  auto domainLength = sub(autoPas.getBoxMax(), autoPas.getBoxMin());
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = mulScalar(domainLength, 0.3);
  std::array<double, 3> searchBoxMin = sub(autoPas.getBoxMin(), searchBoxLengthHalf);
  std::array<double, 3> searchBoxMax = add(autoPas.getBoxMin(), searchBoxLengthHalf);

  auto [particleIDsOwned, particleIDsHalo, particleIDsInBoxOwned, particleIDsInBoxHalo] =
      ForEachTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
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
  ForEachTestHelper::forEachParticleTest(forEachLambda, expectedIDs,
                                         autoPas.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo));
}

/**
 * Test forEachInRegionParallel for every possible combination of [containerOption, cellSizeFactor, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices.
 */
TEST_P(ContainerForEachTest, testForEachInRegionParallel) {
  auto [containerOption, cellSizeFactor, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  using ::autopas::utils::ArrayMath::add;
  using ::autopas::utils::ArrayMath::mulScalar;
  using ::autopas::utils::ArrayMath::sub;
  auto domainLength = sub(autoPas.getBoxMax(), autoPas.getBoxMin());
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = mulScalar(domainLength, 0.3);
  std::array<double, 3> searchBoxMin = sub(autoPas.getBoxMin(), searchBoxLengthHalf);
  std::array<double, 3> searchBoxMax = add(autoPas.getBoxMin(), searchBoxLengthHalf);

  auto [particleIDsOwned, particleIDsHalo, particleIDsInBoxOwned, particleIDsInBoxHalo] =
      ForEachTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
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
  ForEachTestHelper::forEachParticleTest(forEachInRegionLambda, expectedIDs,
                                         autoPas.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo));
}

/**
 * Test forEachParallel for every possible combination of [containerOption, cellSizeFactor, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices.
 */
TEST_P(ContainerForEachTest, testForEachParallel) {
  auto [containerOption, cellSizeFactor, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  using ::autopas::utils::ArrayMath::add;
  using ::autopas::utils::ArrayMath::mulScalar;
  using ::autopas::utils::ArrayMath::sub;
  auto domainLength = sub(autoPas.getBoxMax(), autoPas.getBoxMin());
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = mulScalar(domainLength, 0.3);
  std::array<double, 3> searchBoxMin = sub(autoPas.getBoxMin(), searchBoxLengthHalf);
  std::array<double, 3> searchBoxMax = add(autoPas.getBoxMin(), searchBoxLengthHalf);

  auto [particleIDsOwned, particleIDsHalo, particleIDsInBoxOwned, particleIDsInBoxHalo] =
      ForEachTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
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
  ForEachTestHelper::forEachParticleTest(forEachLambda, expectedIDs,
                                         autoPas.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo));
}

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::Values;
using ::testing::ValuesIn;

static inline auto getTestableContainerOptions() {
#ifdef AUTOPAS_CUDA
  return autopas::ContainerOption::getAllOptions();
#else
  auto containerOptions = autopas::ContainerOption::getAllOptions();
  return containerOptions;
#endif
}

INSTANTIATE_TEST_SUITE_P(Generated, ContainerForEachTest,
                         Combine(ValuesIn(getTestableContainerOptions()), /*cell size factor*/ Values(0.5, 1., 1.5),
                                 /*use const*/ Values(true, false), /*prior force calc*/ Values(true, false),
                                 ValuesIn(autopas::IteratorBehavior::getMostOptions())),
                         ContainerForEachTest::PrintToStringParamName());
