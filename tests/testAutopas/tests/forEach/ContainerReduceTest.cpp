/**
 * @file ContainerReduceTest.cpp
 * @author lgaertner
 * @date 25.08.2021
 */
#include "ContainerReduceTest.h"

#include "ForEachTestHelper.h"
#include "autopas/AutoPasDecl.h"
#include "testingHelpers/EmptyFunctor.h"

extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::iteratePairwise(EmptyFunctor<Molecule> *);

template <typename AutoPasT>
auto ContainerReduceTest::defaultInit(AutoPasT &autoPas, autopas::ContainerOption &containerOption,
                                      double cellSizeFactor) {
  using namespace autopas::utils::ArrayMath::literals;

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

  auto haloBoxMin = autoPas.getBoxMin() - (autoPas.getVerletSkin() + autoPas.getCutoff());
  auto haloBoxMax = autoPas.getBoxMax() + (autoPas.getVerletSkin() + autoPas.getCutoff());

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
 * Test reduceInRegionSequential for every possible combination of [containerOption, cellSizeFactor, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices and
 * checking if the sum of both index lists equal.
 */
TEST_P(ContainerReduceTest, testReduceInRegion) {
  using namespace autopas::utils::ArrayMath::literals;

  auto [containerOption, cellSizeFactor, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  auto domainLength = autoPas.getBoxMax() - autoPas.getBoxMin();
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = domainLength * 0.3;
  std::array<double, 3> searchBoxMin = autoPas.getBoxMin() - searchBoxLengthHalf;
  std::array<double, 3> searchBoxMax = autoPas.getBoxMin() + searchBoxLengthHalf;

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
  std::vector<size_t> particleIDsFound;
  size_t reductionValue = 0ul;

  autoPas.reduceInRegion(
      [&](auto &p, size_t &rv) {
        auto id = p.getID();
        rv += id;
        particleIDsFound.push_back(id);
      },
      reductionValue, searchBoxMin, searchBoxMax, behavior);

  // check that everything was found
  EXPECT_THAT(particleIDsFound, ::testing::UnorderedElementsAreArray(expectedIDs));

  size_t expectedReductionValue = std::accumulate(expectedIDs.begin(), expectedIDs.end(), 0ul);
  EXPECT_EQ(reductionValue, expectedReductionValue);
}

/**
 * Test reduceSequential for every possible combination of [containerOption, cellSizeFactor, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices and
 * checking if the sum of both index lists equal.
 */
TEST_P(ContainerReduceTest, testReduce) {
  using namespace autopas::utils::ArrayMath::literals;

  auto [containerOption, cellSizeFactor, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  auto domainLength = autoPas.getBoxMax() - autoPas.getBoxMin();
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = domainLength * 0.3;
  std::array<double, 3> searchBoxMin = autoPas.getBoxMin() - searchBoxLengthHalf;
  std::array<double, 3> searchBoxMax = autoPas.getBoxMin() + searchBoxLengthHalf;

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
  std::vector<size_t> particleIDsFound;
  size_t reductionValue = 0ul;

  autoPas.reduce(
      [&](auto &p, size_t &rv) {
        auto id = p.getID();
        rv += id;
        particleIDsFound.push_back(id);
      },
      reductionValue, behavior);

  // check that everything was found
  EXPECT_THAT(particleIDsFound, ::testing::UnorderedElementsAreArray(expectedIDs));

  size_t expectedReductionValue = std::accumulate(expectedIDs.begin(), expectedIDs.end(), 0ul);
  EXPECT_EQ(reductionValue, expectedReductionValue);
}

/**
 * Test reduceInRegionParallel for every possible combination of [containerOption, cellSizeFactor, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices and
 * checking if the sum of both index lists equal.
 */
TEST_P(ContainerReduceTest, testReduceInRegionParallel) {
  using namespace autopas::utils::ArrayMath::literals;

  auto [containerOption, cellSizeFactor, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  auto domainLength = autoPas.getBoxMax() - autoPas.getBoxMin();
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = domainLength * 0.3;
  std::array<double, 3> searchBoxMin = autoPas.getBoxMin() - searchBoxLengthHalf;
  std::array<double, 3> searchBoxMax = autoPas.getBoxMin() + searchBoxLengthHalf;

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
  std::vector<size_t> particleIDsFound;
  size_t reductionValue = 0ul;

  autoPas.reduceInRegionParallel(
      [&](auto &p, size_t &rv) {
        auto id = p.getID();
        rv += id;
        particleIDsFound.push_back(id);
      },
      reductionValue, searchBoxMin, searchBoxMax, behavior);

  // check that everything was found
  EXPECT_THAT(particleIDsFound, ::testing::UnorderedElementsAreArray(expectedIDs));

  size_t expectedReductionValue = std::accumulate(expectedIDs.begin(), expectedIDs.end(), 0ul);
  EXPECT_EQ(reductionValue, expectedReductionValue);
}

/**
 * Test reduceParallel for every possible combination of [containerOption, cellSizeFactor, useConstIterator,
 * priorForceCalc, behavior] by adding found particles indices to list and comparing with 'handmade' expectedIndices and
 * checking if the sum of both index lists equal.
 */
TEST_P(ContainerReduceTest, testReduceParallel) {
  using namespace autopas::utils::ArrayMath::literals;

  auto [containerOption, cellSizeFactor, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  auto domainLength = autoPas.getBoxMax() - autoPas.getBoxMin();
  // draw a box around the lower corner of the domain
  auto searchBoxLengthHalf = domainLength * 0.3;
  std::array<double, 3> searchBoxMin = autoPas.getBoxMin() - searchBoxLengthHalf;
  std::array<double, 3> searchBoxMax = autoPas.getBoxMin() + searchBoxLengthHalf;

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
  std::vector<size_t> particleIDsFound;
  size_t reductionValue = 0ul;

  autoPas.reduceParallel(
      [&](auto &p, size_t &rv) {
        auto id = p.getID();
        rv += id;
        particleIDsFound.push_back(id);
      },
      reductionValue, behavior);

  // check that everything was found
  EXPECT_THAT(particleIDsFound, ::testing::UnorderedElementsAreArray(expectedIDs));

  size_t expectedReductionValue = std::accumulate(expectedIDs.begin(), expectedIDs.end(), 0ul);
  EXPECT_EQ(reductionValue, expectedReductionValue);
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

INSTANTIATE_TEST_SUITE_P(Generated, ContainerReduceTest,
                         Combine(ValuesIn(getTestableContainerOptions()), /*cell size factor*/ Values(0.5, 1., 1.5),
                                 /*use const*/ Values(true, false), /*prior force calc*/ Values(true, false),
                                 ValuesIn(autopas::IteratorBehavior::getMostOptions())),
                         ContainerReduceTest::PrintToStringParamName());
