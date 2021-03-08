/**
 * @file RegionParticleIteratorTest.cpp
 * @author F. Gratl
 * @date 08.03.21
 */
#include "RegionParticleIteratorTest.h"

#include "IteratorTestHelper.h"
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayUtils.h"
#include "testingHelpers/EmptyFunctor.h"

using namespace autopas;

template <typename AutoPasT>
auto RegionParticleIteratorTest::defaultInit(AutoPasT &autoPas, autopas::ContainerOption &containerOption,
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

#ifdef AUTOPAS_CUDA
  autoPas.setVerletClusterSize(32);
#endif

  autoPas.init();

  auto haloBoxMin =
      autopas::utils::ArrayMath::subScalar(autoPas.getBoxMin(), autoPas.getVerletSkin() + autoPas.getCutoff());
  auto haloBoxMax =
      autopas::utils::ArrayMath::addScalar(autoPas.getBoxMax(), autoPas.getVerletSkin() + autoPas.getCutoff());

  return std::make_tuple(haloBoxMin, haloBoxMax);
}

TEST_P(RegionParticleIteratorTest, testRegionAroundCorner) {
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
      IteratorTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
  }

  std::vector<size_t> expectedIDs;
  switch (behavior) {
    case autopas::IteratorBehavior::ownedOnly: {
      expectedIDs = particleIDsInBoxOwned;
      break;
    }
    case autopas::IteratorBehavior::haloOnly: {
      expectedIDs = particleIDsInBoxHalo;
      break;
    }
    case autopas::IteratorBehavior::haloAndOwned: {
      expectedIDs = particleIDsInBoxOwned;
      expectedIDs.insert(expectedIDs.end(), particleIDsInBoxHalo.begin(), particleIDsInBoxHalo.end());
      break;
    }
    case autopas::IteratorBehavior::haloOwnedAndDummy: {
      GTEST_FAIL() << "IteratorBehavior::haloOwnedAndDummy should not be tested through this test"
                      " as container behavior with dummy particles is not uniform.";
      break;
    }
  }

  // sanity check: there should be particles in the expected region
  ASSERT_THAT(expectedIDs, ::testing::Not(::testing::IsEmpty()));

  // actual test
  IteratorTestHelper::provideRegionIterator(
      useConstIterator, autoPas, behavior, searchBoxMin, searchBoxMax,
      [&](const auto &autopas, auto &iter) { IteratorTestHelper::findParticles(autoPas, iter, expectedIDs); });
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
  containerOptions.erase(containerOptions.find(autopas::ContainerOption::verletClusterCells));
  return containerOptions;
#endif
}

static inline auto getIteratorBehaviorOptions() {
  auto allOptions = autopas::IteratorBehavior::getAllOptions();
  std::set<autopas::IteratorBehavior> retSet;
  // we ignore dummy particles in the general tests because they can behave differently depending on the container
  std::set<autopas::IteratorBehavior> ignoredOptions = {autopas::IteratorBehavior::haloOwnedAndDummy};
  std::set_difference(allOptions.begin(), allOptions.end(), ignoredOptions.begin(), ignoredOptions.end(),
                      std::inserter(retSet, retSet.begin()));
  return retSet;
}

INSTANTIATE_TEST_SUITE_P(Generated, RegionParticleIteratorTest,
                         Combine(ValuesIn(getTestableContainerOptions()), /*cell size factor*/ Values(0.5, 1., 1.5),
                                 /*use const*/ Values(true, false), /*prior force calc*/ Values(true, false),
                                 ValuesIn(getIteratorBehaviorOptions())),
                         RegionParticleIteratorTest::PrintToStringParamName());
