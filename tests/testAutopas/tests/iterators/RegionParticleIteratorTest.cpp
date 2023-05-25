/**
 * @file RegionParticleIteratorTest.cpp
 * @author F. Gratl
 * @date 08.03.21
 */
#include "RegionParticleIteratorTest.h"

#include "IteratorTestHelper.h"
#include "autopas/AutoPasDecl.h"
#include "testingHelpers/EmptyFunctor.h"
#include "testingHelpers/NumThreadGuard.h"

extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::iteratePairwise(EmptyFunctor<Molecule> *);

template <typename AutoPasT>
auto RegionParticleIteratorTest::defaultInit(AutoPasT &autoPas, const autopas::ContainerOption &containerOption,
                                             double cellSizeFactor) {
  using namespace autopas::utils::ArrayMath::literals;

  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1);
  autoPas.setVerletSkinPerTimestep(0.1);
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

/**
 * 1. Create an AutoPas container with 1000 particles that are placed around its 8 corners.
 * 2. Create a region iterator well around the lower corner of the container
 * 3. Run the region iterator for its full range and track the IDs it encounters
 * 4. Compare the found IDs to the expectations from the initialization.
 */
TEST_P(RegionParticleIteratorTest, testRegionAroundCorner) {
  using namespace autopas::utils::ArrayMath::literals;

  auto [containerOption, cellSizeFactor, useConstIterator, priorForceCalc, behavior] = GetParam();

  // init autopas and fill it with some particles
  autopas::AutoPas<Molecule> autoPas;
  defaultInit(autoPas, containerOption, cellSizeFactor);

  const auto domainLength = autoPas.getBoxMax() - autoPas.getBoxMin();
  // draw a box around the lower corner of the domain
  const auto searchBoxLengthHalf = domainLength * 0.3;
  const auto searchBoxMin = autoPas.getBoxMin() - searchBoxLengthHalf;
  const auto searchBoxMax = autoPas.getBoxMin() + searchBoxLengthHalf;

  // initialize particles and remember which IDs are in the search box
  const auto [particleIDsOwned, particleIDsHalo, particleIDsInBoxOwned, particleIDsInBoxHalo] =
      IteratorTestHelper::fillContainerAroundBoundary(autoPas, searchBoxMin, searchBoxMax);

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
  }

  // set up expectations
  // can't trivially convert this to const + lambda initialization bc behavior is a structured binding
  std::vector<size_t> expectedIDs;
  switch (behavior) {
    case autopas::IteratorBehavior::owned: {
      expectedIDs = particleIDsInBoxOwned;
      break;
    }
    case autopas::IteratorBehavior::halo: {
      expectedIDs = particleIDsInBoxHalo;
      break;
    }
    case autopas::IteratorBehavior::ownedOrHalo: {
      expectedIDs = particleIDsInBoxOwned;
      expectedIDs.insert(expectedIDs.end(), particleIDsInBoxHalo.begin(), particleIDsInBoxHalo.end());
      break;
    }
    default: {
      GTEST_FAIL() << "IteratorBehavior::" << behavior
                   << "  should not be tested through this test!\n"
                      "Container behavior with dummy particles is not uniform.\n"
                      "Using forceSequential is not supported.";
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

static inline auto getTestableContainerOptions() { return autopas::ContainerOption::getAllOptions(); }

INSTANTIATE_TEST_SUITE_P(Generated, RegionParticleIteratorTest,
                         Combine(ValuesIn(getTestableContainerOptions()), /*cell size factor*/ Values(0.5, 1., 1.5),
                                 /*use const*/ Values(true, false), /*prior force calc*/ Values(true, false),
                                 ValuesIn(autopas::IteratorBehavior::getMostOptions())),
                         RegionParticleIteratorTest::PrintToStringParamName());

/**
 * Tests that AutoPas rejects regions where regionMin > regionMax.
 */
TEST_F(RegionParticleIteratorTest, testInvalidBox) {
  // setup
  autopas::AutoPas<Molecule> autoPas{};
  const auto [haloBoxMin, haloBoxMax] = defaultInit(autoPas, autopas::ContainerOption::directSum, 1.);

  // helpers
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::sub;

  // calculate box size
  const std::array<double, 3> haloBoxLength = sub(haloBoxMax, haloBoxMin);
  const std::array<double, 3> haloBoxLength3rd = mulScalar(haloBoxMax, 1. / 3.);

  // calculate points within the domain
  const std::array<double, 3> regionUpperLimit = mulScalar(haloBoxLength3rd, 2.);
  const std::array<double, 3> regionLowerLimit = mulScalar(haloBoxLength3rd, 1.);

  // actual test
  EXPECT_NO_THROW(autoPas.getRegionIterator(regionLowerLimit, regionUpperLimit));
  EXPECT_THROW(autoPas.getRegionIterator(regionUpperLimit, regionLowerLimit),
               autopas::utils::ExceptionHandler::AutoPasException);
}

/**
 * Fills a container with (halo) particles, invokes region iterators with force sequential, and expects all of them to
 * find everything in the search region.
 */
TEST_F(RegionParticleIteratorTest, testForceSequential) {
  // helpers
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::sub;

  /// setup
  autopas::AutoPas<Molecule> autoPas{};
  const auto [haloBoxMin, haloBoxMax] = defaultInit(autoPas, autopas::ContainerOption::linkedCells, 1.);
  const auto haloBoxSize = sub(haloBoxMax, haloBoxMin);

  // define a search box covering the first octant of the domain
  const auto searchBoxMin = haloBoxMin;
  const auto searchBoxMax = mulScalar(haloBoxMax, 0.5);

  const auto particlesPerDimension = 9;
  const auto particleSpacing = div(haloBoxSize, {particlesPerDimension, particlesPerDimension, particlesPerDimension});

  // fill a container with a grid and track what is in the search region
  std::vector<size_t> idsInSearchRegion{};
  std::vector<size_t> idsNotInSearchRegion{};
  size_t id = 0;
  for (int z = 0; z < particlesPerDimension; ++z) {
    for (int y = 0; y < particlesPerDimension; ++y) {
      for (int x = 0; x < particlesPerDimension; ++x, ++id) {
        const auto pos = std::array<double, 3>{x * particleSpacing[0], y * particleSpacing[1], z * particleSpacing[2]};
        const Molecule m{pos, {}, id};
        // depending on the position add the particle as halo or owned
        if (autopas::utils::inBox(pos, autoPas.getBoxMin(), autoPas.getBoxMax())) {
          autoPas.addParticle(m);
        } else {
          autoPas.addHaloParticle(m);
        }
        // depending on the position track the particle id as in or out of the search box
        if (autopas::utils::inBox(pos, searchBoxMin, searchBoxMax)) {
          idsInSearchRegion.push_back(id);
        } else {
          idsNotInSearchRegion.push_back(id);
        }
      }
    }
  }

  /// Actual test: Have several threads iterate the region with force sequential.
  /// All should find everything in idsInSearchRegion.
  const size_t numThreads = 3;
  std::vector<std::vector<size_t>> encounteredIds(numThreads);
  {
    const NumThreadGuard numThreadGuard(numThreads);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
    for (int t = 0; t < numThreads; ++t) {
      for (auto iter = autoPas.getRegionIterator(
               searchBoxMin, searchBoxMax,
               autopas::IteratorBehavior::ownedOrHalo | autopas::IteratorBehavior::forceSequential);
           iter.isValid(); ++iter) {
        encounteredIds[t].push_back(iter->getID());
      }
    }
  }

  /// checks
  for (int t = 0; t < numThreads; ++t) {
    EXPECT_THAT(encounteredIds[t], ::testing::UnorderedElementsAreArray(idsInSearchRegion))
        << "Thread " << t << " did not find the correct IDs.";
  }
}
