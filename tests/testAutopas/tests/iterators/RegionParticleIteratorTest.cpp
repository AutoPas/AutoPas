/**
 * @file RegionParticleIteratorTest.cpp
 * @author F. Gratl
 * @date 08.03.21
 */
#include "RegionParticleIteratorTest.h"

#include "IteratorTestHelper.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/utils/WrapOpenMP.h"
#include "testingHelpers/EmptyPairwiseFunctor.h"
#include "testingHelpers/NumThreadGuard.h"

extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::computeInteractions(EmptyPairwiseFunctor<Molecule> *);

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::Values;
using ::testing::ValuesIn;

// Helper struct to represent region iterator boxes.
struct Box {
  std::array<double, 3> min;
  std::array<double, 3> max;
};

static inline auto getTestableContainerOptions() { return autopas::ContainerOption::getAllOptions(); }

template <typename AutoPasT>
auto RegionParticleIteratorTestBase::defaultInit(AutoPasT &autoPas, const autopas::ContainerOption &containerOption,
                                                 double cellSizeFactor) {
  using namespace autopas::utils::ArrayMath::literals;

  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1);
  autoPas.setVerletSkin(0.2);
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

/**
 * 1. Create an AutoPas container with 1000 particles that are placed around its 8 corners.
 * 2. Create a region iterator well around the lower corner of the container
 * 3. Run the region iterator for its full range and track the IDs it encounters
 * 4. Compare the found IDs to the expectations from the initialization.
 */
TEST_P(RegionParticleIteratorTestOne, testRegionAroundCorner) {
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
    EmptyPairwiseFunctor<Molecule> eFunctor;
    autoPas.computeInteractions(&eFunctor);
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

INSTANTIATE_TEST_SUITE_P(Generated, RegionParticleIteratorTestOne,
                         Combine(ValuesIn(getTestableContainerOptions()), /*cell size factor*/ Values(0.5, 1., 1.5),
                                 /*use const*/ Values(true, false), /*prior force calc*/ Values(true, false),
                                 ValuesIn(autopas::IteratorBehavior::getMostOptions())),
                         RegionParticleIteratorTestOne::PrintToStringParamName());

/**
 * Tests that AutoPas rejects regions where regionMin > regionMax.
 */
TEST_F(RegionParticleIteratorTestOne, testInvalidBox) {
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
TEST_F(RegionParticleIteratorTestOne, testForceSequential) {
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
    AUTOPAS_OPENMP(parallel for)
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

/**
 * This test checks whether a particle that has moved beyond the boundaries of a data structure element but is still
 * stored in this data structure element is found.
 */
TEST_P(RegionParticleIteratorTestTwo, testParticleMisplacement) {
  using namespace autopas::utils::ArrayMath::literals;

  auto containerOption = GetParam();

  autopas::AutoPas<Molecule> autoPas{};
  // Get a 10x10x10 box
  defaultInit(autoPas, containerOption, 1.);
  const auto shortDistance = autoPas.getVerletSkin() * 0.1;

  // This function creates particle positions within the simulation box that are the offset-value distant from the
  // boundary and the positions mirrored at the boundary. In addition, a search box is calculated for each position.
  auto generateParticlePositions = [&autoPas, &shortDistance]() {
    using autopas::utils::Math::isNearRel;
    constexpr size_t numParticles1DTotal = 3;
    constexpr double margin = 1e-10;

    auto cutoff = autoPas.getCutoff();
    auto skin = autoPas.getVerletSkin();
    auto boxMin = autoPas.getBoxMin();
    auto boxMax = autoPas.getBoxMax();
    auto offset = shortDistance;

    // this creates the positions (min + offset) and (max - offset) and positions mirrored at min and max (along one
    // dimension)
    auto generateInteresting1DPositions = [&](double min, double max) -> auto{
      return std::array<std::tuple<double, double>, numParticles1DTotal>{
          {{min + offset, min - offset}, {(max - min) / 2.0, (max - min) / 2.0}, {max - offset, max + offset}}};
    };

    std::vector<std::array<double, 3>> positionsInsideBox;
    std::vector<std::array<double, 3>> positionsMirroredAtBoundary;
    std::vector<std::tuple<Box, Box>> searchRegions;

    for (auto [x, xMirrored] : generateInteresting1DPositions(boxMin[0], boxMax[0])) {
      for (auto [y, yMirrored] : generateInteresting1DPositions(boxMin[1], boxMax[1])) {
        for (auto [z, zMirrored] : generateInteresting1DPositions(boxMin[2], boxMax[2])) {
          if (not(isNearRel(x, (boxMax[0] - boxMin[0]) / 2.0) and isNearRel(y, (boxMax[1] - boxMin[1]) / 2.0) and
                  isNearRel(z, (boxMax[2] - boxMin[2]) / 2.0))) {
            positionsInsideBox.push_back({x, y, z});
            positionsMirroredAtBoundary.push_back({xMirrored, yMirrored, zMirrored});

            // this ensures that a value that is closer than offset to the floored integer, is floored to the next
            // smaller integer
            auto getFlooredValue = [offset](double value) {
              return ((value - floor(value) < (offset - margin)) and (value - floor(value) >= 0)) ? floor(value - 1)
                                                                                                  : floor(value);
            };

            // this ensures that a value that is closer than offset to the ceiled integer, is ceiled to the next greater
            // integer
            auto getCeiledValue = [offset](double value) {
              return ((ceil(value) - value < (offset - margin)) and (ceil(value) - value >= 0)) ? ceil(value + 1)
                                                                                                : ceil(value);
            };

            // store corresponding search boxes for particles and mirrored versions. Note: we add a small margin to the
            // boxes to ensure that only the datastructure element, that actually contains a particle is selected in
            // getParticleImpl
            const Box startRegion{
                {getFlooredValue(x) + margin, getFlooredValue(y) + margin, getFlooredValue(z) + margin},
                {getCeiledValue(x) - margin, getCeiledValue(y) - margin, getCeiledValue(z) - margin}};
            const Box endRegion{{getFlooredValue(xMirrored) + margin, getFlooredValue(yMirrored) + margin,
                                 getFlooredValue(zMirrored) + margin},
                                {getCeiledValue(xMirrored) - margin, getCeiledValue(yMirrored) - margin,
                                 getCeiledValue(zMirrored) - margin}};
            searchRegions.push_back({startRegion, endRegion});
          }
        }
      }
    }
    return std::make_tuple(positionsInsideBox, positionsMirroredAtBoundary, searchRegions);
  };

  const auto [positionsInsideBox, positionsMirroredAtBoundary, searchRegions] = generateParticlePositions();

  std::vector<size_t> IDs;
  for (size_t i{0}; i < positionsInsideBox.size(); ++i) {
    autoPas.addParticle({positionsInsideBox[i], {0., 0., 0.}, i});
    IDs.push_back(i);
  }

  // sanity check for search boxes
  for (const auto searchRegion : searchRegions) {
    const auto startRegion = std::get<0>(searchRegion);
    const auto endRegion = std::get<1>(searchRegion);
    ASSERT_FALSE(autopas::utils::boxesOverlap(startRegion.min, startRegion.max, endRegion.min, endRegion.max))
        << "Search boxes should not overlap. If they do the test is written wrong.";
  }

  // Helper function to check if a region iterator finds a given particle with exptectedID in the given region.
  auto testRegion = [&](const std::array<double, 3> &min, const std::array<double, 3> &max, size_t exptectedID,
                        const std::string &context) {
    size_t numParticlesFound = 0;
    for (auto iter = autoPas.getRegionIterator(min, max, autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      ++numParticlesFound;
      EXPECT_EQ(iter->getID(), exptectedID) << "There should only be one particle with ID 0.\n" << context;
    }
    EXPECT_EQ(numParticlesFound, 1) << "Exactly one particle was inserted in the domain\n" << context;
  };

  // Actual test section
  auto leavingParticles = autoPas.updateContainer();
  // Ensure that particle is found with the given search region before it is moved out of the datastructure
  for (size_t i{0}; i < searchRegions.size(); ++i) {
    const auto startRegion = std::get<0>(searchRegions[i]);
    testRegion(startRegion.min, startRegion.max, IDs[i], "Before particle is moved.");
  }

  // This increments the stepsSinceLastRebuild counter in the LogicHandler and container, which is needed in
  // getParticleImpl of the containers to calculate the correct boxMin and boxMax values
  // Note: It is important to perform this step before the particles are moved. Otherwise they would be sorted into the
  // correct data structure element and the purpose of this test would be lost
  EmptyPairwiseFunctor<Molecule> eFunctor;
  autoPas.computeInteractions(&eFunctor);
  leavingParticles = autoPas.updateContainer();

  // Move the particles outside the simulations box and thus outside the data structure element (e.g: cell) where it is
  // currently stored
  for (auto iter = autoPas.begin(); iter != autoPas.end(); ++iter) {
    iter->setR(positionsMirroredAtBoundary[iter->getID()]);
  }

  // Now we check again with the corresponding search regions if the particles are found.
  for (size_t i{0}; i < searchRegions.size(); ++i) {
    const auto endRegion = std::get<1>(searchRegions[i]);
    testRegion(endRegion.min, endRegion.max, IDs[i], "After particle was moved.");
  }
}

INSTANTIATE_TEST_SUITE_P(Generated, RegionParticleIteratorTestTwo, ValuesIn(getTestableContainerOptions()),
                         RegionParticleIteratorTestTwo::PrintToStringParamName());