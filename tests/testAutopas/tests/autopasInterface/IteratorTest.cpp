/**
 * @file IteratorTest.cpp
 * @author seckler
 * @date 22.07.19
 */

#include "IteratorTest.h"

#include <autopasTools/generators/RandomGenerator.h>
#include <testingHelpers/TouchableParticle.h>

#include "testingHelpers/commonTypedefs.h"

constexpr double cutoff = 1.;
constexpr double skin = 0.2;
constexpr std::array<double, 3> boxMin{0., 0., 0.};
constexpr std::array<double, 3> haloBoxMin{0. - skin - cutoff, 0. - skin - cutoff, 0. - skin - cutoff};
constexpr std::array<double, 3> boxMax{10., 10., 10.};
constexpr std::array<double, 3> haloBoxMax{10. + skin + cutoff, 10. + skin + cutoff, 10. + skin + cutoff};

template <typename AutoPasT>
void defaultInit(AutoPasT &autoPas) {
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(cutoff);
  autoPas.setVerletSkin(skin);
  autoPas.setVerletRebuildFrequency(2);
  autoPas.setNumSamples(2);
  // init autopas
  autoPas.init();
}
template <typename AutoPasT>
void checkRegionIteratorForAllParticles(AutoPasT &autoPas, autopas::IteratorBehavior behavior) {
  for (auto iter1 = autoPas.begin(behavior); iter1.isValid(); ++iter1) {
    unsigned int count = 0;
    auto low = iter1->getR();
    auto up = autopas::utils::ArrayMath::addScalar(low, 1e-10);

    for (auto iter2 = autoPas.getRegionIterator(low, up, behavior); iter2.isValid(); ++iter2) {
      ++count;
      EXPECT_EQ(&(*iter1), &(*iter2));
    }
    EXPECT_EQ(count, 1u);
  }
}

/**
 * Tests the addition and iteration over particles.
 * @param containerOption
 */
void testAdditionAndIteration(testingTuple options) {
  // create AutoPas object
  autopas::AutoPas<Molecule, FMCell> autoPas;

  auto [containerOption, cellSizeOption] = options;

  autoPas.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption});
  autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>({cellSizeOption})));

  defaultInit(autoPas);

  constexpr size_t numParticles1dOwned = 3;
  constexpr size_t numParticles1dTotal = 10;
  auto getPossible1DPositions = [&](double min, double max) -> auto {
    return std::array<double, numParticles1dTotal>{min - cutoff - skin + 1e-10,
                                                   min - cutoff,
                                                   min - skin / 4,
                                                   min,
                                                   min + skin / 4,
                                                   max - skin / 4,
                                                   max,
                                                   max + skin / 4,
                                                   max + cutoff,
                                                   max + cutoff + skin - 1e-10};
    // ensure that all particles are at most skin away from halo!
  };

  size_t id = 0;
  for (auto x : getPossible1DPositions(boxMin[0], boxMax[0])) {
    for (auto y : getPossible1DPositions(boxMin[1], boxMax[1])) {
      for (auto z : getPossible1DPositions(boxMin[2], boxMax[2])) {
        std::array<double, 3> pos{x, y, z};
        Molecule p(pos, {0., 0., 0.}, id, 0);
        ++id;
        // add the two particles!
        if (autopas::utils::inBox(pos, boxMin, boxMax)) {
          autoPas.addParticle(p);
        } else {
          autoPas.addOrUpdateHaloParticle(p);
        }
      }
    }
  }
  {
    size_t count = 0;
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      ++count;
      EXPECT_TRUE(iter->isOwned());
    }

    EXPECT_EQ(count, numParticles1dOwned * numParticles1dOwned * numParticles1dOwned);
  }

  // check number of halo particles
  {
    size_t count = 0;
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
      ++count;
      EXPECT_FALSE(iter->isOwned());
    }

    EXPECT_EQ(count, numParticles1dTotal * numParticles1dTotal * numParticles1dTotal -
                         numParticles1dOwned * numParticles1dOwned * numParticles1dOwned);
  }

  // check number of particles
  {
    size_t count = 0;
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
      ++count;
    }
    EXPECT_EQ(count, numParticles1dTotal * numParticles1dTotal * numParticles1dTotal);
  }

  // check number of halo particles for region iterator
  {
    size_t count = 0;
    for (auto iter = autoPas.getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::haloOnly);
         iter.isValid(); ++iter) {
      ++count;
      EXPECT_FALSE(iter->isOwned());
    }

    EXPECT_EQ(count, numParticles1dTotal * numParticles1dTotal * numParticles1dTotal -
                         numParticles1dOwned * numParticles1dOwned * numParticles1dOwned);
  }

  // check number of particles for region iterator
  {
    size_t count = 0;
    for (auto iter = autoPas.getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::haloAndOwned);
         iter.isValid(); ++iter) {
      ++count;
    }
    EXPECT_EQ(count, numParticles1dTotal * numParticles1dTotal * numParticles1dTotal);
  }

  // check number of owned particles for region iterator
  {
    size_t count = 0;
    for (auto iter = autoPas.getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::ownedOnly);
         iter.isValid(); ++iter) {
      ++count;
    }

    EXPECT_EQ(count, numParticles1dOwned * numParticles1dOwned * numParticles1dOwned);
  }

  // check all particles are in region iterator of their position, ownedOnly
  checkRegionIteratorForAllParticles(autoPas, autopas::IteratorBehavior::ownedOnly);

  // check all particles are in region iterator of their position, haloOnly
  checkRegionIteratorForAllParticles(autoPas, autopas::IteratorBehavior::haloOnly);

  // check all particles are in region iterator of their position, haloAndOwned
  checkRegionIteratorForAllParticles(autoPas, autopas::IteratorBehavior::haloAndOwned);

  // now move particles by at most skin/2, as this is the maximal distance they are allowed to move
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    auto pos = iter->getR();
    for (auto d = 0; d < 3; ++d) {
      if (pos[d] < boxMin[d]) {
        pos[d] += skin / 2;
      } else if (pos[d] >= boxMax[d]) {
        pos[d] -= skin / 2;
      } else if (pos[d] < (boxMax[d] + boxMin[d]) / 2) {
        pos[d] -= skin / 2;
      } else {
        pos[d] += skin / 2;
      }
    }
    iter->setR(pos);
  }

  // check all particles are in region iterator of their position, ownedOnly
  checkRegionIteratorForAllParticles(autoPas, autopas::IteratorBehavior::ownedOnly);

  // check all particles are in region iterator of their position, haloOnly
  checkRegionIteratorForAllParticles(autoPas, autopas::IteratorBehavior::haloOnly);

  // check all particles are in region iterator of their position, haloAndOwned
  checkRegionIteratorForAllParticles(autoPas, autopas::IteratorBehavior::haloAndOwned);
}

/**
 * Tests the equivalence of the range-based for loop with the normal for loop using isValid
 * @param containerOption
 */
void testRangeBasedIterator(testingTuple options) {
  // create AutoPas object
  autopas::AutoPas<Molecule, FMCell> autoPas;

  auto [containerOption, cellSizeOption] = options;

  autoPas.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption});
  autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>({cellSizeOption})));

  defaultInit(autoPas);

  constexpr size_t numParticles1dTotal = 10;
  auto getPossible1DPositions = [&](double min, double max) -> auto {
    return std::array<double, numParticles1dTotal>{min - cutoff - skin + 1e-10,
                                                   min - cutoff,
                                                   min - skin / 4,
                                                   min,
                                                   min + skin / 4,
                                                   max - skin / 4,
                                                   max,
                                                   max + skin / 4,
                                                   max + cutoff,
                                                   max + cutoff + skin - 1e-10};
    // ensure that all particles are at most skin away from halo!
  };

  size_t id = 0;
  for (auto x : getPossible1DPositions(boxMin[0], boxMax[0])) {
    for (auto y : getPossible1DPositions(boxMin[1], boxMax[1])) {
      for (auto z : getPossible1DPositions(boxMin[2], boxMax[2])) {
        std::array<double, 3> pos{x, y, z};
        Molecule p(pos, {0., 0., 0.}, id);
        ++id;
        // add the two particles!
        if (autopas::utils::inBox(pos, boxMin, boxMax)) {
          autoPas.addParticle(p);
        } else {
          autoPas.addOrUpdateHaloParticle(p);
        }
      }
    }
  }

  for (Particle &particle : autoPas) {
    particle.setF({42., 42., 42.});
  }

  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    decltype(iter->getF()) comparison = {42., 42., 42};
    ASSERT_EQ(iter->getF(), comparison);
  }
}

TEST_P(IteratorTest, ParticleAdditionAndIteratorTestNormal) {
  auto options = GetParam();
  testAdditionAndIteration(options);
}

TEST_P(IteratorTest, RangeBasedIterator) {
  auto options = GetParam();
  testRangeBasedIterator(options);
}

/**
 * The main idea of this test is to compare the iterators using openmp with the iterators not using openmp.
 */
void IteratorTest::testOpenMPIterators(autopas::ContainerOption containerOption, double cellSizeFactor,
                                       autopas::IteratorBehavior behavior, bool testRegionIterators) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {8, 8, 8};

  std::array<double, 3> lowCorner = {0, 0, 0};
  std::array<double, 3> highCorner = {3, 3, 3};

  int clusterSize = 64;
  autopas::AutoPas<TouchableParticle, autopas::FullParticleCell<TouchableParticle>> apContainer;

  apContainer.setAllowedContainers({containerOption});
  apContainer.setCellSizeFactor(cellSizeFactor);

  apContainer.setBoxMin(min);
  apContainer.setBoxMax(max);
  apContainer.setCutoff(cutoff);
  apContainer.setVerletSkin(skin);
  apContainer.setVerletClusterSize(clusterSize);

  apContainer.init();

  autopasTools::generators::RandomGenerator::fillWithParticles(apContainer, TouchableParticle({0., 0., 0.}, 0),
                                                               apContainer.getBoxMin(), apContainer.getBoxMax(), 500);
  autopasTools::generators::RandomGenerator::fillWithHaloParticles(
      apContainer, TouchableParticle({0., 0., 0.}, 0), cutoff, 50,
      [](decltype(apContainer) &c, TouchableParticle p) { c.addOrUpdateHaloParticle(p); });

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  {
    // with OpenMP:
    auto begin = testRegionIterators ? apContainer.getRegionIterator(lowCorner, highCorner, behavior)
                                     : apContainer.begin(behavior);
    for (auto iter = begin; iter.isValid(); ++iter) {
      iter->touch();
      if (behavior == autopas::IteratorBehavior::ownedOnly) {
        EXPECT_TRUE(iter->isOwned());
      } else if (behavior == autopas::IteratorBehavior::haloOnly) {
        EXPECT_FALSE(iter->isOwned());
      }
    }
  }
  {
    // without OpenMP:
    auto begin = testRegionIterators ? apContainer.getRegionIterator(lowCorner, highCorner, behavior)
                                     : apContainer.begin(behavior);
    for (auto iter = begin; iter.isValid(); ++iter) {
      EXPECT_EQ(1, iter->getNumTouched());
      if (behavior == autopas::IteratorBehavior::ownedOnly) {
        EXPECT_TRUE(iter->isOwned());
      } else if (behavior == autopas::IteratorBehavior::haloOnly) {
        EXPECT_FALSE(iter->isOwned());
      }
    }
  }
}

/**
 * Compare the OpenMP iterator behavior for owned only.
 */
TEST_P(IteratorTest, testOpenMPIteratorsOwnedOnly) {
  auto [containerOption, cellSizeFactor] = GetParam();
  testOpenMPIterators(containerOption, cellSizeFactor, autopas::IteratorBehavior::ownedOnly, false);
}

/**
 * Compare the OpenMP iterator behavior for halo and owned particles.
 */
TEST_P(IteratorTest, testOpenMPIteratorsHaloAndOwned) {
  auto [containerOption, cellSizeFactor] = GetParam();
  testOpenMPIterators(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloAndOwned, false);
}

/**
 * Compare the OpenMP iterator behavior for halo only.
 */
TEST_P(IteratorTest, testOpenMPIteratorsHaloOnly) {
  auto [containerOption, cellSizeFactor] = GetParam();
  testOpenMPIterators(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloOnly, false);
}

/**
 * Compare the OpenMP RegionIterator behavior for owned only.
 */
TEST_P(IteratorTest, testOpenMPRegionIteratorsOwnedOnly) {
  auto [containerOption, cellSizeFactor] = GetParam();
  testOpenMPIterators(containerOption, cellSizeFactor, autopas::IteratorBehavior::ownedOnly, true);
}

/**
 * Compare the OpenMP RegionIterator behavior for halo and owned particles.
 */
TEST_P(IteratorTest, testOpenMPRegionIteratorsHaloAndOwned) {
  auto [containerOption, cellSizeFactor] = GetParam();
  testOpenMPIterators(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloAndOwned, true);
}

/**
 * Compare the OpenMP RegionIterator behavior for halo only.
 */
TEST_P(IteratorTest, testOpenMPRegionIteratorsHaloOnly) {
  auto [containerOption, cellSizeFactor] = GetParam();
  testOpenMPIterators(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloOnly, true);
}

using ::testing::Combine;
using ::testing::UnorderedElementsAreArray;
using ::testing::Values;
using ::testing::ValuesIn;

INSTANTIATE_TEST_SUITE_P(Generated, IteratorTest,
                         // proper indent
                         Combine(ValuesIn([]() -> std::set<autopas::ContainerOption> {
                                   auto allContainerOptions = autopas::ContainerOption::getAllOptions();
                                   /// @TODO no verletClusterLists yet, so we erase it for now.
                                   allContainerOptions.erase(
                                       allContainerOptions.find(autopas::ContainerOption::verletClusterLists));
                                   return allContainerOptions;
                                 }()),
                                 Values(0.5, 1., 1.5)),
                         IteratorTest::PrintToStringParamName());
