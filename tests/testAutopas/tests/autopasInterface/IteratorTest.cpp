/**
 * @file IteratorTest.cpp
 * @author seckler
 * @date 22.07.19
 */

#include "IteratorTest.h"

#include "autopasTools/generators/RandomGenerator.h"
#include "testingHelpers/EmptyFunctor.h"
#include "testingHelpers/TouchableParticle.h"
#include "testingHelpers/commonTypedefs.h"

constexpr double cutoff = 1.;
constexpr double skin = 0.2;
constexpr std::array<double, 3> boxMin{0., 0., 0.};
constexpr std::array<double, 3> haloBoxMin{0. - skin - cutoff, 0. - skin - cutoff, 0. - skin - cutoff};
constexpr std::array<double, 3> boxMax{10., 10., 10.};
constexpr std::array<double, 3> haloBoxMax{10. + skin + cutoff, 10. + skin + cutoff, 10. + skin + cutoff};

using ::testing::_;

template <typename AutoPasT>
void defaultInit(AutoPasT &autoPas) {
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(cutoff);
  autoPas.setVerletSkin(skin);
  autoPas.setVerletRebuildFrequency(2);
  autoPas.setNumSamples(2);

#ifdef AUTOPAS_CUDA
  autoPas.setVerletClusterSize(32);
#endif
  // init autopas
  autoPas.init();
}

/**
 * Iterate over all particles, generate a region iterator for each that spans a tiny space around them and check if this
 * region iterator finds exactly this particle.
 * @tparam AutoPasT
 * @param autoPas
 * @param behavior
 */
template <typename AutoPasT>
void checkRegionIteratorForAllParticles(AutoPasT &autoPas, autopas::IteratorBehavior behavior) {
  for (auto iter1 = autoPas.begin(behavior); iter1.isValid(); ++iter1) {
    unsigned int count = 0;
    auto low = iter1->getR();
    auto up = autopas::utils::ArrayMath::addScalar(low, 1e-10);

    for (auto iter2 = autoPas.getRegionIterator(low, up, behavior); iter2.isValid(); ++iter2) {
      ++count;
      EXPECT_EQ(&(*iter1), &(*iter2)) << "Wrong particle found";
    }
    EXPECT_EQ(count, 1u) << "The following particle was not found exactly once:\n" << iter1->toString();
  }
}

/**
 * Tests the addition and iteration over particles.
 * @param containerOption
 */
template <bool testConstIterators>
void testAdditionAndIteration(autopas::ContainerOption containerOption, double cellSizeOption, bool priorForceCalc) {
  // create AutoPas object
  autopas::AutoPas<Molecule, FMCell> autoPas;

  // Reference to the AutoPas object to be able to check const iterators.
  std::conditional_t<testConstIterators, const autopas::AutoPas<Molecule, FMCell> &,
                     autopas::AutoPas<Molecule, FMCell> &>
      autoPasRef = autoPas;

  autoPas.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption});
  autoPas.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(containerOption));
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
  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule, FMCell> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
  }
  {
    size_t count = 0;
    for (auto iter = autoPasRef.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
      ++count;
      EXPECT_TRUE(iter->isOwned()) << "Iterator: Found a particle that is not owned while looking only for owned ones!";
    }

    EXPECT_EQ(count, numParticles1dOwned * numParticles1dOwned * numParticles1dOwned)
        << "Iterator: Found incorrect number of owned particles!";
  }

  // check number of halo particles
  {
    size_t count = 0;
    for (auto iter = autoPasRef.begin(autopas::IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
      ++count;
      EXPECT_FALSE(iter->isOwned())
          << "Iterator: Found a particle that is owned while looking only for not owned ones!";
    }

    EXPECT_EQ(count, numParticles1dTotal * numParticles1dTotal * numParticles1dTotal -
                         numParticles1dOwned * numParticles1dOwned * numParticles1dOwned)
        << "Iterator: Found incorrect number of halo only particles!";
  }

  // check number of particles
  {
    size_t count = 0;
    for (auto iter = autoPasRef.begin(autopas::IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
      ++count;
    }
    EXPECT_EQ(count, numParticles1dTotal * numParticles1dTotal * numParticles1dTotal)
        << "Iterator: Found incorrect number of halo + owned particles!";
  }

  // check number of halo particles for region iterator
  {
    size_t count = 0;
    for (auto iter = autoPasRef.getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::haloOnly);
         iter.isValid(); ++iter) {
      ++count;
      EXPECT_FALSE(iter->isOwned())
          << "RegionIterator: Found a particle that is owned while looking only for not owned ones!";
    }

    EXPECT_EQ(count, numParticles1dTotal * numParticles1dTotal * numParticles1dTotal -
                         numParticles1dOwned * numParticles1dOwned * numParticles1dOwned)
        << "RegionIterator: Found incorrect number of halo only particles!";
  }

  // check number of particles for region iterator
  {
    size_t count = 0;
    for (auto iter = autoPasRef.getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::haloAndOwned);
         iter.isValid(); ++iter) {
      ++count;
    }
    EXPECT_EQ(count, numParticles1dTotal * numParticles1dTotal * numParticles1dTotal)
        << "RegionIterator: Found incorrect number of halo + owned particles!";
  }

  // check number of owned particles for region iterator
  {
    size_t count = 0;
    for (auto iter = autoPasRef.getRegionIterator(haloBoxMin, haloBoxMax, autopas::IteratorBehavior::ownedOnly);
         iter.isValid(); ++iter) {
      ++count;
      EXPECT_TRUE(iter->isOwned())
          << "RegionIterator: Found a particle that is not owned while looking only for owned ones!";
    }

    EXPECT_EQ(count, numParticles1dOwned * numParticles1dOwned * numParticles1dOwned)
        << "RegionIterator: Found incorrect number of owned particles!";
  }

  // check all particles are in region iterator of their position, ownedOnly
  checkRegionIteratorForAllParticles(autoPasRef, autopas::IteratorBehavior::ownedOnly);

  // check all particles are in region iterator of their position, haloOnly
  checkRegionIteratorForAllParticles(autoPasRef, autopas::IteratorBehavior::haloOnly);

  // check all particles are in region iterator of their position, haloAndOwned
  checkRegionIteratorForAllParticles(autoPasRef, autopas::IteratorBehavior::haloAndOwned);

  // Now move particles by at most skin/2, as this is the maximal distance they are allowed to move.
  // Don't use autoPasRef here, as we are modifying things!
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
  checkRegionIteratorForAllParticles(autoPasRef, autopas::IteratorBehavior::ownedOnly);

  // check all particles are in region iterator of their position, haloOnly
  checkRegionIteratorForAllParticles(autoPasRef, autopas::IteratorBehavior::haloOnly);

  // check all particles are in region iterator of their position, haloAndOwned
  checkRegionIteratorForAllParticles(autoPasRef, autopas::IteratorBehavior::haloAndOwned);
}

/**
 * Tests the equivalence of the range-based for loop with the normal for loop using isValid
 * @param containerOption
 */
template <bool testConstIterators>
void testRangeBasedIterator(autopas::ContainerOption containerOption, double cellSizeOption, bool priorForceCalc) {
  // create AutoPas object
  autopas::AutoPas<Molecule, FMCell> autoPas;
  // Reference to the AutoPas object to be able to check const iterators.
  std::conditional_t<testConstIterators, const autopas::AutoPas<Molecule, FMCell> &,
                     autopas::AutoPas<Molecule, FMCell> &>
      autoPasRef = autoPas;

  autoPas.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption});
  autoPas.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(containerOption));
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

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule, FMCell> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
  }

  for (Particle &particle : autoPas) {
    particle.setF({42., 42., 42.});
  }

  for (auto &particle : autoPasRef) {
    decltype(particle.getF()) comparison = {42., 42., 42};
    ASSERT_EQ(particle.getF(), comparison);
  }

  for (auto iter = autoPasRef.begin(); iter.isValid(); ++iter) {
    decltype(iter->getF()) comparison = {42., 42., 42};
    ASSERT_EQ(iter->getF(), comparison);
  }
}

TEST_P(IteratorTest, ParticleAdditionAndIteratorTestNormal) {
  auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc] = GetParam();
  if (testConstIterators) {
    testAdditionAndIteration<true>(containerOption, cellSizeFactor, priorForceCalc);
  } else {
    testAdditionAndIteration<false>(containerOption, cellSizeFactor, priorForceCalc);
  }
}

TEST_P(IteratorTest, RangeBasedIterator) {
  auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc] = GetParam();
  if (testConstIterators) {
    testRangeBasedIterator<true>(containerOption, cellSizeFactor, priorForceCalc);
  } else {
    testRangeBasedIterator<false>(containerOption, cellSizeFactor, priorForceCalc);
  }
}

/**
 * The main idea of this test is to compare the iterators using openmp with the iterators not using openmp.
 * If OPENMP is disabled, this tests mainly that no particle is traversed twice.
 */
template <bool testConstIterators>
void IteratorTest::testOpenMPIterators(autopas::ContainerOption containerOption, double cellSizeFactor,
                                       autopas::IteratorBehavior behavior, bool testRegionIterators,
                                       bool priorForceCalc) {
  std::array<double, 3> min = {1, 1, 1};
  std::array<double, 3> max = {8, 8, 8};

  std::array<double, 3> lowCorner = {0, 0, 0};
  std::array<double, 3> highCorner = {3, 3, 3};

  int clusterSize = 64;
  autopas::AutoPas<TouchableParticle, autopas::FullParticleCell<TouchableParticle>> apContainer;
  // Reference to the AutoPas object to be able to check const iterators.
  std::conditional_t<testConstIterators,
                     const autopas::AutoPas<TouchableParticle, autopas::FullParticleCell<TouchableParticle>> &,
                     autopas::AutoPas<TouchableParticle, autopas::FullParticleCell<TouchableParticle>> &>
      autoPasRef = apContainer;

  apContainer.setAllowedContainers({containerOption});
  apContainer.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(containerOption));
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

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<TouchableParticle, autopas::FullParticleCell<TouchableParticle>> eFunctor;
    apContainer.iteratePairwise(&eFunctor);
  }

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
    auto begin = testRegionIterators ? autoPasRef.getRegionIterator(lowCorner, highCorner, behavior)
                                     : autoPasRef.begin(behavior);
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
  auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc] = GetParam();
  if (testConstIterators) {
    testOpenMPIterators<true>(containerOption, cellSizeFactor, autopas::IteratorBehavior::ownedOnly, false,
                              priorForceCalc);
  } else {
    testOpenMPIterators<false>(containerOption, cellSizeFactor, autopas::IteratorBehavior::ownedOnly, false,
                               priorForceCalc);
  }
}

/**
 * Compare the OpenMP iterator behavior for halo and owned particles.
 */
TEST_P(IteratorTest, testOpenMPIteratorsHaloAndOwned) {
  auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc] = GetParam();
  if (testConstIterators) {
    testOpenMPIterators<true>(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloAndOwned, false,
                              priorForceCalc);
  } else {
    testOpenMPIterators<false>(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloAndOwned, false,
                               priorForceCalc);
  }
}

/**
 * Compare the OpenMP iterator behavior for halo only.
 */
TEST_P(IteratorTest, testOpenMPIteratorsHaloOnly) {
  auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc] = GetParam();
  if (testConstIterators) {
    testOpenMPIterators<true>(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloOnly, false,
                              priorForceCalc);
  } else {
    testOpenMPIterators<false>(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloOnly, false,
                               priorForceCalc);
  }
}

/**
 * Compare the OpenMP RegionIterator behavior for owned only.
 */
TEST_P(IteratorTest, testOpenMPRegionIteratorsOwnedOnly) {
  auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc] = GetParam();
  if (testConstIterators) {
    testOpenMPIterators<true>(containerOption, cellSizeFactor, autopas::IteratorBehavior::ownedOnly, true,
                              priorForceCalc);
  } else {
    testOpenMPIterators<false>(containerOption, cellSizeFactor, autopas::IteratorBehavior::ownedOnly, true,
                               priorForceCalc);
  }
}

/**
 * Compare the OpenMP RegionIterator behavior for halo and owned particles.
 */
TEST_P(IteratorTest, testOpenMPRegionIteratorsHaloAndOwned) {
  auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc] = GetParam();
  if (testConstIterators) {
    testOpenMPIterators<true>(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloAndOwned, true,
                              priorForceCalc);
  } else {
    testOpenMPIterators<false>(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloAndOwned, true,
                               priorForceCalc);
  }
}

/**
 * Compare the OpenMP RegionIterator behavior for halo only.
 */
TEST_P(IteratorTest, testOpenMPRegionIteratorsHaloOnly) {
  auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc] = GetParam();
  if (testConstIterators) {
    testOpenMPIterators<true>(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloOnly, true,
                              priorForceCalc);
  } else {
    testOpenMPIterators<false>(containerOption, cellSizeFactor, autopas::IteratorBehavior::haloOnly, true,
                               priorForceCalc);
  }
}

/**
 * The main idea of this test is to check whether deletion through a RegionIterator doesn't do stupid things.
 */
template <bool testConstIterators>
void testRegionIteratorDeletion(autopas::ContainerOption containerOption, double cellSizeFactor, bool priorForceCalc) {
  // create AutoPas object
  autopas::AutoPas<Molecule, FMCell> autoPas;
  // Reference to the AutoPas object to be able to check const iterators.
  std::conditional_t<testConstIterators, const autopas::AutoPas<Molecule, FMCell> &,
                     autopas::AutoPas<Molecule, FMCell> &>
      autoPasRef = autoPas;

  autoPas.setAllowedContainers(std::set<autopas::ContainerOption>{containerOption});
  autoPas.setAllowedTraversals(autopas::compatibleTraversals::allCompatibleTraversals(containerOption));
  autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>({cellSizeFactor})));

  defaultInit(autoPas);

  size_t numParticles = 200;

  srand(42);
  for (size_t id = 0; id < numParticles; ++id) {
    auto pos = autopasTools::generators::RandomGenerator::randomPosition(haloBoxMin, haloBoxMax);
    Molecule p(pos, {0., 0., 0.}, id);
    // add the particle!
    if (autopas::utils::inBox(pos, boxMin, boxMax)) {
      autoPas.addParticle(p);
    } else {
      autoPas.addOrUpdateHaloParticle(p);
    }
  }

  if (priorForceCalc) {
    // the prior force calculation is partially wanted as this sometimes changes the state of the internal containers.
    EmptyFunctor<Molecule, FMCell> eFunctor;
    autoPas.iteratePairwise(&eFunctor);
  }

  std::set<unsigned long> beforeParticles;
  for (auto &particle : autoPasRef) {
    beforeParticles.insert(particle.getID());
  }
  ASSERT_EQ(beforeParticles.size(), numParticles);

  // delete some particles
  std::set<unsigned long> deletedParticles;
  for (auto it = autoPas.getRegionIterator({0., 0., 3.}, {10., 10., 4.5}); it != autoPas.end(); ++it) {
    deletedParticles.insert(it->getID());
    autoPas.deleteParticle(it);
  }

  // calculate the particles that should still be there.
  std::set<unsigned long> shouldBeParticles;
  for (auto id : beforeParticles) {
    if (deletedParticles.find(id) == deletedParticles.end()) {
      // Has not been deleted!
      shouldBeParticles.insert(id);
    }
  }

  // check whether they are still there.
  std::vector<unsigned long> afterParticles;
  for (auto &particle : autoPasRef) {
    EXPECT_NE(shouldBeParticles.find(particle.getID()), shouldBeParticles.end()) << "id:" << particle.getID();
    afterParticles.push_back(particle.getID());
  }

  EXPECT_EQ(afterParticles.size(), shouldBeParticles.size());
}

TEST_P(IteratorTest, RegionIteratorDeletion) {
  auto [containerOption, cellSizeFactor, testConstIterators, priorForceCalc] = GetParam();
  if (testConstIterators) {
    testRegionIteratorDeletion<true>(containerOption, cellSizeFactor, priorForceCalc);
  } else {
    testRegionIteratorDeletion<false>(containerOption, cellSizeFactor, priorForceCalc);
  }
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

INSTANTIATE_TEST_SUITE_P(Generated, IteratorTest,
                         Combine(ValuesIn(getTestableContainerOptions()), Values(0.5, 1., 1.5), Values(true, false),
                                 Values(true, false)),
                         IteratorTest::PrintToStringParamName());
