/**
 * @file ParticleIteratorTest.cpp
 * @author F. Gratl
 * @date 08.03.21
 */

#include "ParticleIteratorTest.h"

#include "IteratorTestHelper.h"
#include "autopas/iterators/ParticleIterator.h"
#include "autopas/options/IteratorBehavior.h"
#include "testingHelpers/commonTypedefs.h"

using namespace autopas;
using namespace autopas::internal;

void ParticleIteratorTest::testAllParticlesFoundPattern(const std::vector<size_t> &cellsWithParticles,
                                                        size_t numThreads, size_t numAdditionalParticleVectors) {
  constexpr size_t particlesPerCell = 4;
  auto cells = IteratorTestHelper::generateCellsWithPattern(10, cellsWithParticles, particlesPerCell);
  std::vector<std::vector<Molecule>> additionalParticles(numAdditionalParticleVectors);
  auto particlesTotal = cellsWithParticles.size() * particlesPerCell;
  // fill every additional particle buffer with particles
  for (auto &additionalParticleVector : additionalParticles) {
    for (size_t i = 0; i < particlesPerCell; ++i) {
      additionalParticleVector.push_back(Molecule({1., 2., 3.}, {0, 0, 0}, particlesTotal++, 0));
    }
  }

  std::vector<size_t> foundParticles;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel num_threads(numThreads) reduction(vecMerge : foundParticles)
#endif
  {
    for (auto iter = ParticleIterator<Molecule, FMCell, true>(
             &cells, 0, nullptr, IteratorBehavior::ownedOrHalo,
             numAdditionalParticleVectors > 0 ? &additionalParticles : nullptr);
         iter.isValid(); ++iter) {
      auto particleID = iter->getID();
      foundParticles.push_back(particleID);
    }
  }

  // expect all indices from 1..N to be found
  std::vector<size_t> expectedIndices(particlesTotal);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 0);
  ASSERT_THAT(foundParticles, ::testing::UnorderedElementsAreArray(expectedIndices));
}

TEST_P(ParticleIteratorTest, testFullIterator_EFEFFEEFEF) {
  auto [numThreads, numAdditionalParticleVectors] = GetParam();
  // Empty Full Empty Full Full Empty Empty Full Empty Full
  testAllParticlesFoundPattern({1ul, 3ul, 4ul, 7ul, 9ul}, numThreads, numAdditionalParticleVectors);
}

TEST_P(ParticleIteratorTest, testFullIterator_FEFEEFFEFE) {
  auto [numThreads, numAdditionalParticleVectors] = GetParam();
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  testAllParticlesFoundPattern({0ul, 2ul, 5ul, 6ul, 8ul}, numThreads, numAdditionalParticleVectors);
}

TEST_F(ParticleIteratorTest, testFullIterator_deletion) {
  // Full Empty Full Empty Empty Full Full Empty Full Empty
  const std::vector<size_t> cellsToFill = {0ul, 2ul, 5ul, 6ul, 8ul};
  const size_t numParticlesToAddPerCell = 4;
  auto data = IteratorTestHelper::generateCellsWithPattern(10, cellsToFill, numParticlesToAddPerCell);

  int numFoundParticles = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : numFoundParticles)
#endif
  {
    ParticleIterator<Molecule, FMCell, true> iter(&data);
    for (; iter.isValid(); ++iter, ++numFoundParticles) {
      autopas::internal::deleteParticle(iter);
    }
  }
  ASSERT_EQ(numFoundParticles, cellsToFill.size() * numParticlesToAddPerCell);

  numFoundParticles = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : numFoundParticles)
#endif
  {
    ParticleIterator<Molecule, FMCell, true> iter(&data);
    for (; iter.isValid(); ++iter) {
      ++numFoundParticles;
    }
  }
  ASSERT_EQ(numFoundParticles, 0);
}

#ifdef AUTOPAS_OPENMP
const std::vector<size_t> threadNumsToTest{1, 2, 4};
#else
// no need to test more thread counts when openmp is not enabled
const std::vector<size_t> threadNumsToTest{1};
#endif

INSTANTIATE_TEST_SUITE_P(Generated, ParticleIteratorTest,
                         ::testing::Combine(::testing::ValuesIn(threadNumsToTest), ::testing::Values(0, 1, 2, 4)),
                         ParticleIteratorTest::PrintToStringParamName());

template <class Iter>
auto ParticleIteratorTest::iteratorsBehaveEqually(Iter &iter1, Iter &iter2) {
  for (; iter1.isValid() and iter2.isValid(); ++iter1, ++iter2) {
    EXPECT_EQ(*iter1, *iter2);
  }
}

TEST_F(ParticleIteratorTest, testCopyConstructor) {
  constexpr size_t particlesPerCell = 4;
  auto cells = IteratorTestHelper::generateCellsWithPattern(10, {1ul, 3ul, 4ul, 7ul, 9ul}, particlesPerCell);

  std::vector<std::vector<Molecule>> additionalVectors(3);
  additionalVectors[1].emplace_back(Molecule({0., 0., 0.}, {0., 0., 0.}, 1337));

  constexpr bool modifyable = true;
  autopas::internal::ParticleIterator<Molecule, FMCell, modifyable> iter(
      &cells, 0, nullptr, IteratorBehavior::ownedOrHalo, &additionalVectors);

  auto iterCopy{iter};
  iteratorsBehaveEqually(iter, iterCopy);
}

TEST_F(ParticleIteratorTest, testCopyAssignment) {
  constexpr size_t particlesPerCell = 4;
  auto cells = IteratorTestHelper::generateCellsWithPattern(10, {1ul, 3ul, 4ul, 7ul, 9ul}, particlesPerCell);

  std::vector<std::vector<Molecule>> additionalVectors(3);
  additionalVectors[1].emplace_back(Molecule({0., 0., 0.}, {0., 0., 0.}, 1337));

  constexpr bool modifyable = true;
  autopas::internal::ParticleIterator<Molecule, FMCell, modifyable> iter(
      &cells, 0, nullptr, IteratorBehavior::ownedOrHalo, &additionalVectors);

  auto iterCopy = iter;
  iteratorsBehaveEqually(iter, iterCopy);
}

/**
 * Generates an iterator in a parallel region but iterates with only one and expects to find everything.
 * @note This behavior is needed by VerletClusterLists::updateHaloParticle().
 */
TEST_F(ParticleIteratorTest, testForceSequential) {
  constexpr size_t particlesPerCell = 1;
  auto cells = IteratorTestHelper::generateCellsWithPattern(4, {0ul, 1ul, 2ul, 3ul}, particlesPerCell);
  auto particlesTotal = cells.size() * particlesPerCell;
  std::vector<size_t> expectedIndices(particlesTotal);
  std::iota(expectedIndices.begin(), expectedIndices.end(), 0);

  constexpr size_t numAdditionalVectors = 3;
  std::vector<std::vector<Molecule>> additionalVectors(numAdditionalVectors);

  size_t particleId = expectedIndices.size() + 100;
  for (auto &vector : additionalVectors) {
    vector.emplace_back(Molecule({0., 0., 0.}, {0., 0., 0.}, particleId));
    expectedIndices.push_back(particleId);
    ++particleId;
  }

#pragma omp parallel
  {
    std::vector<size_t> foundParticles;
    constexpr bool modifyable = true;
    autopas::internal::ParticleIterator<Molecule, FMCell, modifyable> iter(
        &cells, 0, nullptr, IteratorBehavior::ownedOrHalo | IteratorBehavior::forceSequential, &additionalVectors);
    for (; iter.isValid(); ++iter) {
      foundParticles.push_back(iter->getID());
    }
    EXPECT_THAT(foundParticles, ::testing::UnorderedElementsAreArray(expectedIndices));
  }
}

/**
 * Make sure the iterator does not crash and immediately is marked invalid when offset behind the last cell.
 */
TEST_F(ParticleIteratorTest, testNothing) {
  std::vector<FMCell> cells{5};
  ParticleIterator<Molecule, FMCell, true> iter{&cells, cells.size()};
  EXPECT_FALSE(iter.isValid());
}

/**
 * Same as testNothing but also with additional particle vectors.
 */
TEST_F(ParticleIteratorTest, testNothingWithAdditionalVectors) {
  std::vector<FMCell> cells{5};
  std::vector<std::vector<Molecule>> additionalPartilcleVectors(3);
  ParticleIterator<Molecule, FMCell, true> iter{&cells, cells.size(), nullptr, IteratorBehavior::ownedOrHaloOrDummy,
                                                &additionalPartilcleVectors};
  EXPECT_FALSE(iter.isValid());
}