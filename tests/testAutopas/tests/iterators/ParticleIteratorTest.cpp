/**
 * @file ParticleIteratorTest.cpp
 * @author tchipev
 * @date 19.01.18
 */

#include "ParticleIteratorTest.h"

#include "autopas/iterators/ParticleIterator.h"
#include "autopas/options/IteratorBehavior.h"
#include "testingHelpers/commonTypedefs.h"

using namespace autopas;
using namespace autopas::internal;

std::vector<FMCell> ParticleIteratorTest::generateCellsWithPattern(const size_t numCells,
                                                                   const std::vector<size_t> &cellsToFill,
                                                                   const size_t particlesPerCell) {
  std::vector<FMCell> cells(numCells);
  size_t numParticlesAdded = 0;
  for (auto cellId : cellsToFill) {
    for (size_t i = 0; i < particlesPerCell; ++i) {
      auto iAsDouble = static_cast<double>(i);
      Molecule m({iAsDouble, iAsDouble, iAsDouble}, {0, 0, 0}, numParticlesAdded++, 0);
      cells[cellId].addParticle(m);
    }
  }
  return cells;
}

void ParticleIteratorTest::testAllParticlesFoundPattern(const std::vector<size_t> &cellsWithParticles,
                                                        size_t numThreads, size_t numAdditionalParticleVectors) {
  constexpr size_t particlesPerCell = 4;
  auto cells = generateCellsWithPattern(10, cellsWithParticles, particlesPerCell);
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
             &cells, 0, nullptr, IteratorBehavior::haloAndOwned,
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
  auto data = generateCellsWithPattern(10, cellsToFill, numParticlesToAddPerCell);

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
                         ::testing::Combine(::testing::ValuesIn(threadNumsToTest), ::testing::Values(0, 1, 2, 4)));

template <class Iter>
auto ParticleIteratorTest::iteratorsBehaveEqually(Iter &iter1, Iter &iter2) {
  for (; iter1.isValid() and iter2.isValid(); ++iter1, ++iter2) {
    EXPECT_EQ(*iter1, *iter2);
  }
}

TEST_F(ParticleIteratorTest, testCopyConstructor) {
  constexpr size_t particlesPerCell = 4;
  auto cells = generateCellsWithPattern(10, {1ul, 3ul, 4ul, 7ul, 9ul}, particlesPerCell);

  std::vector<std::vector<Molecule>> additionalVectors(3);
  additionalVectors[1].emplace_back(Molecule({0., 0., 0.}, {0., 0., 0.}, 1337));

  constexpr bool modifyable = true;
  autopas::internal::ParticleIterator<Molecule, FMCell, modifyable> iter(
      &cells, 0, nullptr, IteratorBehavior::haloAndOwned, &additionalVectors);

  auto iterCopy{iter};
  iteratorsBehaveEqually(iter, iterCopy);
}

TEST_F(ParticleIteratorTest, testCopyAssignment) {
  constexpr size_t particlesPerCell = 4;
  auto cells = generateCellsWithPattern(10, {1ul, 3ul, 4ul, 7ul, 9ul}, particlesPerCell);

  std::vector<std::vector<Molecule>> additionalVectors(3);
  additionalVectors[1].emplace_back(Molecule({0., 0., 0.}, {0., 0., 0.}, 1337));

  constexpr bool modifyable = true;
  autopas::internal::ParticleIterator<Molecule, FMCell, modifyable> iter(
      &cells, 0, nullptr, IteratorBehavior::haloAndOwned, &additionalVectors);

  auto iterCopy = iter;
  iteratorsBehaveEqually(iter, iterCopy);
}