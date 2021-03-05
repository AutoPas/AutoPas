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

/**
 * Test the iterator behavior for owning and halo molecules.
 * this function expects that two molecules are already added to container:
 * mol should be added as normal particle.
 * haloMol as haloParticle.
 * @tparam Container
 * @tparam Molecule
 * @param container should have already an added mol (as owning molecule) and
 * haloMol (als halo molecule)
 * @param mol
 * @param haloMol
 */
template <class Container, class Molecule>
void testContainerIteratorBehavior(Container &container, Molecule &mol, Molecule &haloMol) {
  // default
  int count = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : count)
#endif
  for (auto iter = container.begin(); iter.isValid(); ++iter) {
    count++;
  }
  EXPECT_EQ(count, 2);

  // haloAndOwned (same as default)
  count = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : count)
#endif
  for (auto iter = container.begin(IteratorBehavior::haloAndOwned); iter.isValid(); ++iter) {
    count++;
  }
  EXPECT_EQ(count, 2);

  // owned only
  count = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : count)
#endif
  for (auto iter = container.begin(IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    count++;
    EXPECT_EQ(iter->getID(), mol.getID());
  }
  EXPECT_EQ(count, 1);

  // halo only
  count = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : count)
#endif
  for (auto iter = container.begin(IteratorBehavior::haloOnly); iter.isValid(); ++iter) {
    count++;
    EXPECT_EQ(iter->getID(), haloMol.getID());
  }
  EXPECT_EQ(count, 1);
}

TEST_F(ParticleIteratorTest, testIteratorBehaviorDirectSum) {
  DirectSum<Molecule> ds({0., 0., 0.}, {10., 10., 10.}, 3, 0.);
  Molecule mol({1., 1., 1.}, {0., 0., 0.}, 1);
  ds.addParticle(mol);
  Molecule haloMol({-1., 1., 1.}, {0., 0., 0.}, 2, 0);
  ds.addHaloParticle(haloMol);

  testContainerIteratorBehavior(ds, mol, haloMol);
}

TEST_F(ParticleIteratorTest, testIteratorBehaviorLinkedCells) {
  LinkedCells<Molecule> linkedCells({0., 0., 0.}, {10., 10., 10.}, 3, 0., 1.);
  Molecule mol({1., 1., 1.}, {0., 0., 0.}, 1);
  linkedCells.addParticle(mol);
  Molecule haloMol({-1., 1., 1.}, {0., 0., 0.}, 2, 0);
  linkedCells.addHaloParticle(haloMol);

  testContainerIteratorBehavior(linkedCells, mol, haloMol);
}

TEST_F(ParticleIteratorTest, testIteratorBehaviorVerletLists) {
  VerletLists<Molecule> verletLists({0., 0., 0.}, {10., 10., 10.}, 3, 0.);
  Molecule mol({1., 1., 1.}, {0., 0., 0.}, 1, 0);
  verletLists.addParticle(mol);
  Molecule haloMol({-1., 1., 1.}, {0., 0., 0.}, 2, 0);
  verletLists.addHaloParticle(haloMol);

  // test normally
  testContainerIteratorBehavior(verletLists, mol, haloMol);

  // swap everything around, test if it still valid :)
  haloMol.setR({1., 1., 1.});
  mol.setR({-1., 1., 1.});
  verletLists.begin(IteratorBehavior::ownedOnly)->setR({-1., 1., 1.});
  verletLists.begin(IteratorBehavior::haloOnly)->setR({1., 1., 1.});

  testContainerIteratorBehavior(verletLists, mol, haloMol);
}

#ifdef AUTOPAS_OPENMP
const std::vector<size_t> threadNumsToTest{1, 2, 4};
#else
// no need to test more thread counts when openmp is not enabled
const std::vector<size_t> threadNumsToTest{1};
#endif

INSTANTIATE_TEST_SUITE_P(Generated, ParticleIteratorTest,
                         ::testing::Combine(::testing::ValuesIn(threadNumsToTest), ::testing::Values(0, 1, 2, 4)));