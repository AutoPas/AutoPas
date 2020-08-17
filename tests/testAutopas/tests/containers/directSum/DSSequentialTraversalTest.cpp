/**
 * @file DSSequentialTraversalTest.cpp
 * @author F. Gratl
 * @date 11/23/18
 */

#include "DSSequentialTraversalTest.h"

#include "autopas/containers/directSum/traversals/DSSequentialTraversal.h"
#include "autopasTools/generators/RandomGenerator.h"

using ::testing::_;
using ::testing::AtLeast;

TEST_F(DSSequentialTraversalTest, testTraversalAoS) { testTraversal(false); }

TEST_F(DSSequentialTraversalTest, testTraversalSoA) { testTraversal(true); }

void DSSequentialTraversalTest::testTraversal(bool useSoA) {
  size_t numParticles = 20;
  size_t numHaloParticles = 10;

  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(2);
  autopas::Particle defaultParticle;

  Particle particle;
  for (size_t i = 0; i < numParticles + numHaloParticles; ++i) {
    particle.setID(i);
    // first particles go in domain cell rest to halo cell
    if (i < numParticles) {
      particle.setR(autopasTools::generators::RandomGenerator::randomPosition({0, 0, 0}, {10, 10, 10}));
      cells[0].addParticle(particle);
    } else {
      particle.setR(autopasTools::generators::RandomGenerator::randomPosition({10, 10, 10}, {20, 20, 20}));
      cells[1].addParticle(particle);
    }
  }

  if (useSoA) {
    autopas::DSSequentialTraversal<FPCell, MFunctor, autopas::DataLayoutOption::soa, true> traversal(
        &functor, std::numeric_limits<double>::max());
    // domain SoA with itself
    EXPECT_CALL(functor, SoAFunctorSingle(_, true)).Times(1);
    // domain SoA with halo
    EXPECT_CALL(functor, SoAFunctorPair(_, _, true)).Times(1);
    std::for_each(cells.begin(), cells.end(), [](auto &c) { c._particleSoABuffer.resizeArrays(2); });
    traversal.setCellsToTraverse(cells);
    traversal.traverseParticlePairs();
  } else {
    autopas::DSSequentialTraversal<FPCell, MFunctor, autopas::DataLayoutOption::aos, true> traversal(
        &functor, std::numeric_limits<double>::max());
    // interactions in main cell + interactions with halo.
    size_t expectedFunctorCalls = numParticles * (numParticles - 1) / 2 + numParticles * numHaloParticles;
    EXPECT_CALL(functor, AoSFunctor(_, _, true)).Times((int)expectedFunctorCalls);
    traversal.setCellsToTraverse(cells);
    traversal.traverseParticlePairs();
  }
}

#if defined(AUTOPAS_CUDA)
TEST_F(DSSequentialTraversalTest, testTraversalCuda) {
  size_t numParticles = 20;
  size_t numHaloParticles = 10;

  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(2);
  autopas::Particle defaultParticle;

  Particle particle;
  for (size_t i = 0; i < numParticles + numHaloParticles; ++i) {
    particle.setID(i);
    // first particles go in domain cell rest to halo cell
    if (i < numParticles) {
      particle.setR(autopasTools::generators::RandomGenerator::randomPosition({0, 0, 0}, {10, 10, 10}));
      cells[0].addParticle(particle);
    } else {
      particle.setR(autopasTools::generators::RandomGenerator::randomPosition({10, 10, 10}, {20, 20, 20}));
      cells[1].addParticle(particle);
    }
  }

  autopas::DirectSumTraversal<FPCell, MFunctor, autopas::DataLayoutOption::cuda, true> traversal(
      &functor, 100 /*big enough to contain both particles*/);
  // domain SoA with itself
  EXPECT_CALL(functor, CudaFunctor(_, true)).Times(1);
  // domain SoA with halo
  EXPECT_CALL(functor, CudaFunctor(_, _, true)).Times(1);
  traversal.setCellsToTraverse(cells);
  traversal.traverseParticlePairs();
}

#endif
