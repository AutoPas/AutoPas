/**
 * @file DSSequentialTraversalTest.cpp
 * @author F. Gratl
 * @date 11/23/18
 */

#include "DSSequentialTraversalTest.h"

#include "autopas/containers/directSum/traversals/DSSequentialTraversal.h"
#include "autopasTools/generators/UniformGenerator.h"

using ::testing::_;
using ::testing::AtLeast;

TEST_F(DSSequentialTraversalTest, testTraversalPairwiseAoS) { testTraversalPairwise(false); }

TEST_F(DSSequentialTraversalTest, testTraversalPairwiseSoA) { testTraversalPairwise(true); }

TEST_F(DSSequentialTraversalTest, testTraversalTriwiseAoS) { testTraversalTriwise(false); }

TEST_F(DSSequentialTraversalTest, testTraversalTriwiseSoA) { testTraversalTriwise(true); }

std::vector<FPCell> DSSequentialTraversalTest::fillParticleCells(size_t numParticles, size_t numHaloParticlesPerCell) {
  autopas::ParticleBaseFP64 particle;
  std::vector<FPCell> cells(7);
  std::mt19937 generator(42);

  // helper function to randomly place particles in the specified cell
  auto addParticlesToCell = [&](const size_t cellID, const size_t num, const std::array<double, 3> boxMin,
                                const std::array<double, 3> boxMax) {
    for (size_t i = 0; i < num; i++) {
      particle.setR(autopasTools::generators::UniformGenerator::randomPosition(generator, boxMin, boxMax));
      cells[cellID].addParticle(particle);
    }
  };

  // Add particles to all cells
  addParticlesToCell(0, numParticles, {0, 0, 0}, {10, 10, 10});
  addParticlesToCell(1, numHaloParticlesPerCell, {-10, -10, -10}, {0, 20, 20});
  addParticlesToCell(2, numHaloParticlesPerCell, {10, -10, -10}, {20, 20, 20});
  addParticlesToCell(3, numHaloParticlesPerCell, {0, -10, -10}, {10, 0, 20});
  addParticlesToCell(4, numHaloParticlesPerCell, {0, 10, -10}, {10, 20, 20});
  addParticlesToCell(5, numHaloParticlesPerCell, {0, 0, -10}, {10, 10, 0});
  addParticlesToCell(6, numHaloParticlesPerCell, {0, 0, 10}, {10, 10, 20});

  return cells;
}

void DSSequentialTraversalTest::testTraversalPairwise(bool useSoA) {
  const size_t numParticles = 20;
  const size_t numHaloParticlesPerCell = 2;
  const size_t numHaloParticles = numHaloParticlesPerCell * 6;

  auto cells = fillParticleCells(numParticles, numHaloParticlesPerCell);

  MPairwiseFunctor functor;

  if (useSoA) {
    autopas::DSSequentialTraversal<FPCell, MPairwiseFunctor> traversal(&functor, std::numeric_limits<double>::max(),
                                                                       autopas::DataLayoutOption::soa, true);
    // domain SoA with itself
    EXPECT_CALL(functor, SoAFunctorSingle(_, true)).Times(1);
    // domain SoA with halo domains
    EXPECT_CALL(functor, SoAFunctorPair(_, _, true)).Times(6);
    std::for_each(cells.begin(), cells.end(), [](auto &c) { c._particleSoABuffer.resizeArrays(2); });
    traversal.setCellsToTraverse(cells);
    traversal.traverseParticles();
  } else {
    autopas::DSSequentialTraversal<FPCell, MPairwiseFunctor> traversal(&functor, std::numeric_limits<double>::max(),
                                                                       autopas::DataLayoutOption::aos, true);
    // interactions in main cell + interactions with halo cells.
    size_t expectedFunctorCalls = numParticles * (numParticles - 1) / 2 + numParticles * numHaloParticles;
    EXPECT_CALL(functor, AoSFunctor(_, _, true)).Times((int)expectedFunctorCalls);
    traversal.setCellsToTraverse(cells);
    traversal.traverseParticles();
  }
}

void DSSequentialTraversalTest::testTraversalTriwise(bool useSoA) {
  const size_t numParticles = 20;
  const size_t numHaloParticlesPerCell = 2;
  const size_t numHaloParticles = numHaloParticlesPerCell * 6;

  auto cells = fillParticleCells(numParticles, numHaloParticlesPerCell);

  MTriwiseFunctor functor;

  if (useSoA) {
    autopas::DSSequentialTraversal<FPCell, MTriwiseFunctor> traversal(&functor, std::numeric_limits<double>::max(),
                                                                      autopas::DataLayoutOption::soa, true);
    // domain SoA with itself
    EXPECT_CALL(functor, SoAFunctorSingle(_, true)).Times(1);
    // domain SoA with halo domains
    EXPECT_CALL(functor, SoAFunctorPair(_, _, true)).Times(6);
    // domain SoA with halo domains
    EXPECT_CALL(functor, SoAFunctorTriple(_, _, _, true)).Times(15);
    std::for_each(cells.begin(), cells.end(), [](auto &c) { c._particleSoABuffer.resizeArrays(2); });
    traversal.setCellsToTraverse(cells);
    traversal.traverseParticles();
  } else {
    autopas::DSSequentialTraversal<FPCell, MTriwiseFunctor> traversal(&functor, std::numeric_limits<double>::max(),
                                                                      autopas::DataLayoutOption::aos, true);
    // interactions in main cell + interactions with halo cells.
    size_t expectedFunctorCalls = numParticles * (numParticles - 1) * (numParticles - 2) / 6 +
                                  numParticles * (numParticles - 1) * numHaloParticles / 2 +
                                  numParticles * numHaloParticles * (numHaloParticles - 1) / 2;
    EXPECT_CALL(functor, AoSFunctor(_, _, _, true)).Times((int)expectedFunctorCalls);
    traversal.setCellsToTraverse(cells);
    traversal.traverseParticles();
  }
}