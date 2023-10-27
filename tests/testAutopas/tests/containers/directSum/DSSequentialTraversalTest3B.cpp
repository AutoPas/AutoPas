/**
 * @file DSSequentialTraversalTest3B.cpp
 * @author muehlhaeusser
 * @date 26.10.23
 */

#include "DSSequentialTraversalTest3B.h"

#include "autopas/containers/directSum/traversals/DSSequentialTraversal3B.h"
#include "autopasTools/generators/RandomGenerator.h"

using ::testing::_;
using ::testing::AtLeast;

TEST_F(DSSequentialTraversalTest3B, testTraversalAoS) { testTraversal(false); }

TEST_F(DSSequentialTraversalTest3B, testTraversalSoA) { testTraversal(true); }

void DSSequentialTraversalTest3B::testTraversal(bool useSoA) {
  size_t numParticles = 10;
  size_t numHaloParticles = 5;

  MTriwiseFunctor functor;
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
    autopas::DSSequentialTraversal3B<FPCell, MTriwiseFunctor, autopas::DataLayoutOption::soa, true> traversal(
        &functor, std::numeric_limits<double>::max());
    // domain SoA with itself
    EXPECT_CALL(functor, SoAFunctorSingle(_, true)).Times(1);
    // domain SoA with halo
    EXPECT_CALL(functor, SoAFunctorPair(_, _, true)).Times(1);
    std::for_each(cells.begin(), cells.end(), [](auto &c) { c._particleSoABuffer.resizeArrays(2); });
    traversal.setCellsToTraverse(cells);
    traversal.traverseParticleTriplets();
  } else {
    autopas::DSSequentialTraversal3B<FPCell, MTriwiseFunctor, autopas::DataLayoutOption::aos, true> traversal(
        &functor, std::numeric_limits<double>::max());
    // interactions in main cell + interactions with halo.
    size_t expectedFunctorCalls = numParticles * (numParticles - 1) * (numParticles - 2) / 6 + numParticles * (numParticles - 1) * numHaloParticles / 2 + numParticles * numHaloParticles * (numHaloParticles - 1) / 2;
    EXPECT_CALL(functor, AoSFunctor(_, _, _, true)).Times((int)expectedFunctorCalls);
    traversal.setCellsToTraverse(cells);
    traversal.traverseParticleTriplets();
  }
}
