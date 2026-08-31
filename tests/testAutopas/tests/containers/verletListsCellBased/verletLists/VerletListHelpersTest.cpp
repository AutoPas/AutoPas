/**
 * @file VerletListHelpersTest.cpp
 * @author muehlhaeusser
 * @date 2026-07-26
 */

#include "VerletListHelpersTest.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/utils/WrapOpenMP.h"
#include "testingHelpers/NumThreadGuard.h"

using namespace autopas;

/**
 * Tests that the counter functor behaves safely across multiple threads without data races
 * or lost increments.
 */
TEST_F(VerletListHelpersTest, CounterFunctorThreadSafety) {
  constexpr size_t numParticles = 1000;
  generateDenseParticles(numParticles);

  // Initialize cache-line padded atomic counters
  std::vector<typename Helpers::VerletListCounterFunctor::PaddedAtomic> counts(numParticles);
  Helpers::VerletListCounterFunctor functor(counts, particleToIndex, interactionLength);

  {
    // Force execution across 4 threads
    NumThreadGuard numThreadGuard(4);
        AUTOPAS_OPENMP(parallel for schedule(dynamic))
        for (size_t i = 0; i < numParticles; ++i) {
          for (size_t j = i + 1; j < numParticles; ++j) {
            // newton3 = true means particle `i` stores `j`, but `j` does not store `i`
            functor.AoSFunctor(cell[i], cell[j], true);
          }
        }
  }

  // Validation: particle i should have interacted with exactly (numParticles - 1 - i) particles
  for (size_t i = 0; i < numParticles; ++i) {
    EXPECT_EQ(counts[i].value.load(std::memory_order_seq_cst), numParticles - 1 - i)
        << "Lost increments detected for particle " << i << " due to race condition!";
  }
}

/**
 * Tests that the filler functor can safely write indices directly into a pre-allocated CRS
 * array across multiple threads without overwriting slots or throwing out-of-bounds.
 */
TEST_F(VerletListHelpersTest, FillerFunctorThreadSafety) {
  constexpr size_t numParticles = 1000;
  generateDenseParticles(numParticles);

  typename Helpers::NeighborListCRS neighborList;

  // Simulate Prefix Sum
  neighborList.offsets.resize(numParticles + 1, 0);
  for (size_t i = 0; i < numParticles; ++i) {
    neighborList.offsets[i + 1] = neighborList.offsets[i] + (numParticles - 1 - i);
  }
  neighborList.indices.resize(neighborList.offsets.back());

  // Initialize fill cursors to the beginning of each particle's offset
  std::vector<typename Helpers::VerletListCounterFunctor::PaddedAtomic> fillPos(numParticles);
  for (size_t i = 0; i < numParticles; ++i) {
    fillPos[i].value.store(neighborList.offsets[i], std::memory_order_relaxed);
  }

  Helpers::VerletListFillerFunctor filler(neighborList, fillPos, particleToIndex, interactionLength);

  {
    // Force execution across 4 threads
    NumThreadGuard numThreadGuard(4);
        AUTOPAS_OPENMP(parallel for schedule(dynamic))
        for (size_t i = 0; i < numParticles; ++i) {
          for (size_t j = i + 1; j < numParticles; ++j) {
            filler.AoSFunctor(cell[i], cell[j], true);
          }
        }
  }

  // Validation
  for (size_t i = 0; i < numParticles; ++i) {
    auto startIter = neighborList.indices.begin() + neighborList.offsets[i];
    auto endIter = neighborList.indices.begin() + neighborList.offsets[i + 1];

    // Copy and sort the stored neighbor indices (fetch_add insertion order is non-deterministic in parallel)
    std::vector<size_t> neighbors(startIter, endIter);
    std::sort(neighbors.begin(), neighbors.end());

    ASSERT_EQ(neighbors.size(), numParticles - 1 - i);
    for (size_t j = 0; j < neighbors.size(); ++j) {
      EXPECT_EQ(neighbors[j], i + 1 + j) << "Corrupted CRS array insertion detected!";
    }
  }
}

/**
     * Tests that when newton3=false, the counter functor increments the neighbor count for
  BOTH particles.
     */
TEST_F(VerletListHelpersTest, CounterFunctorNewton3False) {
  cell.addParticle(ParticleType({0., 0., 0.}, {0., 0., 0.}, 0));
  cell.addParticle(ParticleType({0., 0., 0.}, {0., 0., 0.}, 1));
  particleToIndex[&cell[0]] = 0;
  particleToIndex[&cell[1]] = 1;

  std::vector<typename Helpers::VerletListCounterFunctor::PaddedAtomic> counts(2);
  Helpers::VerletListCounterFunctor functor(counts, particleToIndex, interactionLength);

  // newton3 = false -> should register for both 0 and 1
  functor.AoSFunctor(cell[0], cell[1], false);

  EXPECT_EQ(counts[0].value.load(), 1);
  EXPECT_EQ(counts[1].value.load(), 1);
}

/**
 * Tests that when newton3=false, the filler functor inserts the neighbor index for
BOTH particles.
 */
TEST_F(VerletListHelpersTest, FillerFunctorNewton3False) {
  cell.addParticle(ParticleType({0., 0., 0.}, {0., 0., 0.}, 0));
  cell.addParticle(ParticleType({0., 0., 0.}, {0., 0., 0.}, 1));
  particleToIndex[&cell[0]] = 0;
  particleToIndex[&cell[1]] = 1;

  typename Helpers::NeighborListCRS neighborList;
  neighborList.offsets = {0, 1, 2};  // 1 neighbor each
  neighborList.indices.resize(2);

  std::vector<typename Helpers::VerletListCounterFunctor::PaddedAtomic> fillPos(2);
  fillPos[0].value.store(0);
  fillPos[1].value.store(1);

  Helpers::VerletListFillerFunctor filler(neighborList, fillPos, particleToIndex, interactionLength);
  filler.AoSFunctor(cell[0], cell[1], false);  // newton3 = false

  EXPECT_EQ(neighborList.indices[0], 1) << "Particle 0 should have particle 1 in its list.";
  EXPECT_EQ(neighborList.indices[1], 0) << "Particle 1 should have particle 0 in its list.";
}

/**
 * Tests that particles outside the interaction length are completely ignored.
 */
TEST_F(VerletListHelpersTest, DistanceFiltering) {
  // Particles placed further apart than interactionLength (1.0)
  cell.addParticle(ParticleType({0., 0., 0.}, {0., 0., 0.}, 0));
  cell.addParticle(ParticleType({2.0, 0., 0.}, {0., 0., 0.}, 1));
  particleToIndex[&cell[0]] = 0;
  particleToIndex[&cell[1]] = 1;

  std::vector<typename Helpers::VerletListCounterFunctor::PaddedAtomic> counts(2);
  Helpers::VerletListCounterFunctor functor(counts, particleToIndex, interactionLength);

  functor.AoSFunctor(cell[0], cell[1], true);

  // Neither should be counted
  EXPECT_EQ(counts[0].value.load(), 0);
  EXPECT_EQ(counts[1].value.load(), 0);
}

/**
 * Tests that Dummy particles are skipped and do not generate interactions.
 */
TEST_F(VerletListHelpersTest, DummyParticleFiltering) {
  // Two particles at origin (distance 0.0), but one is a dummy
  ParticleType p0({0., 0., 0.}, {0., 0., 0.}, 0);
  ParticleType p1({0., 0., 0.}, {0., 0., 0.}, 1);
  p1.setOwnershipState(autopas::OwnershipState::dummy);

  cell.addParticle(p0);
  cell.addParticle(p1);
  particleToIndex[&cell[0]] = 0;
  particleToIndex[&cell[1]] = 1;

  std::vector<typename Helpers::VerletListCounterFunctor::PaddedAtomic> counts(2);
  Helpers::VerletListCounterFunctor functor(counts, particleToIndex, interactionLength);

  functor.AoSFunctor(cell[0], cell[1], true);

  // Since particle 1 is a dummy, the interaction should be aborted entirely.
  EXPECT_EQ(counts[0].value.load(), 0);
}

/**
 * Tests that PairVerletListCounterFunctor behaves safely across multiple threads without data races.
 */
TEST_F(VerletListHelpersTest, PairCounterFunctorThreadSafety) {
  constexpr size_t numParticles = 50;
  generateDenseParticles(numParticles);

  std::vector<typename Helpers::VerletListCounterFunctor::PaddedAtomic> counts(numParticles);
  Helpers::PairVerletListCounterFunctor functor(counts, particleToIndex, interactionLength);

  {
    NumThreadGuard numThreadGuard(4);
    AUTOPAS_OPENMP(parallel for schedule(dynamic))
    for (size_t i = 0; i < numParticles; ++i) {
      for (size_t j = i + 1; j < numParticles; ++j) {
        for (size_t k = j + 1; k < numParticles; ++k) {
          functor.AoSFunctor(cell[i], cell[j], cell[k], false);
        }
      }
    }
  }

  for (size_t i = 0; i < numParticles; ++i) {
    const size_t remaining = numParticles - 1 - i;
    const size_t expectedPairs = remaining >= 2 ? (remaining * (remaining - 1)) / 2 : 0;
    EXPECT_EQ(counts[i].value.load(std::memory_order_seq_cst), expectedPairs)
        << "Lost increments detected for particle " << i << " due to race condition!";
  }
}

/**
 * Tests that PairVerletListFillerFunctor can safely write index pairs into a pre-allocated CRS array.
 */
TEST_F(VerletListHelpersTest, PairFillerFunctorThreadSafety) {
  constexpr size_t numParticles = 50;
  generateDenseParticles(numParticles);

  typename Helpers::NeighborPairsListCRS neighborPairsList;
  neighborPairsList.offsets.resize(numParticles + 1, 0);
  for (size_t i = 0; i < numParticles; ++i) {
    const size_t remaining = numParticles - 1 - i;
    const size_t count = remaining >= 2 ? (remaining * (remaining - 1)) / 2 : 0;
    neighborPairsList.offsets[i + 1] = neighborPairsList.offsets[i] + count;
  }
  neighborPairsList.pairs.resize(neighborPairsList.offsets.back());

  std::vector<typename Helpers::VerletListCounterFunctor::PaddedAtomic> fillPos(numParticles);
  for (size_t i = 0; i < numParticles; ++i) {
    fillPos[i].value.store(neighborPairsList.offsets[i], std::memory_order_relaxed);
  }

  Helpers::PairVerletListFillerFunctor filler(neighborPairsList, fillPos, particleToIndex, interactionLength);

  {
    NumThreadGuard numThreadGuard(4);
    AUTOPAS_OPENMP(parallel for schedule(dynamic))
    for (size_t i = 0; i < numParticles; ++i) {
      for (size_t j = i + 1; j < numParticles; ++j) {
        for (size_t k = j + 1; k < numParticles; ++k) {
          filler.AoSFunctor(cell[i], cell[j], cell[k], false);
        }
      }
    }
  }

  for (size_t i = 0; i < numParticles; ++i) {
    const size_t remaining = numParticles - 1 - i;
    const size_t expectedPairs = remaining >= 2 ? (remaining * (remaining - 1)) / 2 : 0;
    EXPECT_EQ(neighborPairsList.count(i), expectedPairs);
  }
}