/**
 * @file NeighborIdentificationFunctorTest.cpp
 * @author S. Newcome
 * @date 22.07.2026
 */

#include <gtest/gtest.h>

#include <vector>

#include "autopas/AutoPas.h"
#include "autopas/baseFunctors/InteractionListGeneratorFunctor.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/AlignedAllocator.h"
#include "neighborIdentificationLibrary/NeighborIdentificationFunctor.h"
#include "testingHelpers/NumThreadGuard.h"

/**
 * This is a simple test suite that checks that the public-facing NeighborIdentificationFunctor matches the internal
 * InteractionListGeneratorFunctor, except for the functor name.
 */

namespace {
using NeighborListAoSType = std::unordered_map<autopas::ParticleBaseFP64 *, std::vector<autopas::ParticleBaseFP64 *>>;
using AoSPolicy = autopas::AoSNeighborListPolicy<autopas::ParticleBaseFP64>;
using InternalFunctor = autopas::InteractionListGeneratorFunctor<autopas::ParticleBaseFP64, AoSPolicy>;
using NeighborIdFunctor = autopas::NeighborIdentificationFunctor<autopas::ParticleBaseFP64>;
using Cell = autopas::FullParticleCell<autopas::ParticleBaseFP64>;
using VerletListSoA = std::vector<size_t, autopas::AlignedAllocator<size_t>>;

/**
 * Interaction length used by all tests. Particles are placed well within it.
 */
constexpr double interactionLength = 1.0;

/**
 * Runs AoSFunctor on the first two particles of the cell and returns the resulting neighbor lists.
 * @tparam Functor_T
 * @param cell
 */
template <class Functor_T>
NeighborListAoSType runAoS(Cell &cell) {
  NeighborListAoSType map;
  for (auto &particle : cell) {
    map[&particle];
  }
  AoSPolicy policy(map);
  Functor_T functor(policy, interactionLength, /*gatherNewton3Lists*/ false);
  functor.AoSFunctor(cell[0], cell[1], /*newton3*/ true);
  return map;
}

/**
 * Runs SoAFunctorSingle over the whole cell and returns the resulting neighbor lists.
 * @tparam Functor_T
 * @param cell
 */
template <class Functor_T>
NeighborListAoSType runSoASingle(Cell &cell) {
  NeighborListAoSType map;
  for (auto &particle : cell) {
    map[&particle];
  }
  AoSPolicy policy(map);
  Functor_T functor(policy, interactionLength, /*gatherNewton3Lists*/ false);
  functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functor.SoAFunctorSingle(cell._particleSoABuffer, /*newton3*/ true);
  return map;
}

/**
 * Runs SoAFunctorPair over the two cells and returns the resulting neighbor lists.
 * @tparam Functor_T
 * @param cell1
 * @param cell2
 */
template <class Functor_T>
NeighborListAoSType runSoAPair(Cell &cell1, Cell &cell2) {
  NeighborListAoSType map;
  for (auto &particle : cell1) {
    map[&particle];
  }
  for (auto &particle : cell2) {
    map[&particle];
  }
  AoSPolicy policy(map);
  Functor_T functor(policy, interactionLength, /*gatherNewton3Lists*/ false);
  functor.SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, /*newton3*/ true);
  return map;
}

/**
 * Runs SoAFunctorVerlet for the first particle against the second and returns the resulting neighbor lists.
 * @tparam Functor_T
 * @param cell
 */
template <class Functor_T>
NeighborListAoSType runSoAVerlet(Cell &cell) {
  NeighborListAoSType map;
  for (auto &particle : cell) {
    map[&particle];
  }
  AoSPolicy policy(map);
  Functor_T functor(policy, interactionLength, /*gatherNewton3Lists*/ false);
  functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);
  const VerletListSoA neighborsOfZero{1};
  functor.SoAFunctorVerlet(cell._particleSoABuffer, 0, neighborsOfZero, /*newton3*/ true);
  return map;
}

/**
 * The public functor reports its own name rather than the internal functor's.
 */
TEST(NeighborIdentificationFunctor_Test, ReportsPublicName) {
  NeighborListAoSType map;
  AoSPolicy policy(map);
  NeighborIdFunctor functor(policy, interactionLength, /*gatherNewton3Lists*/ false);
  EXPECT_EQ(functor.getName(), "NeighborIdentificationFunctor");
}

/**
 * AoSFunctor produces the same neighbor lists through the public subclass as through the internal functor.
 */
TEST(NeighborIdentificationFunctor_Test, AoSFunctorMatchesInternal) {
  Cell cell;
  cell.addParticle(autopas::ParticleBaseFP64({0., 0., 0.}, {0., 0., 0.}, 0));
  cell.addParticle(autopas::ParticleBaseFP64({0.5, 0., 0.}, {0., 0., 0.}, 1));
  EXPECT_EQ(runAoS<InternalFunctor>(cell), runAoS<NeighborIdFunctor>(cell));
}

/**
 * SoAFunctorSingle produces the same neighbor lists through the public subclass as through the internal functor.
 */
TEST(NeighborIdentificationFunctor_Test, SoAFunctorSingleMatchesInternal) {
  Cell cell;
  cell.addParticle(autopas::ParticleBaseFP64({0., 0., 0.}, {0., 0., 0.}, 0));
  cell.addParticle(autopas::ParticleBaseFP64({0.5, 0., 0.}, {0., 0., 0.}, 1));
  EXPECT_EQ(runSoASingle<InternalFunctor>(cell), runSoASingle<NeighborIdFunctor>(cell));
}

/**
 * SoAFunctorPair produces the same neighbor lists through the public subclass as through the internal functor.
 */
TEST(NeighborIdentificationFunctor_Test, SoAFunctorPairMatchesInternal) {
  Cell cell1, cell2;
  cell1.addParticle(autopas::ParticleBaseFP64({0., 0., 0.}, {0., 0., 0.}, 0));
  cell2.addParticle(autopas::ParticleBaseFP64({0.5, 0., 0.}, {0., 0., 0.}, 1));
  EXPECT_EQ((runSoAPair<InternalFunctor>(cell1, cell2)), (runSoAPair<NeighborIdFunctor>(cell1, cell2)));
}

/**
 * SoAFunctorVerlet produces the same neighbor lists through the public subclass as through the internal functor.
 */
TEST(NeighborIdentificationFunctor_Test, SoAFunctorVerletMatchesInternal) {
  Cell cell;
  cell.addParticle(autopas::ParticleBaseFP64({0., 0., 0.}, {0., 0., 0.}, 0));
  cell.addParticle(autopas::ParticleBaseFP64({0.5, 0., 0.}, {0., 0., 0.}, 1));
  EXPECT_EQ(runSoAVerlet<InternalFunctor>(cell), runSoAVerlet<NeighborIdFunctor>(cell));
}

// ---------------------------------------------------------------------------------------------------------------------
// Neighbor List initialization test
// ---------------------------------------------------------------------------------------------------------------------

/**
 * Tests that NeighborIdentificationFunctor::initializeNeighborList delegates to the policy.
 */
TEST(NeighborIdentificationFunctor_Test, InitializeNeighborListSeedsAndClearsMap) {
  NeighborListAoSType neighborLists;
  AoSPolicy policy(neighborLists);
  NeighborIdFunctor functor(policy, interactionLength, /*gatherNewton3Lists*/ false);

  autopas::AutoPas<autopas::ParticleBaseFP64> autoPas;
  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1.0);
  autoPas.setVerletSkin(0.2);
  autoPas.init();

  autoPas.addParticle(autopas::ParticleBaseFP64({0., 0., 0.}, {0., 0., 0.}, 0));
  autoPas.addParticle(autopas::ParticleBaseFP64({0., 0., 0.}, {0., 0., 0.}, 1));
  autoPas.addParticle(autopas::ParticleBaseFP64({0., 0., 0.}, {0., 0., 0.}, 2));

  functor.initializeNeighborList(autoPas.begin());

  // One empty list per particle, keyed on the particle's address.
  size_t numParticles = 0;
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(neighborLists.count(&(*iter)), 1u)
        << "Particle with ID " << (*iter).getID() << " should have exactly one entry in the map, but found "
        << neighborLists.count(&(*iter)) << ".";
    EXPECT_TRUE(neighborLists.at(&(*iter)).empty());
    ++numParticles;
  }
  EXPECT_EQ(numParticles, 3u);
  EXPECT_EQ(neighborLists.size(), numParticles);

  // Populate a list and insert a stale key, then reinitialize: the second call must clear both.
  autopas::ParticleBaseFP64 staleParticle({0., 0., 0.}, {0., 0., 0.}, 99);
  neighborLists[&staleParticle];
  neighborLists.at(&(*autoPas.begin())).push_back(&staleParticle);

  functor.initializeNeighborList(autoPas.begin());

  EXPECT_EQ(neighborLists.count(&staleParticle), 0u);
  EXPECT_EQ(neighborLists.size(), numParticles);
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_TRUE(neighborLists.at(&(*iter)).empty());
  }
}

// ---------------------------------------------------------------------------------------------------------------------
// End-to-end test: Actual compilation and correct usage with AutoPas
// ---------------------------------------------------------------------------------------------------------------------

/**
 * End-to-end test that AutoPas compiles and runs with the NeighborIdentificationFunctor and that
 * computeInteractions() automatically initializes the functor's neighbor list map for the current particles, so the
 * caller does not have to seed it.
 *
 * This also acts as a regression test for a problem which arose during the development of the functor: between the
 * seeding of the map and the actual application of the functor, no particle is allowed to be moved in memory, or
 * the seeding becomes invalid. Originally, the balancing of the buffer vectors violated this, which had to be
 * corrected, hence why this test adds to the buffer.
 */
TEST(NeighborIdentificationFunctor_Test, ComputeInteractionsInitializesMap) {
  // Two threads => two per-thread buffers that get balanced.
  const NumThreadGuard numThreadGuard(2);

  autopas::AutoPas<autopas::ParticleBaseFP64> autoPas;
  autoPas.setBoxMin({0., 0., 0.});
  autoPas.setBoxMax({10., 10., 10.});
  autoPas.setCutoff(1.0);
  autoPas.setVerletSkin(0.2);
  autoPas.init();

  NeighborListAoSType neighborLists;
  AoSPolicy policy(neighborLists);
  NeighborIdFunctor functor(policy, /*interactionLength*/ 1.0, /*gatherNewton3Lists*/ false);

  // Two particles added into the container.
  autoPas.addParticle(autopas::ParticleBaseFP64({1.0, 1.0, 1.0}, {0., 0., 0.}, 0));
  autoPas.addParticle(autopas::ParticleBaseFP64({1.5, 1.0, 1.0}, {0., 0., 0.}, 1));

  // Test compute interactions
  autoPas.computeInteractions(&functor);

  // The two container particles must each have ended up in the other's list.
  {
    std::vector<autopas::ParticleBaseFP64 *> particlePtrs;
    for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
      particlePtrs.push_back(&(*iter));
    }
    ASSERT_EQ(particlePtrs.size(), 2u);
    ASSERT_EQ(neighborLists.size(), 2u);
    const auto &list0 = neighborLists.at(particlePtrs[0]);
    const auto &list1 = neighborLists.at(particlePtrs[1]);
    EXPECT_NE(std::ranges::find(list0, particlePtrs[1]), list0.end());
    EXPECT_NE(std::ranges::find(list1, particlePtrs[0]), list1.end());
  }

  // Two further owned particles (inside the box) and two halo particles (outside the box) now land in the buffers.
  // Specifically, they get added to the thread 0 buffer.
  autoPas.addParticle(autopas::ParticleBaseFP64({8.0, 8.0, 8.0}, {0., 0., 0.}, 2));
  autoPas.addParticle(autopas::ParticleBaseFP64({8.5, 8.0, 8.0}, {0., 0., 0.}, 3));
  autoPas.addHaloParticle(autopas::ParticleBaseFP64({-0.5, 1.0, 1.0}, {0., 0., 0.}, 4));
  autoPas.addHaloParticle(autopas::ParticleBaseFP64({10.5, 8.0, 8.0}, {0., 0., 0.}, 5));

  // Second pass balances the buffers (relocating buffer particles) and re-initializes the map. Must not crash, throw,
  // or leave dangling references.
  autoPas.computeInteractions(&functor);

  // Every current particle (container + buffers, owned + halo) must have an entry in the freshly initialized map.
  size_t numParticles = 0;
  for (auto iter = autoPas.begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(neighborLists.count(&(*iter)), 1u) << "Particle missing from the initialized neighbor list map.";
    ++numParticles;
  }
  EXPECT_EQ(numParticles, 6u);
  EXPECT_EQ(neighborLists.size(), 6u);
}

}  // namespace