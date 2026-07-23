/**
 * @file NeighborIdentificationFunctorTest.cpp
 * @author S. Newcome
 * @date 22.07.2026
 */

#include <gtest/gtest.h>

#include <vector>

#include "autopas/baseFunctors/InteractionListGeneratorFunctor.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/AlignedAllocator.h"
#include "neighborIdentification/NeighborIdentificationFunctor.h"

/**
 * This is a simple test suite that checks that the public-facing NeighborIdentificationFunctor matches the internal
 * InteractionListGeneratorFunctor, except for the functor name.
 */

namespace {

using InternalFunctor = autopas::InteractionListGeneratorFunctor<autopas::ParticleBaseFP64, false>;
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
Functor_T::NeighborListAoSType runAoS(Cell &cell) {
  typename Functor_T::NeighborListAoSType map;
  for (auto &particle : cell) {
    map[&particle];
  }
  Functor_T functor(map, interactionLength, /*gatherNewton3Lists*/ false);
  functor.AoSFunctor(cell[0], cell[1], /*newton3*/ true);
  return map;
}

/**
 * Runs SoAFunctorSingle over the whole cell and returns the resulting neighbor lists.
 * @tparam Functor_T
 * @param cell
 */
template <class Functor_T>
Functor_T::NeighborListAoSType runSoASingle(Cell &cell) {
  typename Functor_T::NeighborListAoSType map;
  for (auto &particle : cell) {
    map[&particle];
  }
  Functor_T functor(map, interactionLength, /*gatherNewton3Lists*/ false);
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
Functor_T::NeighborListAoSType runSoAPair(Cell &cell1, Cell &cell2) {
  typename Functor_T::NeighborListAoSType map;
  for (auto &particle : cell1) {
    map[&particle];
  }
  for (auto &particle : cell2) {
    map[&particle];
  }
  Functor_T functor(map, interactionLength, /*gatherNewton3Lists*/ false);
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
Functor_T::NeighborListAoSType runSoAVerlet(Cell &cell) {
  typename Functor_T::NeighborListAoSType map;
  for (auto &particle : cell) {
    map[&particle];
  }
  Functor_T functor(map, interactionLength, /*gatherNewton3Lists*/ false);
  functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);
  const VerletListSoA neighborsOfZero{1};
  functor.SoAFunctorVerlet(cell._particleSoABuffer, 0, neighborsOfZero, /*newton3*/ true);
  return map;
}

/**
 * The public functor reports its own name rather than the internal functor's.
 */
TEST(NeighborIdentificationFunctor_Test, ReportsPublicName) {
  NeighborIdFunctor::NeighborListAoSType map;
  NeighborIdFunctor functor(map, interactionLength, /*gatherNewton3Lists*/ false);
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

}  // namespace