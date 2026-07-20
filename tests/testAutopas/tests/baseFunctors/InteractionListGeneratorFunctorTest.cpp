/**
 * @file InteractionListGeneratorFunctorTest.cpp
 * @author S. Newcome
 * @date 16.07.2026
 */

#include "InteractionListGeneratorFunctorTest.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "autopas/baseFunctors/InteractionListGeneratorFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/WrapOpenMP.h"
#include "generators/src/GridGenerator.h"
#include "testingHelpers/NumThreadGuard.h"
#include "testingHelpers/commonTypedefs.h"

namespace {

using autopas::OwnershipState;
using InteractionListGenFunc = autopas::InteractionListGeneratorFunctor<Molecule>;
using InteractionListTypeAoS = InteractionListGenFunc::NeighborListAoSType;
using VerletListTypeSoA = std::vector<size_t, autopas::AlignedAllocator<size_t>>;

/**
 * Interaction length used by all tests.
 */
constexpr double cutoffLength = 2.5;

/**
 * Number of threads used by the thread safety tests.
 */
constexpr int numThreads = 2;

/**
 * A few thread safety tests need a coloring scheme for n3 use for four "elements" (AoS: particles, SoA: cells). We
 * split the six pairs of elements into three rounds where two threads can work independently.
 */
using ColoringRound = std::array<std::pair<size_t, size_t>, numThreads>;
constexpr std::array<ColoringRound, 3> coloringRounds{ColoringRound{{{0, 1}, {2, 3}}}, ColoringRound{{{0, 2}, {1, 3}}},
                                                      ColoringRound{{{0, 3}, {1, 2}}}};

/**
 * Adds particles at the given positions to the cell. Particle IDs are set to idOffset + the particle's index so that
 * the generated neighbor lists can be identified by ID; idOffset keeps IDs unique when several cells are used
 * together.
 * @param cell The cell to add the particles to.
 * @param positions Position of each particle; one particle is created per entry.
 * @param ownerships OwnershipState for each particle. Must have the same size as positions.
 * @param idOffset Value added to each particle's index to form its ID. Defaults to 0.
 */
void fillCell(FMCell &cell, const std::vector<std::array<double, 3>> &positions,
              const std::vector<OwnershipState> &ownerships, const size_t idOffset = 0) {
  for (size_t i = 0; i < positions.size(); ++i) {
    Molecule m(positions[i], {0., 0., 0.}, i + idOffset);
    m.setOwnershipState(ownerships[i]);
    cell.addParticle(m);
  }
}

/**
 * Expands an ownership variant to one OwnershipState per particle by distributing the variant's states in a
 * round-robin.
 * @param variant The ownership variant to expand.
 * @param numParticles Number of particles to produce states for.
 * @return The OwnershipState of each particle.
 */
std::vector<OwnershipState> expandOwnerships(const OwnershipVariant &variant, size_t numParticles) {
  std::vector<OwnershipState> ownerships(numParticles);
  for (size_t i = 0; i < numParticles; ++i) {
    ownerships[i] = variant.states[i % variant.states.size()];
  }
  return ownerships;
}

/**
 * Pre-populates the neighbor list interactionList with an empty list for every particle in the cell (keys are
 * &cell[i]). This ensures the interactionList is never rehashed during the (potentially parallel) functor application,
 * and should mirror uses of this function. Must be called for each cell used in the test with the same interactionList.
 * @param interactionList The neighbor list map to populate.
 * @param cell The cell whose particles are used as keys.
 */
void initMap(InteractionListTypeAoS &interactionList, FMCell &cell) {
  for (size_t i = 0; i < cell.size(); ++i) {
    interactionList[&cell[i]];
  }
}

/**
 * Creates the four cells of 3x3x3 particles used by the thread safety tests. The cells are laid out next to each other
 * in the x-y plane and their particle IDs continue from cell to cell, so that every particle has a unique ID.
 * @return The created cells.
 */
std::vector<FMCell> makeGridCells() {
  constexpr size_t numCells = 4;
  constexpr std::array<size_t, 3> particlesPerDim{3, 3, 3};
  constexpr size_t particlesPerCell = particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2];

  std::vector<FMCell> cells(numCells);
  for (size_t cellID = 0; cellID < numCells; ++cellID) {
    const Molecule defaultParticle({0., 0., 0.}, {0., 0., 0.}, cellID * particlesPerCell);
    // Cells are placed in a 2x2 arrangement
    const size_t cellsX = cellID / 2;
    const size_t cellsY = cellID % 2;
    const std::array<double, 3> offset{static_cast<double>(3 * cellsX), static_cast<double>(3 * cellsY), 0.};
    autopasTools::generators::GridGenerator::fillWithParticles(cells[cellID], particlesPerDim, defaultParticle,
                                                               {1., 1., 1.}, offset);
  }
  return cells;
}

/**
 * Collects the IDs of the neighbors stored for a particle, sorted ascending, for comparison with the expected IDs.
 * Sorting makes the comparisons independent of the order in which neighbors were added.
 * @param map The neighbor list map to read from.
 * @param cell The cell owning the particle.
 * @param index Index of the particle within the cell.
 * @return The sorted IDs of the particles in &cell[index]'s neighbor list.
 */
std::vector<size_t> neighborIDs(const InteractionListTypeAoS &map, FMCell &cell, size_t index) {
  const auto &neighbors = map.at(&cell[index]);
  std::vector<size_t> ids;
  ids.reserve(neighbors.size());
  std::ranges::transform(neighbors, std::back_inserter(ids), [](const auto *neighbor) { return neighbor->getID(); });
  std::ranges::sort(ids);
  return ids;
}

/**
 * Rewrites the neighbor lists such that every stored pair sits in the list of its lower-ID particle. When gathering
 * newton3 lists the functor stores each pair in exactly one particle's list but makes no guarantee in which, so
 * gathered lists are sorted into this canonical form before comparison against expected lists written accordingly.
 * @param map The neighbor list map to rewrite.
 */
void sortIntoLowerIDLists(InteractionListTypeAoS &map) {
  for (auto &[particlePtr, neighborList] : map) {
    for (auto *neighborPtr : neighborList) {
      if (neighborPtr->getID() < particlePtr->getID()) {
        map.at(neighborPtr).push_back(particlePtr);
      }
    }
    std::erase_if(neighborList,
                  [id = particlePtr->getID()](const auto *neighborPtr) { return neighborPtr->getID() < id; });
  }
}

/**
 * Compares the neighbor list of every particle in the cell against the expected IDs.
 * @param map The neighbor list map to read from.
 * @param cell The cell whose particles' lists are checked.
 * @param expected Expected neighbor list of each particle, given as ascending particle IDs. One entry per particle.
 */
void compareWithExpectedLists(const InteractionListTypeAoS &map, FMCell &cell, const ExpectedLists &expected) {
  ASSERT_EQ(expected.size(), cell.size());
  for (size_t i = 0; i < cell.size(); ++i) {
    EXPECT_EQ(neighborIDs(map, cell, i), expected[i]) << "Wrong neighbor list for particle with ID " << cell[i].getID();
  }
}

/**
 * Compares the neighbor list of every particle in the cell against the list generated for the corresponding particle of
 * the reference cell. Used by the thread safety tests to compare a parallel application against a serial one.
 * @param map The neighbor list map produced by the parallel application.
 * @param cell The cell used by the parallel application.
 * @param referenceMap The neighbor list map produced by the serial application.
 * @param referenceCell The cell used by the serial application.
 */
void expectListsMatch(const InteractionListTypeAoS &map, FMCell &cell, const InteractionListTypeAoS &referenceMap,
                      FMCell &referenceCell) {
  for (size_t i = 0; i < cell.size(); ++i) {
    EXPECT_EQ(neighborIDs(map, cell, i), neighborIDs(referenceMap, referenceCell, i))
        << "Parallel and serial neighbor lists differ for particle with ID " << cell[i].getID();
  }
}

/**
 * Applies the AoSFunctor to every pair of particles in the cell, exactly as an AutoPas traversal would: if newton3 is
 * used, once per pair (i < j) with newton3 enabled, otherwise once for each direction of the pair with newton3
 * disabled.
 * @param functor The functor to apply.
 * @param cell The cell whose particle pairs are visited.
 * @param n3Used Whether newton3 is used.
 */
void applyAoSOverAllPairs(InteractionListGenFunc &functor, FMCell &cell, bool n3Used) {
  for (size_t i = 0; i < cell.size(); ++i) {
    for (size_t j = i + 1; j < cell.size(); ++j) {
      if (n3Used) {
        functor.AoSFunctor(cell[i], cell[j], /*newton3*/ true);
      } else {
        functor.AoSFunctor(cell[i], cell[j], /*newton3*/ false);
        functor.AoSFunctor(cell[j], cell[i], /*newton3*/ false);
      }
    }
  }
}

/**
 * Removes dummy particles from expected neighbor lists. Used to create lists of expected lists with dummies from lists
 * which assume all particles are owned.
 * @param expected Expected neighbor lists, assuming particles are owned/halo, indexed by particle ID.
 * @param ownerships OwnershipState of each particle.
 * @return The pruned lists.
 */
ExpectedLists pruneDummies(ExpectedLists expected, const std::vector<OwnershipState> &ownerships) {
  for (size_t i = 0; i < expected.size(); ++i) {
    if (ownerships[i] == OwnershipState::dummy) {
      expected[i].clear();
    } else {
      std::erase_if(expected[i], [&](size_t id) { return ownerships[id] == OwnershipState::dummy; });
    }
  }
  return expected;
}

/**
 * Reduces bidirectional expected neighbor lists to gathered newton3 lists. We keep the interaction where the particle
 * with the higher ID is in the list of the particle with the lower.
 * @param expected Expected bidirectional neighbor lists, indexed by particle ID.
 * @return The pruned lists.
 */
ExpectedLists pruneToGatheredN3Lists(ExpectedLists expected) {
  for (size_t i = 0; i < expected.size(); ++i) {
    std::erase_if(expected[i], [i](size_t id) { return id < i; });
  }
  return expected;
}

/**
 * A vector of particles used by most test cases. In-range pairs (0,1), (0,2), (1,2), (1,3), (2,3). Particle 3 is out of
 * range of particle 0, and particle 4 is isolated. No pair sits exactly at the cutoff.
 */
const std::vector<std::array<double, 3>> standardPositions{
    {0., 0., 0.}, {1., 1., 0.}, {1., 1., 2.}, {2., 3., 1.}, {10., 10., 10.}};

/**
 * A vector of particles such that no particle is in range of any other.
 */
const std::vector<std::array<double, 3>> spreadPositions{
    {0., 0., 0.}, {2., 2., 2.}, {4., 4., 4.}, {6., 6., 6.}, {8., 8., 8.}};

/**
 * Ownership variants used to build the test case cross products, distributed over the particles in a round-robin: all
 * owned, a mix with halos (which must participate like owned particles), a mix with dummies (whose pairs must all be
 * skipped), and all three states.
 */
const std::vector<OwnershipVariant> ownershipVariants{
    {"AllOwned", {OwnershipState::owned}},
    {"SomeHalo", {OwnershipState::halo, OwnershipState::owned}},
    {"SomeDummy", {OwnershipState::owned, OwnershipState::dummy}},
    {"Mixed", {OwnershipState::dummy, OwnershipState::halo, OwnershipState::owned}}};

/**
 * The three valid gatherNewton3Lists / newton3 combinations. gatherNewton3Lists=true requires newton3=true, so that
 * combination is not covered. Used by all test case cross products.
 */
const std::vector<N3Mode> n3Modes{{"N3Lists", /*gatherN3Lists*/ true, /*n3Used*/ true},
                                  {"NoN3ListsN3Generation", /*gatherN3Lists*/ false, /*n3Used*/ true},
                                  {"NoN3ListsNoN3Generation", /*gatherN3Lists*/ false, /*n3Used*/ false}};

/**
 * Expected neighbor lists for particles placed at the `standardPositions`.
 */
const ExpectedLists bidirectionalExpectedNeighborsStandard{{1, 2}, {0, 2, 3}, {0, 1, 3}, {1, 2}, {}};

/**
 * Expected neighbor lists for particles placed at the `spreadPositions`. Is simply a vector of empty lists.
 */
const ExpectedLists bidirectionalExpectedNeighborsOutOfRange(spreadPositions.size());

/**
 * Full Verlet lists of the particles placed at the `standardPositions`, given as particle indices. Assumes a Verlet
 * skin of 1.3. This is much larger than in practice, but means that pair (0,3) is within the Verlet radius whilst
 * outside the cutoff.
 */
const ExpectedLists verletListsStandard{{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}, {}};

/**
 * Full Verlet lists of the particles placed at the `spreadPositions`, given as particle indices. Built with the same
 * Verlet skin as `verletListsStandard`. No pair given below is within the cutoff, so all should be filtered out by the
 * functor.
 */
const ExpectedLists verletListsOutOfRange{{1}, {0, 2}, {1, 3}, {2, 4}, {3}};

/**
 * Builds the test cases shared by the AoSFunctor, SoAFunctorSingle and SoAFunctorVerlet tests as the cross product of
 * the particle positions, the newton3 modes and the ownership variants.
 * @return The generated test cases.
 */
std::vector<ListGenerationParams> getListGenerationParams() {
  std::vector<ListGenerationParams> params;
  for (const bool outOfRange : {false, true}) {
    const auto &positions = outOfRange ? spreadPositions : standardPositions;
    const auto &verletLists = outOfRange ? verletListsOutOfRange : verletListsStandard;
    for (const auto &mode : n3Modes) {
      for (const auto &ownershipVariant : ownershipVariants) {
        const auto ownerships = expandOwnerships(ownershipVariant, positions.size());
        auto expectedLists = pruneDummies(
            outOfRange ? bidirectionalExpectedNeighborsOutOfRange : bidirectionalExpectedNeighborsStandard, ownerships);
        if (mode.gatherN3Lists) {
          expectedLists = pruneToGatheredN3Lists(expectedLists);
        }
        // With newton3 used, each pair is handled by a single call, so every pair may only appear in one particle's
        // Verlet list. (This mimics the InteractionListGenFunc used to generate the Verlet Lists being called with
        // gatherN3Lists=true)
        const auto prunedVerletLists = mode.n3Used ? pruneToGatheredN3Lists(verletLists) : verletLists;

        params.push_back({.name = mode.name + ownershipVariant.name + (outOfRange ? "AllOutOfRange" : ""),
                          .positions = positions,
                          .ownerships = ownerships,
                          .n3Mode = mode,
                          .expectedNeighbors = expectedLists,
                          .verletLists = prunedVerletLists});
      }
    }
  }
  return params;
}

/**
 * Builds the SoAFunctorPair test cases as the cross product of the second cell's positions (in / out of range), the
 * newton3 modes and the ownership variants (applied to both cells' particles). Without newton3 only the first
 * cell's lists are written.
 * @return The generated test cases.
 */
std::vector<PairListGenerationParams> getPairListGenerationParams() {
  std::vector<PairListGenerationParams> params;

  const std::vector<std::array<double, 3>> cell1Positions{{0., 0., 0.}, {1., 1., 1.}, {0., 2., 0.}};
  const std::vector<std::array<double, 3>> inRangeCell2Positions{{2., 0., 0.}, {2., 2., 0.}, {5., 0., 0.}};
  const std::vector<std::array<double, 3>> outOfRangeCell2Positions{{10., 0., 0.}, {15., 0., 0.}, {0., 4., 4.}};

  // Cell 2's particles are given IDs continuing from cell 1's, so lists referring to them use those IDs.
  const size_t cell2IDOffset = cell1Positions.size();

  // Expected lists of the particles in cell 1's interactions from cell 2...
  const ExpectedLists baseExpectedListsCell1{{cell2IDOffset}, {cell2IDOffset, cell2IDOffset + 1}, {cell2IDOffset + 1}};
  // and vice-versa
  const ExpectedLists baseExpectedListsCell2{{0, 1}, {1, 2}, {}};

  // Empty lists, for the cases where a cell gets no interactions at all.
  const ExpectedLists noNeighborsCell1(cell1Positions.size());
  const ExpectedLists noNeighborsCell2(inRangeCell2Positions.size());

  for (const bool cell2OutOfRange : {false, true}) {
    for (const auto &mode : n3Modes) {
      for (const auto &ownershipVariant : ownershipVariants) {
        const auto ownerships =
            expandOwnerships(ownershipVariant, cell1Positions.size() + inRangeCell2Positions.size());
        const auto cell2Begin = ownerships.begin() + static_cast<std::ptrdiff_t>(cell2IDOffset);
        const std::vector<OwnershipState> ownerships1(ownerships.begin(), cell2Begin);
        const std::vector<OwnershipState> ownerships2(cell2Begin, ownerships.end());

        // Prune dummies: a dummy particle's own list stays empty and a dummy is never listed as a neighbor. The lists
        // of one cell are indexed by the position within that cell but hold the IDs of the other cell's particles.
        auto expectedListsCell1 = cell2OutOfRange ? noNeighborsCell1 : baseExpectedListsCell1;
        for (size_t i = 0; i < expectedListsCell1.size(); ++i) {
          if (ownerships1[i] == OwnershipState::dummy) {
            expectedListsCell1[i].clear();
          } else {
            std::erase_if(expectedListsCell1[i],
                          [&](size_t id) { return ownerships2[id - cell2IDOffset] == OwnershipState::dummy; });
          }
        }

        // Cell 2's lists are only written if there is something in range and the interaction is written back to both
        // particles, i.e. newton3 is used but newton3 lists are not gathered.
        const bool cell2Written = not cell2OutOfRange and mode.n3Used and not mode.gatherN3Lists;
        auto expectedListsCell2 = cell2Written ? baseExpectedListsCell2 : noNeighborsCell2;
        for (size_t i = 0; i < expectedListsCell2.size(); ++i) {
          if (ownerships2[i] == OwnershipState::dummy) {
            expectedListsCell2[i].clear();
          } else {
            std::erase_if(expectedListsCell2[i], [&](size_t id) { return ownerships1[id] == OwnershipState::dummy; });
          }
        }

        params.push_back({.name = mode.name + ownershipVariant.name + (cell2OutOfRange ? "AllOutOfRange" : ""),
                          .positions1 = cell1Positions,
                          .positions2 = cell2OutOfRange ? outOfRangeCell2Positions : inRangeCell2Positions,
                          .ownerships1 = ownerships1,
                          .ownerships2 = ownerships2,
                          .n3Mode = mode,
                          .expectedNeighbors1 = expectedListsCell1,
                          .expectedNeighbors2 = expectedListsCell2});
      }
    }
  }
  return params;
}

}  // namespace

// ---------------------------------------------------------------------------------------------------------------------
// AoSFunctor
// ---------------------------------------------------------------------------------------------------------------------

/**
 * Tests that applying the AoSFunctor over all particle pairs produces the expected neighbor lists.
 */
TEST_P(InteractionListGeneratorFunctorListGenTest, AoSFunctor) {
  const auto &params = GetParam();
  FMCell cell;
  fillCell(cell, params.positions, params.ownerships);
  InteractionListTypeAoS map;
  initMap(map, cell);

  InteractionListGenFunc functor(map, cutoffLength, params.n3Mode.gatherN3Lists);
  applyAoSOverAllPairs(functor, cell, params.n3Mode.n3Used);

  if (params.n3Mode.gatherN3Lists) {
    sortIntoLowerIDLists(map);
  }
  compareWithExpectedLists(map, cell, params.expectedNeighbors);
}

/**
 * Tests the directional semantics of a single AoSFunctor call on one in-range pair, which the parameterized tests above
 * cannot observe:
 * - gatherN3Lists=false, newton3=false: only the first particle's list is filled.
 * - gatherN3Lists=false, newton3=true: both particles' lists are filled by the one call.
 * - gatherN3Lists=true, newton3=true: the pair is stored exactly once.
 */
TEST_F(InteractionListGeneratorFunctorTest, AoSSingleCallDirectionality) {
  const auto run = [](bool gatherN3Lists, bool newton3) {
    FMCell cell;
    fillCell(cell, {{0., 0., 0.}, {1., 0., 0.}}, std::vector<OwnershipState>(2, OwnershipState::owned));
    InteractionListTypeAoS map;
    initMap(map, cell);
    InteractionListGenFunc functor(map, cutoffLength, gatherN3Lists);
    functor.AoSFunctor(cell[0], cell[1], newton3);
    if (gatherN3Lists) {
      // Sort such that higher ID particles are in the lists of lower ID particles with N3 lists, to avoid any
      // arbitrariness from the functor
      sortIntoLowerIDLists(map);
    }
    return std::make_pair(neighborIDs(map, cell, 0), neighborIDs(map, cell, 1));
  };

  EXPECT_EQ(run(/*gatherN3Lists*/ false, /*newton3*/ false),
            std::make_pair(std::vector<size_t>{1}, std::vector<size_t>{}));
  EXPECT_EQ(run(/*gatherN3Lists*/ false, /*newton3*/ true),
            std::make_pair(std::vector<size_t>{1}, std::vector<size_t>{0}));
  EXPECT_EQ(run(/*gatherN3Lists*/ true, /*newton3*/ true),
            std::make_pair(std::vector<size_t>{1}, std::vector<size_t>{}));
}

// ---------------------------------------------------------------------------------------------------------------------
// SoAFunctorSingle
// ---------------------------------------------------------------------------------------------------------------------

/**
 * Tests that SoAFunctorSingle produces the expected neighbor lists for all particles of one cell.
 */
TEST_P(InteractionListGeneratorFunctorListGenTest, SoAFunctorSingle) {
  const auto &params = GetParam();
  FMCell cell;
  fillCell(cell, params.positions, params.ownerships);
  InteractionListTypeAoS map;
  initMap(map, cell);

  InteractionListGenFunc functor(map, cutoffLength, params.n3Mode.gatherN3Lists);
  functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functor.SoAFunctorSingle(cell._particleSoABuffer, params.n3Mode.n3Used);

  if (params.n3Mode.gatherN3Lists) {
    sortIntoLowerIDLists(map);
  }
  compareWithExpectedLists(map, cell, params.expectedNeighbors);
}

// ---------------------------------------------------------------------------------------------------------------------
// SoAFunctorPair
// ---------------------------------------------------------------------------------------------------------------------

/**
 * Tests that SoAFunctorPair produces the expected neighbor lists for all particles of both cells.
 */
TEST_P(InteractionListGeneratorFunctorSoAPairTest, NeighborLists) {
  const auto &params = GetParam();
  FMCell cell1, cell2;
  fillCell(cell1, params.positions1, params.ownerships1, /*idOffset*/ 0);
  fillCell(cell2, params.positions2, params.ownerships2, /*idOffset*/ params.positions1.size());
  InteractionListTypeAoS map;
  initMap(map, cell1);
  initMap(map, cell2);

  InteractionListGenFunc functor(map, cutoffLength, params.n3Mode.gatherN3Lists);
  functor.SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, params.n3Mode.n3Used);

  if (params.n3Mode.gatherN3Lists) {
    sortIntoLowerIDLists(map);
  }
  compareWithExpectedLists(map, cell1, params.expectedNeighbors1);
  compareWithExpectedLists(map, cell2, params.expectedNeighbors2);
}

INSTANTIATE_TEST_SUITE_P(Generated, InteractionListGeneratorFunctorSoAPairTest,
                         ::testing::ValuesIn(getPairListGenerationParams()), ParamNameGenerator());

// ---------------------------------------------------------------------------------------------------------------------
// SoAFunctorVerlet
// ---------------------------------------------------------------------------------------------------------------------

/**
 * Tests that SoAFunctorVerlet produces the expected neighbor lists.
 */
TEST_P(InteractionListGeneratorFunctorListGenTest, SoAFunctorVerlet) {
  const auto &params = GetParam();
  FMCell cell;
  fillCell(cell, params.positions, params.ownerships);
  InteractionListTypeAoS map;
  initMap(map, cell);

  InteractionListGenFunc functor(map, cutoffLength, params.n3Mode.gatherN3Lists);
  functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);
  for (size_t i = 0; i < cell.size(); ++i) {
    functor.SoAFunctorVerlet(cell._particleSoABuffer, i,
                             VerletListTypeSoA(params.verletLists[i].begin(), params.verletLists[i].end()),
                             params.n3Mode.n3Used);
  }

  if (params.n3Mode.gatherN3Lists) {
    sortIntoLowerIDLists(map);
  }
  compareWithExpectedLists(map, cell, params.expectedNeighbors);
}

/**
 * Test cases shared by the AoSFunctor, SoAFunctorSingle and SoAFunctorVerlet tests above.
 */
INSTANTIATE_TEST_SUITE_P(Generated, InteractionListGeneratorFunctorListGenTest,
                         ::testing::ValuesIn(getListGenerationParams()), ParamNameGenerator());

// ---------------------------------------------------------------------------------------------------------------------
// Unsupported option combination and capability queries
// ---------------------------------------------------------------------------------------------------------------------

/**
 * Tests that AoSFunctor, SoAFunctorPair, and SoAFunctorVerlet throw when asked to gather newton3 lists without
 * newton3, which is the unsupported gatherNewton3Lists=true / newton3=false combination. SoAFunctorSingle is
 * deliberately excluded: it ignores its newton3 argument and always behaves as if newton3 were enabled.
 */
TEST_F(InteractionListGeneratorFunctorTest, N3ListsWithoutN3Throws) {
  FMCell cell1, cell2;
  fillCell(cell1, standardPositions, std::vector<OwnershipState>(5, OwnershipState::owned), /*idOffset*/ 0);
  fillCell(cell2, standardPositions, std::vector<OwnershipState>(5, OwnershipState::owned), /*idOffset*/ 5);
  InteractionListTypeAoS map;
  initMap(map, cell1);
  initMap(map, cell2);

  InteractionListGenFunc functor(map, cutoffLength, /*gatherNewton3Lists*/ true);
  functor.SoALoader(cell1, cell1._particleSoABuffer, 0, /*skipSoAResize*/ false);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0, /*skipSoAResize*/ false);

  EXPECT_THROW(functor.AoSFunctor(cell1[0], cell1[1], /*newton3*/ false),
               autopas::utils::ExceptionHandler::AutoPasException);
  EXPECT_THROW(functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, /*newton3*/ false),
               autopas::utils::ExceptionHandler::AutoPasException);
  EXPECT_THROW(functor.SoAFunctorVerlet(cell1._particleSoABuffer, 0, VerletListTypeSoA{1, 2, 3, 4}, /*newton3*/ false),
               autopas::utils::ExceptionHandler::AutoPasException);
}

/**
 * Tests that the functor always allows newton3 but only allows non-newton3 when not gathering newton3 lists,
 * matching the throw behavior of the functor calls.
 */
TEST_F(InteractionListGeneratorFunctorTest, AllowsNonNewton3DependsOnGatherN3Lists) {
  InteractionListTypeAoS map;
  InteractionListGenFunc noN3ListsFunctor(map, cutoffLength, /*gatherNewton3Lists*/ false);
  InteractionListGenFunc n3ListsFunctor(map, cutoffLength, /*gatherNewton3Lists*/ true);

  EXPECT_TRUE(noN3ListsFunctor.allowsNewton3());
  EXPECT_TRUE(noN3ListsFunctor.allowsNonNewton3());
  EXPECT_TRUE(n3ListsFunctor.allowsNewton3());
  EXPECT_FALSE(n3ListsFunctor.allowsNonNewton3());
}

// ---------------------------------------------------------------------------------------------------------------------
// Full flow: list generation + LJFunctor application vs. direct LJFunctor application
// ---------------------------------------------------------------------------------------------------------------------

/**
 * Full-flow test: generates interaction lists with InteractionListGeneratorFunctor on a 3x3x3 grid of particles,
 * applies LJFunctor's AoSFunctor to every list entry. As all interactions returned by InteractionListGeneratorFunctor
 * should be within the cutoff, each application should lead to force contributions (this is checked). Finally, the
 * resulting forces are compared against directly applying LJFunctor's AoSFunctor to all particle pairs.
 */
TEST_P(InteractionListGeneratorFunctorFullFlowTest, ForcesMatchDirectLJApplication) {
  const auto &[mode, ownershipVariant] = GetParam();

  // 3x3x3 grid with spacing 1, with the case's OwnershipStates assigned to the particles in a round-robin.
  constexpr std::array<size_t, 3> particlesPerDim{3, 3, 3};
  const auto ownerships =
      expandOwnerships(ownershipVariant, particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2]);

  FMCell cell, referenceCell;
  for (FMCell *cellPtr : {&cell, &referenceCell}) {
    autopasTools::generators::GridGenerator::fillWithParticles(*cellPtr, particlesPerDim, Molecule{}, {1., 1., 1.},
                                                               {0., 0., 0.});
    for (size_t i = 0; i < cellPtr->size(); ++i) {
      (*cellPtr)[i].setOwnershipState(ownerships[i]);
    }
  }

  // Generate the interaction lists.
  InteractionListTypeAoS map;
  initMap(map, cell);
  InteractionListGenFunc listGenFunctor(map, cutoffLength, mode.gatherN3Lists);
  applyAoSOverAllPairs(listGenFunctor, cell, mode.n3Used);

  LJFunctorType<> ljFunctor(cutoffLength);
  ljFunctor.setParticleProperties(/*epsilon24*/ 24., /*sigmaSquared*/ 1.);

  // Apply the LJ functor to the generated lists.
  for (auto &[particlePtr, neighborList] : map) {
    for (auto *neighborPtr : neighborList) {
      const auto particleForceBefore = particlePtr->getF();
      const auto neighborForceBefore = neighborPtr->getF();
      ljFunctor.AoSFunctor(*particlePtr, *neighborPtr, mode.gatherN3Lists);
      EXPECT_NE(particlePtr->getF(), particleForceBefore)
          << "Force of particle " << particlePtr->getID() << " unchanged by interaction with particle "
          << neighborPtr->getID();
      if (mode.gatherN3Lists) {
        EXPECT_NE(neighborPtr->getF(), neighborForceBefore)
            << "Force of particle " << neighborPtr->getID() << " unchanged by interaction with particle "
            << particlePtr->getID();
      }
    }
  }

  // Create reference forces
  for (size_t i = 0; i < referenceCell.size(); ++i) {
    for (size_t j = i + 1; j < referenceCell.size(); ++j) {
      ljFunctor.AoSFunctor(referenceCell[i], referenceCell[j], /*newton3*/ true);
    }
  }

  // Check forces match
  for (size_t i = 0; i < cell.size(); ++i) {
    for (size_t d = 0; d < 3; ++d) {
      EXPECT_NEAR(cell[i].getF()[d], referenceCell[i].getF()[d], 1e-10)
          << "Force mismatch for particle " << i << " in dimension " << d;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Generated, InteractionListGeneratorFunctorFullFlowTest,
                         ::testing::Combine(::testing::ValuesIn(n3Modes), ::testing::ValuesIn(ownershipVariants)),
                         ParamNameGenerator());

// ---------------------------------------------------------------------------------------------------------------------
// Thread safety
// ---------------------------------------------------------------------------------------------------------------------

/**
 * Tests that AoSFunctor has no data races by generating the particle's lists within a cell with particle-level
 * parallelization. This test will not reliably fail in cases where there are data races, but should reliably fail with
 * a thread sanitizer enabled.
 *
 * Without newton3 a call only writes the first particle's list, so the work can simply be partitioned by that particle.
 * With newton3, we use the four-element colouring scheme above.
 */
TEST_P(IntListGenFuncThreadSafeTest, AoSFunctor) {
  const auto &mode = GetParam();
  constexpr std::array<size_t, 3> particlesPerDim{2, 2, 1};

  // Serial reference.
  FMCell referenceCell;
  autopasTools::generators::GridGenerator::fillWithParticles(referenceCell, particlesPerDim, Molecule{}, {1., 1., 1.},
                                                             {0., 0., 0.});
  InteractionListTypeAoS referenceMap;
  initMap(referenceMap, referenceCell);
  InteractionListGenFunc referenceFunctor(referenceMap, cutoffLength, mode.gatherN3Lists);
  for (size_t i = 0; i < referenceCell.size(); ++i) {
    for (size_t j = 0; j < referenceCell.size(); ++j) {
      if (mode.n3Used ? i < j : i != j) {
        referenceFunctor.AoSFunctor(referenceCell[i], referenceCell[j], mode.n3Used);
      }
    }
  }

  // Parallel application.
  FMCell cell;
  autopasTools::generators::GridGenerator::fillWithParticles(cell, particlesPerDim, Molecule{}, {1., 1., 1.},
                                                             {0., 0., 0.});
  InteractionListTypeAoS map;
  initMap(map, cell);
  InteractionListGenFunc functor(map, cutoffLength, mode.gatherN3Lists);
  {
    const NumThreadGuard numThreadGuard(numThreads);
    if (mode.n3Used) {
      for (const auto &round : coloringRounds) {
        AUTOPAS_OPENMP(parallel for)
        for (auto [i, j] : round) {
          functor.AoSFunctor(cell[i], cell[j], /*newton3*/ true);
        }
      }
    } else {
      AUTOPAS_OPENMP(parallel for)
      for (size_t i = 0; i < cell.size(); ++i) {
        for (size_t j = 0; j < cell.size(); ++j) {
          if (i != j) {
            functor.AoSFunctor(cell[i], cell[j], /*newton3*/ false);
          }
        }
      }
    }
  }

  expectListsMatch(map, cell, referenceMap, referenceCell);
}

/**
 * Tests that SoAFunctorSingle has no data races by having different threads generate the lists of several cells
 * that share one interaction list map. With a thread sanitizer enabled, this should produce errors if the functor is
 * not thread safe.
 */
TEST_P(IntListGenFuncThreadSafeTest, SoAFunctorSingle) {
  const auto &mode = GetParam();

  // Serial reference
  auto referenceCells = makeGridCells();
  InteractionListTypeAoS referenceMap;
  for (auto &referenceCell : referenceCells) {
    initMap(referenceMap, referenceCell);
  }
  InteractionListGenFunc referenceFunctor(referenceMap, cutoffLength, mode.gatherN3Lists);
  for (auto &referenceCell : referenceCells) {
    referenceFunctor.SoALoader(referenceCell, referenceCell._particleSoABuffer, 0, /*skipSoAResize*/ false);
    referenceFunctor.SoAFunctorSingle(referenceCell._particleSoABuffer, mode.n3Used);
  }

  // Parallel application
  auto cells = makeGridCells();
  InteractionListTypeAoS map;
  for (auto &cell : cells) {
    initMap(map, cell);
  }
  InteractionListGenFunc functor(map, cutoffLength, mode.gatherN3Lists);
  for (auto &cell : cells) {
    functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);
  }
  {
    const NumThreadGuard numThreadGuard(numThreads);
    AUTOPAS_OPENMP(parallel for)
    for (size_t cellIndex = 0; cellIndex < cells.size(); ++cellIndex) {
      functor.SoAFunctorSingle(cells[cellIndex]._particleSoABuffer, mode.n3Used);
    }
  }

  for (size_t cellIndex = 0; cellIndex < cells.size(); ++cellIndex) {
    expectListsMatch(map, cells[cellIndex], referenceMap, referenceCells[cellIndex]);
  }
}

/**
 * Tests that SoAFunctorPair has no data races by generating the lists of several cells into one shared map in parallel.
 * With a thread sanitizer enabled, this should produce errors if the functor is not thread safe.
 *
 * Without newton3, we parallelize per outer-loop cell (i.e. the one being written to). With newton3, we use the
 * coloring scheme defined above.
 */
TEST_P(IntListGenFuncThreadSafeTest, SoAFunctorPair) {
  const auto &mode = GetParam();

  // Serial reference
  auto referenceCells = makeGridCells();
  InteractionListTypeAoS referenceMap;
  for (auto &referenceCell : referenceCells) {
    initMap(referenceMap, referenceCell);
  }
  InteractionListGenFunc referenceFunctor(referenceMap, cutoffLength, mode.gatherN3Lists);
  for (auto &referenceCell : referenceCells) {
    referenceFunctor.SoALoader(referenceCell, referenceCell._particleSoABuffer, 0, /*skipSoAResize*/ false);
  }
  for (size_t i = 0; i < referenceCells.size(); ++i) {
    for (size_t j = 0; j < referenceCells.size(); ++j) {
      if (mode.n3Used ? i < j : i != j) {
        referenceFunctor.SoAFunctorPair(referenceCells[i]._particleSoABuffer, referenceCells[j]._particleSoABuffer,
                                        mode.n3Used);
      }
    }
  }

  // Parallel application
  auto cells = makeGridCells();
  InteractionListTypeAoS map;
  for (auto &cell : cells) {
    initMap(map, cell);
  }
  InteractionListGenFunc functor(map, cutoffLength, mode.gatherN3Lists);
  for (auto &cell : cells) {
    functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);
  }
  {
    const NumThreadGuard numThreadGuard(numThreads);
    if (mode.n3Used) {
      for (const auto &round : coloringRounds) {
        AUTOPAS_OPENMP(parallel for)
        for (auto [i, j] : round) {
          functor.SoAFunctorPair(cells[i]._particleSoABuffer, cells[j]._particleSoABuffer, /*newton3*/ true);
        }
      }
    } else {
      AUTOPAS_OPENMP(parallel for)
      for (size_t i = 0; i < cells.size(); ++i) {
        for (size_t j = 0; j < cells.size(); ++j) {
          if (i != j) {
            functor.SoAFunctorPair(cells[i]._particleSoABuffer, cells[j]._particleSoABuffer, /*newton3*/ false);
          }
        }
      }
    }
  }

  for (size_t cellIndex = 0; cellIndex < cells.size(); ++cellIndex) {
    expectListsMatch(map, cells[cellIndex], referenceMap, referenceCells[cellIndex]);
  }
}

/**
 * Tests that SoAFunctorVerlet has no data races when used with cells and a coloring scheme. We mimic pairwise verlet
 * lists (one list per particle per own/neighboring cell), but should cover any case where verlet list is designed
 * around coloring. With a thread sanitizer enabled, this should produce errors if the functor is not thread safe.
 *
 * We test both the internal cell interactions and the cell pair interactions.
 */
TEST_P(IntListGenFuncThreadSafeTest, SoAFunctorVerletPairwiseLists) {
  const auto &mode = GetParam();

  // Four cells, as coloringRounds pairs up four elements.
  constexpr size_t numCells = 4;
  constexpr std::array<size_t, 3> particlesPerDim{2, 2, 2};
  constexpr size_t particlesPerCell = particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2];

  // We create Pairwise Verlet Lists, but assume every particle is in every particle's Verlet list.
  std::vector<std::vector<std::vector<VerletListTypeSoA>>> verletLists(
      numCells,
      std::vector<std::vector<VerletListTypeSoA>>(numCells, std::vector<VerletListTypeSoA>(particlesPerCell)));
  for (size_t c1 = 0; c1 < numCells; ++c1) {
    for (size_t c2 = 0; c2 < numCells; ++c2) {
      for (size_t p = 0; p < particlesPerCell; ++p) {
        const size_t first = c1 * particlesPerCell + p;
        for (size_t q = 0; q < particlesPerCell; ++q) {
          const size_t second = c2 * particlesPerCell + q;
          if (mode.n3Used ? first < second : first != second) {
            verletLists[c1][c2][p].push_back(second);
          }
        }
      }
    }
  }

  // Serial reference. SoAFunctorVerlet addresses all particles through one SoA, so unlike makeGridCells all cells'
  // particles go into a single "cell" (from which one SoA will be loaded) here, in a 2x2 arrangement 3 apart.
  FMCell referenceCell;
  for (size_t cellID = 0; cellID < numCells; ++cellID) {
    const Molecule defaultParticle({0., 0., 0.}, {0., 0., 0.}, cellID * particlesPerCell);
    const size_t cellsX = cellID / 2;
    const size_t cellsY = cellID % 2;
    autopasTools::generators::GridGenerator::fillWithParticles(
        referenceCell, particlesPerDim, defaultParticle, {1., 1., 1.},
        {static_cast<double>(3 * cellsX), static_cast<double>(3 * cellsY), 0.});
  }
  InteractionListTypeAoS referenceMap;
  initMap(referenceMap, referenceCell);
  InteractionListGenFunc referenceFunctor(referenceMap, cutoffLength, mode.gatherN3Lists);
  referenceFunctor.SoALoader(referenceCell, referenceCell._particleSoABuffer, 0, /*skipSoAResize*/ false);
  for (size_t c1 = 0; c1 < numCells; ++c1) {
    for (size_t c2 = 0; c2 < numCells; ++c2) {
      for (size_t p = 0; p < particlesPerCell; ++p) {
        referenceFunctor.SoAFunctorVerlet(referenceCell._particleSoABuffer, c1 * particlesPerCell + p,
                                          verletLists[c1][c2][p], mode.n3Used);
      }
    }
  }

  // Parallel application.
  FMCell cell;
  for (size_t cellID = 0; cellID < numCells; ++cellID) {
    const Molecule defaultParticle({0., 0., 0.}, {0., 0., 0.}, cellID * particlesPerCell);
    const size_t cellsX = cellID / 2;
    const size_t cellsY = cellID % 2;
    autopasTools::generators::GridGenerator::fillWithParticles(
        cell, particlesPerDim, defaultParticle, {1., 1., 1.},
        {static_cast<double>(3 * cellsX), static_cast<double>(3 * cellsY), 0.});
  }
  InteractionListTypeAoS map;
  initMap(map, cell);
  InteractionListGenFunc functor(map, cutoffLength, mode.gatherN3Lists);
  functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);
  {
    const NumThreadGuard numThreadGuard(numThreads);

    AUTOPAS_OPENMP(parallel for)
    for (size_t c = 0; c < numCells; ++c) {
      for (size_t p = 0; p < particlesPerCell; ++p) {
        functor.SoAFunctorVerlet(cell._particleSoABuffer, c * particlesPerCell + p, verletLists[c][c][p], mode.n3Used);
      }
    }

    if (mode.n3Used) {
      for (const auto &round : coloringRounds) {
        AUTOPAS_OPENMP(parallel for)
        for (auto [c1, c2] : round) {
          for (size_t p = 0; p < particlesPerCell; ++p) {
            functor.SoAFunctorVerlet(cell._particleSoABuffer, c1 * particlesPerCell + p, verletLists[c1][c2][p],
                                     /*newton3*/ true);
          }
        }
      }
    } else {
      AUTOPAS_OPENMP(parallel for)
      for (size_t c1 = 0; c1 < numCells; ++c1) {
        for (size_t c2 = 0; c2 < numCells; ++c2) {
          if (c1 != c2) {
            for (size_t p = 0; p < particlesPerCell; ++p) {
              functor.SoAFunctorVerlet(cell._particleSoABuffer, c1 * particlesPerCell + p, verletLists[c1][c2][p],
                                       /*newton3*/ false);
            }
          }
        }
      }
    }
  }

  expectListsMatch(map, cell, referenceMap, referenceCell);
}

INSTANTIATE_TEST_SUITE_P(Generated, IntListGenFuncThreadSafeTest, ::testing::ValuesIn(n3Modes), ParamNameGenerator());

/**
 * Tests that SoAFunctorVerlet has no data races when newton3 is disabled and particles are parallelized at a particle
 * level. With a thread sanitizer enabled, this should produce errors if the functor is not thread safe. Tests gathering
 * n3 and non-n3 lists.
 */
TEST_F(InteractionListGeneratorFunctorTest, ThreadSafeSoAFunctorVerlet) {
  constexpr std::array<size_t, 3> particlesPerDim{2, 2, 1};
  constexpr size_t numParticles = particlesPerDim[0] * particlesPerDim[1] * particlesPerDim[2];

  for (const bool gatherN3Lists : {false, true}) {
    SCOPED_TRACE(gatherN3Lists ? "gatherNewton3Lists = true" : "gatherNewton3Lists = false");

    // Assume every other particle is in the Verlet list of every particle (adjusted for n3).
    std::vector<VerletListTypeSoA> verletLists(numParticles);
    for (size_t i = 0; i < numParticles; ++i) {
      for (size_t j = 0; j < numParticles; ++j) {
        if (gatherN3Lists ? i < j : i != j) {
          verletLists[i].push_back(j);
        }
      }
    }

    // Serial reference.
    FMCell referenceCell;
    autopasTools::generators::GridGenerator::fillWithParticles(referenceCell, particlesPerDim, Molecule{}, {1., 1., 1.},
                                                               {0., 0., 0.});
    InteractionListTypeAoS referenceMap;
    initMap(referenceMap, referenceCell);
    InteractionListGenFunc referenceFunctor(referenceMap, cutoffLength, gatherN3Lists);
    referenceFunctor.SoALoader(referenceCell, referenceCell._particleSoABuffer, 0, /*skipSoAResize*/ false);
    for (size_t i = 0; i < referenceCell.size(); ++i) {
      referenceFunctor.SoAFunctorVerlet(referenceCell._particleSoABuffer, i, verletLists[i], /*newton3*/ gatherN3Lists);
    }

    // Parallel application.
    FMCell cell;
    autopasTools::generators::GridGenerator::fillWithParticles(cell, particlesPerDim, Molecule{}, {1., 1., 1.},
                                                               {0., 0., 0.});
    InteractionListTypeAoS map;
    initMap(map, cell);
    InteractionListGenFunc functor(map, cutoffLength, gatherN3Lists);
    functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize*/ false);
    {
      const NumThreadGuard numThreadGuard(numThreads);
      AUTOPAS_OPENMP(parallel for)
      for (size_t i = 0; i < cell.size(); ++i) {
        functor.SoAFunctorVerlet(cell._particleSoABuffer, i, verletLists[i], /*newton3*/ gatherN3Lists);
      }
    }

    expectListsMatch(map, cell, referenceMap, referenceCell);
  }
}