/**
 * @file PseudoVerletListsTest.cpp
 * @author Lars Doll
 * @date 20.12.2025
 */

#include "PseudoVerletListsTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/pseudoVerletLists/PseudoVerletLists.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/utils/ArrayMath.h"
#include "testingHelpers/commonTypedefs.h"

/**
 *  This test verifies that for each cell, there is an entry with 13 sortedCellViews in OrientationList.
 *  Assumes cellSizeFactor = 1.0
 */
TEST_F(PseudoVerletListsTest, OrientationListsHaveCorrectSize) {
  constexpr double cutoff = 1.0;
  constexpr double skin = 0.2;
  constexpr double cellSizeFactor = 1.0;

  autopas::PseudoVerletLists<ParticleFP64> container({0., 0., 0.}, {10., 10., 10.}, cutoff, skin, cellSizeFactor);

  ParticleFP64 p({1., 1., 1.}, {0., 0., 0.}, 0);
  container.addParticle(p);

  container.rebuildNeighborLists(nullptr);

  const auto &orientationLists = container.getOrientationList();

  // checks if there are as many orientationLists as Cells
  ASSERT_EQ(orientationLists.size(), container.getCells().size());

  // checks if there are 13 directions for each cell
  for (size_t cellId = 0; cellId < orientationLists.size(); ++cellId) {
    EXPECT_EQ(orientationLists[cellId].size(), 13) << "Wrong number of directions in cell " << cellId;
  }
}

/**
 * This test verifies that the particles are stored in the correct order.
 * To do this, it checks in each direction whether the projected distance of the particles is monotonically increasing.
 */
TEST_F(PseudoVerletListsTest, OrientationListsAreSorted) {
  constexpr double cutoff = 2.0;
  constexpr double skin = 0.2;
  constexpr double cellSizeFactor = 1.0;

  autopas::PseudoVerletLists<ParticleFP64> container({0., 0., 0.}, {2.5, 2.5, 2.5}, cutoff, skin, cellSizeFactor);

  container.addParticle(ParticleFP64({1.5, 1.0, 1.0}, {0., 0., 0.}, 0));
  container.addParticle(ParticleFP64({1.0, 1.5, 1.0}, {0., 0., 0.}, 1));
  container.addParticle(ParticleFP64({1.0, 1.0, 1.5}, {0., 0., 0.}, 2));
  container.addParticle(ParticleFP64({0.5, 1.0, 1.0}, {0., 0., 0.}, 3));
  container.addParticle(ParticleFP64({1.0, 0.5, 1.0}, {0., 0., 0.}, 4));
  container.addParticle(ParticleFP64({1.0, 1.0, 0.5}, {0., 0., 0.}, 5));
  container.addParticle(ParticleFP64({0.7, 0.7, 1.0}, {0., 0., 0.}, 6));
  container.addParticle(ParticleFP64({0.7, 1.0, 0.7}, {0., 0., 0.}, 7));
  container.addParticle(ParticleFP64({1.0, 0.7, 0.7}, {0., 0., 0.}, 8));
  container.addParticle(ParticleFP64({0.8, 0.8, 0.8}, {0., 0., 0.}, 9));

  container.rebuildNeighborLists(nullptr);

  const auto &orientationLists = container.getOrientationList();

  // iterate over all cells
  for (size_t cellId = 0; cellId < orientationLists.size(); ++cellId) {
    const auto &cellOrientationLists = orientationLists[cellId];

    // check each SortedCellView
    for (size_t dir = 0; dir < cellOrientationLists.size(); ++dir) {
      const auto &sortedView = cellOrientationLists[dir];

      // views with 0 or 1 particles are trivially sorted
      if (sortedView.size() < 2) continue;

      // check monotonic increase of projected positions
      for (size_t i = 0; i + 1 < sortedView.size(); ++i) {
        EXPECT_LE(sortedView._particles[i].first, sortedView._particles[i + 1].first)
            << "SortedCellView not sorted in cell " << cellId << ", direction " << dir << ", at index " << i;
      }
    }
  }
}

/**
 * This test verifies that a particle exists only in its cell
 * and its corresponding sortedCellViews and cannot be found in any other lists.
 */
TEST_F(PseudoVerletListsTest, ParticleAppearsOnlyInItsCellSortedViews) {
  constexpr double cutoff = 1.0;
  constexpr double skin = 0.2;
  constexpr double cellSizeFactor = 1.0;

  autopas::PseudoVerletLists<ParticleFP64> container({0., 0., 0.}, {10., 10., 10.}, cutoff, skin, cellSizeFactor);

  ParticleFP64 p({1.5, 1.0, 1.0}, {0., 0., 0.}, 0);
  container.addParticle(p);

  container.rebuildNeighborLists(nullptr);

  const auto &orientationLists = container.getOrientationList();
  const auto &cells = container.getCells();

  // find the cell with the particle
  size_t particleCellId = 0;
  for (size_t i = 0; i < cells.size(); ++i) {
    for (auto &part : cells[i]) {
      if (part.getID() == p.getID()) {
        particleCellId = i;
        break;
      }
    }
  }

  // check if the particle is in all sortedCellViews of the cell
  const auto &cellViews = orientationLists[particleCellId];
  bool particleFound = false;
  for (const auto &view : cellViews) {
    for (const auto &projPair : view._particles) {
      if (projPair.second->getID() == p.getID()) {
        particleFound = true;
        break;
      }
    }
    if (particleFound) break;
  }

  EXPECT_TRUE(particleFound) << "The Particle should exist in all SortedCellViews";

  // check that the particle is in no other cell
  for (size_t i = 0; i < orientationLists.size(); ++i) {
    if (i == particleCellId) continue;
    for (const auto &view : orientationLists[i]) {
      for (const auto &projPair : view._particles) {
        EXPECT_NE(projPair.second->getID(), p.getID()) << "Particle exists in SortedCellViews of different cells";
      }
    }
  }
}

/**
 * This test verifies that the directions used in PseudoVerletLists have been set correctly assuming a cellSizeFactor
 * of 1.
 */
TEST_F(PseudoVerletListsTest, DirectionsAreCorrectAndNormalized) {
  constexpr double cutoff = 1.0;
  constexpr double skin = 0.2;
  constexpr double cellSizeFactor = 1.0;

  const autopas::PseudoVerletLists<ParticleFP64> container({0., 0., 0.}, {10., 10., 10.}, cutoff, skin, cellSizeFactor);

  const auto directions = container.getDirections();

  // raw directions as specified
  static constexpr std::array<std::array<double, 3>, 13> rawDirections = {{
      {{1, 0, 0}},
      {{-1, 1, 0}},
      {{0, 1, 0}},
      {{1, 1, 0}},
      {{-1, -1, 1}},
      {{0, -1, 1}},
      {{1, -1, 1}},
      {{-1, 0, 1}},
      {{0, 0, 1}},
      {{1, 0, 1}},
      {{-1, 1, 1}},
      {{0, 1, 1}},
      {{1, 1, 1}},
  }};

  ASSERT_EQ(directions.size(), rawDirections.size());

  for (size_t i = 0; i < rawDirections.size(); ++i) {
    const auto expected = autopas::utils::ArrayMath::normalize(rawDirections[i]);

    const auto &dir = directions[i];

    for (size_t d = 0; d < 3; ++d) {
      constexpr double tol = 1e-12;
      EXPECT_NEAR(dir[d], expected[d], tol) << "Direction mismatch at index " << i << ", component " << d;
    }
  }
}