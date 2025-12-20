/**
* @file PseudoVerletListsTest.cpp
 * @author Lars Doll
 * @date 20.12.2025
 */

#include "PseudoVerletListsTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "testingHelpers/commonTypedefs.h"

TEST_F(PseudoVerletListsTest, OrientationListsHaveCorrectSize) {
  const double cutoff = 1.0;
  const double skin = 0.2;
  const double cellSizeFactor = 1.0;

  autopas::PseudoVerletLists<ParticleFP64> container(
      {0., 0., 0.}, {10., 10., 10.}, cutoff, skin, cellSizeFactor);

  // add at least one particle so cells exist
  ParticleFP64 p({1., 1., 1.}, {0., 0., 0.}, 0);
  container.addParticle(p);

  // rebuild neighbor (orientation) lists
  container.rebuildNeighborLists(nullptr);

  const auto &orientationLists = container.getOrientationLists();

  // one orientation list per cell
  ASSERT_EQ(orientationLists.size(), container.getCells().size());

  // each cell must have exactly 13 SortedCellViews
  for (size_t cellId = 0; cellId < orientationLists.size(); ++cellId) {
    EXPECT_EQ(orientationLists[cellId].size(), 13)
        << "Wrong number of orientation lists in cell " << cellId;
  }
}

TEST_F(PseudoVerletListsTest, OrientationListsAreSorted) {
  const double cutoff = 2.0;
  const double skin = 0.2;
  const double cellSizeFactor = 1.0;

  autopas::PseudoVerletLists<ParticleFP64> container(
      {0., 0., 0.}, {5., 5., 5.}, cutoff, skin, cellSizeFactor);

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

  const auto &orientationLists = container.getOrientationLists();

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
        EXPECT_LE(sortedView._particles[i].first,
                  sortedView._particles[i + 1].first)
            << "SortedCellView not sorted in cell " << cellId
            << ", direction " << dir
            << ", at index " << i;
      }
    }
  }
}

TEST_F(PseudoVerletListsTest, ParticleAppearsOnlyInItsCellSortedViews) {
  const double cutoff = 1.0;
  const double skin = 0.2;
  const double cellSizeFactor = 1.0;

  autopas::PseudoVerletLists<ParticleFP64> container(
      {0., 0., 0.}, {10., 10., 10.}, cutoff, skin, cellSizeFactor);

  // füge ein Partikel hinzu
  ParticleFP64 p({1.5, 1.0, 1.0}, {0., 0., 0.}, 0);
  container.addParticle(p);

  // Orientation Lists rebuilden
  container.rebuildNeighborLists(nullptr);  // Traversal wird hier nicht verwendet

  const auto &orientationLists = container.getOrientationLists();
  const auto &cells = container.getCells();

  // finde die Zelle, in der das Partikel ist
  size_t particleCellId = 0;
  for (size_t i = 0; i < cells.size(); ++i) {
    for (auto &part : cells[i]) {
      if (part.getID() == p.getID()) {
        particleCellId = i;
        break;
      }
    }
  }

  // prüfe nur die SortedCellViews dieser Zelle
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

  EXPECT_TRUE(particleFound) << "Das Partikel sollte in den SortedCellViews seiner eigenen Zelle erscheinen!";

  // optional: prüfe, dass es in keiner anderen Zelle erscheint
  for (size_t i = 0; i < orientationLists.size(); ++i) {
    if (i == particleCellId) continue;
    for (const auto &view : orientationLists[i]) {
      for (const auto &projPair : view._particles) {
        EXPECT_NE(projPair.second->getID(), p.getID())
            << "Partikel taucht fälschlicherweise in einer SortedCellView einer anderen Zelle auf!";
      }
    }
  }
}
