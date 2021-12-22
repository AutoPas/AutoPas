/**
 * @file KokkosDirectSumTest.cpp
 * @author lgaertner
 * @date 01.12.21
 */

#include "KokkosDirectSumTest.h"

KokkosDirectSumTest::KokkosDirectSumTest() = default;

TEST_F(KokkosDirectSumTest, testUpdateContainer) {
  auto kokkosDirectSum = autopas::KokkosDirectSum<autopas::Particle>({0.0, 0.0, 0.0}, {5.0, 5.0, 5.0}, 0.0, 1.0);

  autopas::Particle p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  autopas::Particle p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  autopas::Particle p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  autopas::Particle p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  autopas::Particle p5({3.5, 3.5, 2.5}, {0, 0, 0}, 4);
  autopas::Particle p6({5.5, 5.5, 5.5}, {0, 0, 0}, 5);

  kokkosDirectSum.addParticleImpl(p1);
  kokkosDirectSum.addParticleImpl(p2);
  kokkosDirectSum.addHaloParticleImpl(p3);
  kokkosDirectSum.addHaloParticleImpl(p4);
  kokkosDirectSum.addParticleImpl(p5);
  kokkosDirectSum.addParticleImpl(p6);
  EXPECT_EQ(kokkosDirectSum.getNumParticles(), 6);
  EXPECT_TRUE(kokkosDirectSum.getIsDirty());

  auto outParticles = kokkosDirectSum.updateContainer(false);

  kokkosDirectSum.resortContainerAndDeleteDummies();

  // after updating the container, halo particles should be removed, owned particles remaining bunched together
  size_t summedIndices = 0;
  kokkosDirectSum.reduce([&](autopas::Particle &p, size_t &r) { r += p.getID(); }, summedIndices,
                         autopas::IteratorBehavior::ownedOrHaloOrDummy);

  EXPECT_EQ(summedIndices, 5);
  EXPECT_EQ(kokkosDirectSum.getNumParticles(), 3);

  EXPECT_EQ(outParticles.size(), 1);

  // check cells[] metadata
  Kokkos::View<autopas::KokkosParticleCell<autopas::Particle> *> cellsMirror = kokkosDirectSum.getCellsHost();
  EXPECT_EQ(cellsMirror[0].begin, 0);
  EXPECT_EQ(cellsMirror[0].cellSize, 3);
  EXPECT_EQ(cellsMirror[1].begin, 3);
  EXPECT_EQ(cellsMirror[1].cellSize, 0);
}

TEST_F(KokkosDirectSumTest, testUpdateHaloParticle) {
  auto kokkosDirectSum = autopas::KokkosDirectSum<autopas::Particle>({0.0, 0.0, 0.0}, {5.0, 5.0, 5.0}, 0.0, 1.0);

  autopas::Particle p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  autopas::Particle p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  autopas::Particle p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  autopas::Particle p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  autopas::Particle p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  autopas::Particle p6({3.5, 3.5, 3.5}, {0, 0, 0}, 5);

  kokkosDirectSum.addParticleImpl(p1);
  kokkosDirectSum.addParticleImpl(p2);
  kokkosDirectSum.addHaloParticleImpl(p3);
  kokkosDirectSum.addHaloParticleImpl(p4);
  kokkosDirectSum.addParticleImpl(p5);
  EXPECT_EQ(kokkosDirectSum.getNumParticles(), 5);
  EXPECT_TRUE(kokkosDirectSum.getIsDirty());

  kokkosDirectSum.updateContainer(false);
  bool updated = kokkosDirectSum.updateHaloParticle(p3);
  EXPECT_TRUE(updated);

  updated = kokkosDirectSum.updateHaloParticle(p6);
  EXPECT_FALSE(updated);

  kokkosDirectSum.resortContainerAndDeleteDummies();

  EXPECT_EQ(kokkosDirectSum.getNumParticles(), 4);
  EXPECT_FALSE(kokkosDirectSum.getIsDirty());

  // check cells[] metadata
  Kokkos::View<autopas::KokkosParticleCell<autopas::Particle> *> cellsMirror = kokkosDirectSum.getCellsHost();
  EXPECT_EQ(cellsMirror[0].begin, 0);
  EXPECT_EQ(cellsMirror[0].cellSize, 3);
  EXPECT_EQ(cellsMirror[1].begin, 3);
  EXPECT_EQ(cellsMirror[1].cellSize, 1);
}
