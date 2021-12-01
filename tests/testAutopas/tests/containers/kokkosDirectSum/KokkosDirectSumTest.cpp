/**
 * @file KokkosDirectSumTest.cpp
 * @author lgaertner
 * @date 01.12.21
 */

#include "KokkosDirectSumTest.h"

KokkosDirectSumTest::KokkosDirectSumTest() = default;

TEST_F(KokkosDirectSumTest, testUpdateContainer) {
  auto kokkosDirectSum = autopas::KokkosDirectSum<autopas::Particle>({0.0,0.0,0.0}, {5.0,5.0,5.0}, 0.0, 1.0);

  autopas::Particle p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  autopas::Particle p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  autopas::Particle p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  autopas::Particle p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  autopas::Particle p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);

  kokkosDirectSum.addParticleImpl(p1);
  kokkosDirectSum.addParticleImpl(p2);
  kokkosDirectSum.addHaloParticleImpl(p3);
  kokkosDirectSum.addHaloParticleImpl(p4);
  kokkosDirectSum.addParticleImpl(p5);
  EXPECT_EQ(kokkosDirectSum.getNumParticles(), 5);

  kokkosDirectSum.updateContainer(false);

  //after updating the container, halo particles should be removed, owned particles remaining bunched together
  size_t summedIndices = 0;
  kokkosDirectSum.reduce([&](autopas::Particle &p, size_t &r){r += p.getID();}, summedIndices, autopas::IteratorBehavior::ownedOrHaloOrDummy);

  EXPECT_EQ(summedIndices, 5);
  EXPECT_EQ(kokkosDirectSum.getNumParticles(), 3);
}