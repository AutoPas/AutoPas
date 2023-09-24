/**
 * @file FullParticleCellTest.cpp
 * @author seckler
 * @date 10.09.19
 */

#include "FullParticleCellTest.h"

#include "testingHelpers/commonTypedefs.h"

TEST_F(FullParticleCellTest, testRangeBasedLoop) {
  autopas::FullParticleCell<Particle> cell({1., 1., 1.});

  Particle p({.1, .2, .3}, {0., 0., 0.}, 0);
  cell.addParticle(p);

  Particle p2({.5, .2, .2}, {0., 0., 0.}, 1);
  cell.addParticle(p2);

  EXPECT_EQ(cell.numParticles(), 2);

  for (Particle &particle : cell) {
    particle.setF({42., 42., 42.});
  }

  for (auto p : cell) {
    decltype(p.getF()) comparison = {42., 42., 42};
    ASSERT_EQ(p.getF(), comparison);
  }
}