/**
 * @file FullParticleCellTest.cpp
 * @author seckler
 * @date 10.09.19
 */

#include "FullParticleCellTest.h"

#include "testingHelpers/commonTypedefs.h"

TEST_F(FullParticleCellTest, testRangeBasedLoop) {
  autopas::FullParticleCell<ParticleFP64> cell({1., 1., 1.});

  ParticleFP64 p({.1, .2, .3}, {0., 0., 0.}, 0);
  cell.addParticle(p);

  ParticleFP64 p2({.5, .2, .2}, {0., 0., 0.}, 1);
  cell.addParticle(p2);

  EXPECT_EQ(cell.size(), 2);

  for (ParticleFP64 &particle : cell) {
    particle.setF({42., 42., 42.});
  }

  for (auto p : cell) {
    decltype(p.getF()) comparison = {42., 42., 42};
    ASSERT_EQ(p.getF(), comparison);
  }
}