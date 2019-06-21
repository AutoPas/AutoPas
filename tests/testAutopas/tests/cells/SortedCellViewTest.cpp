/**
 * @file SortedCellViewTest.cpp
 * @author C. Menges
 * @date 26.05.2019
 */

#include "SortedCellViewTest.h"
#include "autopas/cells/FullParticleCell.h"

TEST_F(SortedCellViewTest, testParticleAccess) {
  auto fpc = autopas::FullParticleCell<Particle>();
  Particle p1 = Particle();
  fpc.addParticle(p1);
  auto fspc = autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>>(fpc, {1., 0., 0.});
  EXPECT_EQ(fspc._particles.size(), 1);
  std::array<double, 3> force{3.1416, 2.7183, 9.8067};
  fspc._particles.front().second->addF(force);
  EXPECT_THAT(fpc._particles.front().getF(), testing::ContainerEq(force));
}