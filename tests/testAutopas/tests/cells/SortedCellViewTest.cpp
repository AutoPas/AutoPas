/**
 * @file SortedCellViewTest.cpp
 * @author C. Menges
 * @date 26.05.2019
 */

#include "SortedCellViewTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/SortedCellView.h"
#include "testingHelpers/commonTypedefs.h"

namespace SortedCellViewTest {

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

TEST_F(SortedCellViewTest, testParticleSorting) {
  auto fpc = autopas::FullParticleCell<Particle>();
  Particle p1 = Particle({0., 3., 0.}, {0., 0., 0.}, 0);
  fpc.addParticle(p1);
  Particle p2 = Particle({2., 1., 2.}, {0., 0., 0.}, 2);
  fpc.addParticle(p2);
  Particle p3 = Particle({1., 2., 1.}, {0., 0., 0.}, 1);
  fpc.addParticle(p3);
  Particle p4 = Particle({3., 0., 3.}, {0., 0., 0.}, 3);
  fpc.addParticle(p4);

  unsigned int id = 0u;

  {
    auto fspc = autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>>(fpc, {1., 0., 0.});
    EXPECT_EQ(fspc.numParticles(), 4);

    for (auto &p : fspc._particles) {
      EXPECT_DOUBLE_EQ(p.first, static_cast<double>(id));
      EXPECT_EQ(p.second->getID(), id);
      id++;
    }
  }
  {
    auto fspc = autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>>(fpc, {0., 1., 0.});
    EXPECT_EQ(fspc.numParticles(), 4);

    for (auto &p : fspc._particles) {
      id--;
      EXPECT_DOUBLE_EQ(p.first, static_cast<double>(3u - id));
      EXPECT_EQ(p.second->getID(), id);
    }
  }
}
}  // end namespace SortedCellViewTest
