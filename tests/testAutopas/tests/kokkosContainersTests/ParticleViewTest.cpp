/**
 * @file ParticleViewTest.cpp
 *
 * @author lgaertner
 * @date 01.12.2021
 */

#include "ParticleViewTest.h"

ParticleViewTest::ParticleViewTest() = default;

TEST_F(ParticleViewTest, testAddingParticles) {
  auto particleView = ParticleView<autopas::Particle>();

  EXPECT_EQ(particleView.getCapacity(), 8);
  EXPECT_EQ(particleView.getSize(), 0);

  autopas::Particle p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  p1.setOwnershipState(autopas::OwnershipState::owned);
  autopas::Particle p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  p2.setOwnershipState(autopas::OwnershipState::owned);
  autopas::Particle p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  p3.setOwnershipState(autopas::OwnershipState::owned);
  autopas::Particle p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  p4.setOwnershipState(autopas::OwnershipState::halo);
  autopas::Particle p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  p5.setOwnershipState(autopas::OwnershipState::owned);
  autopas::Particle p6({2.5, 1.5, 1.5}, {0, 0, 0}, 5);
  p6.setOwnershipState(autopas::OwnershipState::halo);
  autopas::Particle p7({2.5, 2.5, 2.5}, {0, 0, 0}, 6);
  p7.setOwnershipState(autopas::OwnershipState::halo);
  autopas::Particle p8({2.5, 2.5, 2.5}, {0, 0, 0}, 7);
  p8.setOwnershipState(autopas::OwnershipState::owned);
  autopas::Particle p9({2.5, 2.5, 2.5}, {0, 0, 0}, 8);
  p9.setOwnershipState(autopas::OwnershipState::owned);

  particleView.addParticle(p1);
  particleView.addParticle(p2);
  particleView.addParticle(p3);
  particleView.addParticle(p4);
  particleView.addParticle(p5);

  EXPECT_EQ(particleView.getCapacity(), 8);
  EXPECT_EQ(particleView.getSize(), 5);

  particleView.addParticle(p6);
  particleView.addParticle(p7);
  particleView.addParticle(p8);
  particleView.addParticle(p9);

  EXPECT_EQ(particleView.getCapacity(), 16);
  EXPECT_EQ(particleView.getSize(), 9);
}

TEST_F(ParticleViewTest, testBinningParticles) {
  auto particleView = ParticleView<autopas::Particle>();
  auto particleRange = Kokkos::RangePolicy<>(0ul, 5ul);

  Kokkos::View<autopas::KokkosParticleCell<autopas::Particle> *> cells("testBinningParticles::cells", 2);
  Kokkos::View<autopas::KokkosParticleCell<autopas::Particle> *>::HostMirror cellsHostMirror =
      Kokkos::create_mirror_view(cells);

  autopas::Particle p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  p1.setOwnershipState(autopas::OwnershipState::owned);
  autopas::Particle p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  p2.setOwnershipState(autopas::OwnershipState::halo);
  autopas::Particle p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  p3.setOwnershipState(autopas::OwnershipState::owned);
  autopas::Particle p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  p4.setOwnershipState(autopas::OwnershipState::halo);
  autopas::Particle p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  p5.setOwnershipState(autopas::OwnershipState::owned);

  particleView.addParticle(p1);
  particleView.addParticle(p2);
  particleView.addParticle(p3);
  particleView.addParticle(p4);
  particleView.addParticle(p5);

  // expected ownerships before binning : {o,h,o,h,o}
  auto particlesBeforeBinning = particleView.getParticles();
  int resultBeforeBinning = 1;
  Kokkos::parallel_reduce(
      "testBinningParticles::fillCorrectOwnerships", particleRange,
      KOKKOS_LAMBDA(const size_t &i, int &r) {
        if (i == 0 or i == 2 or i == 4) {
          r = r and particlesBeforeBinning[i].isOwned();
        } else {
          r = r and particlesBeforeBinning[i].isHalo();
        }
      },
      Kokkos::LAnd<int>(resultBeforeBinning));
  EXPECT_TRUE(resultBeforeBinning);

  auto particleBinningLambda = [](autopas::Particle &p) -> size_t { return p.isOwned() ? 0 : 1; };
  particleView.binParticles(particleBinningLambda, cells, "testBinningParticles");

  // expected ownerships after binning : {o,o,o,h,h}
  auto particlesAfterBinning = particleView.getParticles();
  int resultAfterBinning = 1;
  Kokkos::parallel_reduce(
      "testBinningParticles::fillCorrectOwnerships", particleRange,
      KOKKOS_LAMBDA(const size_t &i, int &r) {
        auto p = particlesBeforeBinning[i];
        if (i < 3) {
          r = r and particlesAfterBinning[i].isOwned();
        } else {
          r = r and particlesAfterBinning[i].isHalo();
        }
      },
      Kokkos::LAnd<int>(resultAfterBinning));
  EXPECT_TRUE(resultAfterBinning);

  // also expect unchanged size and capacity
  EXPECT_EQ(particleView.getSize(), 5);
  EXPECT_EQ(particleView.getCapacity(), 8);

  // check indices in 'cells'
  Kokkos::deep_copy(cellsHostMirror, cells);

  EXPECT_EQ(cellsHostMirror[0].begin, 0);
  EXPECT_EQ(cellsHostMirror[0].cellSize, 3);

  EXPECT_EQ(cellsHostMirror[1].begin, 3);
  EXPECT_EQ(cellsHostMirror[1].cellSize, 2);
}

TEST_F(ParticleViewTest, testBinningParticlesAndDeleteDummy) {
  auto particleView = ParticleView<autopas::Particle>();
  auto particleRange = Kokkos::RangePolicy<>(0ul, 5ul);

  Kokkos::View<autopas::KokkosParticleCell<autopas::Particle> *> cells("testBinningParticles::cells", 2);
  Kokkos::View<autopas::KokkosParticleCell<autopas::Particle> *>::HostMirror cellsHostMirror =
      Kokkos::create_mirror_view(cells);

  autopas::Particle p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  p1.setOwnershipState(autopas::OwnershipState::owned);
  autopas::Particle p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  p2.setOwnershipState(autopas::OwnershipState::dummy);
  autopas::Particle p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  p3.setOwnershipState(autopas::OwnershipState::dummy);
  autopas::Particle p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  p4.setOwnershipState(autopas::OwnershipState::halo);
  autopas::Particle p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  p5.setOwnershipState(autopas::OwnershipState::owned);

  particleView.addParticle(p1);
  particleView.addParticle(p2);
  particleView.addParticle(p3);
  particleView.addParticle(p4);
  particleView.addParticle(p5);

  // expected ownerships before binning : {o,d,d,h,o}
  auto particlesBeforeBinning = particleView.getParticles();
  int resultBeforeBinning = 1;
  Kokkos::parallel_reduce(
      "testBinningParticles::fillCorrectOwnerships", particleRange,
      KOKKOS_LAMBDA(const size_t &i, int &r) {
        if (i == 0 or i == 4) {
          r = r and particlesBeforeBinning[i].isOwned();
        } else if (i == 3) {
          r = r and particlesBeforeBinning[i].isHalo();
        } else {
          r = r and particlesBeforeBinning[i].isDummy();
        }
      },
      Kokkos::LAnd<int>(resultBeforeBinning));
  EXPECT_TRUE(resultBeforeBinning);

  auto particleBinningLambda = [](autopas::Particle &p) -> size_t { return p.isOwned() ? 0 : 1; };
  particleView.binParticles(particleBinningLambda, cells, "testBinningParticles");

  // expected ownerships after binning : {o,o,h}
  auto particlesAfterBinning = particleView.getParticles();
  Kokkos::RangePolicy<> newParticleRange(0ul, 3ul);
  int resultAfterBinning = 1;
  Kokkos::parallel_reduce(
      "testBinningParticles::fillCorrectOwnerships", newParticleRange,
      KOKKOS_LAMBDA(const size_t &i, int &r) {
        auto p = particlesBeforeBinning[i];
        if (i < 2) {
          r = r and particlesAfterBinning[i].isOwned();
        } else {
          r = r and particlesAfterBinning[i].isHalo();
        }
      },
      Kokkos::LAnd<int>(resultAfterBinning));
  EXPECT_TRUE(resultAfterBinning);

  // also expect unchanged size and capacity
  EXPECT_EQ(particleView.getSize(), 3);
  EXPECT_EQ(particleView.getCapacity(), 8);

  // check indices in 'cells'
  Kokkos::deep_copy(cellsHostMirror, cells);

  EXPECT_EQ(cellsHostMirror[0].begin, 0);
  EXPECT_EQ(cellsHostMirror[0].cellSize, 2);

  EXPECT_EQ(cellsHostMirror[1].begin, 2);
  EXPECT_EQ(cellsHostMirror[1].cellSize, 1);
}
