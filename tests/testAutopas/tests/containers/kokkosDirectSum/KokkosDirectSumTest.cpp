/**
 * @file LinkedCellsTest.cpp
 * @author lgaertner
 * @date 01.12.21
 */

#include "KokkosDirectSumTest.h"

TYPED_TEST_SUITE_P(KokkosDirectSumTest);

TYPED_TEST_P(KokkosDirectSumTest, testUpdateContainer) {
  decltype(this->_kokkosDirectSum) kokkosDirectSum({0., 0., 0.}, {3., 3., 3.}, 1., 0.);

  autopas::Particle p1({0.5, 0.5, 0.5}, {0, 0, 0}, 0);
  p1.setOwnershipState(autopas::OwnershipState::owned);
  autopas::Particle p2({1.5, 1.5, 1.5}, {0, 0, 0}, 1);
  p2.setOwnershipState(autopas::OwnershipState::halo);
  autopas::Particle p3({1.6, 1.5, 1.5}, {0, 0, 0}, 2);
  p3.setOwnershipState(autopas::OwnershipState::halo);
  autopas::Particle p4({2.5, 1.5, 1.5}, {0, 0, 0}, 3);
  p4.setOwnershipState(autopas::OwnershipState::owned);
  autopas::Particle p5({2.5, 2.5, 2.5}, {0, 0, 0}, 4);
  p5.setOwnershipState(autopas::OwnershipState::owned);

  kokkosDirectSum.addParticle(p1);
  kokkosDirectSum.addParticle(p2);
  kokkosDirectSum.addParticle(p3);
  kokkosDirectSum.addParticle(p4);
  kokkosDirectSum.addParticle(p5);

  kokkosDirectSum.updateContainer(false);
}

REGISTER_TYPED_TEST_SUITE_P(KokkosDirectSumTest, testUpdateContainer);

// INSTANTIATE_TYPED_TEST_SUITE_P(GeneratedTyped, KokkosDirectSumTest);
