/**
 * @file FlopCounterTest.cpp
 * @author F. Gratl
 * @date 01.06.18
 */

#include "FlopCounterTest.h"

/**
 * Generates a square of four particles, iterates over it with the FlopCounter and checks its values
 * @param dataLayoutOption
 */
void FlopCounterTest::test(autopas::DataLayoutOption dataLayoutOption) {
  AutoPas<Particle, FPCell> autoPas;

  autoPas.init({0, 0, 0}, {3, 3, 3}, 1, 0, 1, {autopas::ContainerOptions::directSum});

  std::vector<Particle> molVec{Particle({1, 1, 1}, {0, 0, 0}, 0), Particle({1, 1, 2}, {0, 0, 0}, 1),
                               Particle({1, 2, 1}, {0, 0, 0}, 2), Particle({1, 2, 2}, {0, 0, 0}, 3)};

  for (auto &m : molVec) {
    autoPas.addParticle(m);
  }

  autopas::FlopCounterFunctor<Particle, FPCell> flopCounterFunctor(autoPas.getContainer()->getCutoff());

  autoPas.iteratePairwise(&flopCounterFunctor, dataLayoutOption);

  // every particle checks the distance to all others. Only half of the calculations are made due to Newton 3.
  auto expectedDistanceCalculations = molVec.size() * (molVec.size() - 1) / 2;
  ASSERT_EQ(expectedDistanceCalculations, flopCounterFunctor.getDistanceCalculations());

  // in theory each particle has two in range but only one kernel call because of Newton 3.
  auto expectedKernelCalls = molVec.size();
  ASSERT_EQ(expectedKernelCalls, flopCounterFunctor.getKernelCalls());

  // distance calculations cost 8 flops
  auto expectedFlops = expectedDistanceCalculations * 8 + expectedKernelCalls;
  ASSERT_EQ(expectedFlops, flopCounterFunctor.getFlops(1));

  // two out of three particles are in range
  auto expectedHitRate = 2. / 3.;
  ASSERT_NEAR(expectedHitRate, flopCounterFunctor.getHitRate(), 1e-14);
}

TEST_F(FlopCounterTest, testFlopCounterAoS4Mol) { test(autopas::DataLayoutOption::aos); }

TEST_F(FlopCounterTest, testFlopCounterSoA4Mol) { test(autopas::DataLayoutOption::soa); }
