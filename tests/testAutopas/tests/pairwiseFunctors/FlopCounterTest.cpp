/**
 * @file FlopCounterTest.cpp
 * @author F. Gratl
 * @date 01.06.18
 */

#include "FlopCounterTest.h"

#include "autopas/AutoPas.h"
#include "autopas/pairwiseFunctors/FlopCounterFunctor.h"

/**
 * Generates a square of four particles, iterates over it with the FlopCounter and checks its values
 * @param dataLayoutOption
 */
void FlopCounterTest::test(autopas::DataLayoutOption dataLayoutOption) {
  autopas::AutoPas<Particle> autoPas;

  autoPas.setBoxMin({0, 0, 0});
  autoPas.setBoxMax({3, 3, 3});
  autoPas.setCutoff(1);
  autoPas.setAllowedContainers({autopas::ContainerOption::directSum});
  autoPas.setAllowedTraversals({autopas::TraversalOption::ds_sequential});
  autoPas.setAllowedNewton3Options({autopas::Newton3Option::enabled});
  autoPas.init();

  // made s.t. it is somewhat orthogonal to the projection line for single cells (1,1,1) in SortedCellView
  // -> projection values are all small enough s.t. flopCounterFunctor will always be called
  // SortedCellView is possibly used since DS uses one cell -> CellFunctor
  std::vector<Particle> molVec{Particle({1, 1.14, 1}, {0, 0, 0}, 0), Particle({0.83, 0.54, 1.78}, {0, 0, 0}, 1),
                               Particle({1.63, 0, 1.53}, {0, 0, 0}, 2), Particle({1.8, 0.6, 0.75}, {0, 0, 0}, 3)};
  /*std::vector<Particle> molVec{Particle({1, 1, 1}, {0, 0, 0}, 0), Particle({1, 1, 2}, {0, 0, 0}, 1),
                               Particle({1, 2, 1}, {0, 0, 0}, 2), Particle({1, 2, 2}, {0, 0, 0}, 3)};*/

  for (auto &m : molVec) {
    autoPas.addParticle(m);
  }

  autopas::FlopCounterFunctor<Particle> flopCounterFunctor(autoPas.getCutoff());

  autoPas.iteratePairwise(&flopCounterFunctor);

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

TEST_F(FlopCounterTest, testFlopCounterAoSOpenMP) {
  bool newton3 = true;
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p3({0., 2., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p4({0.1, 2.2, 0.3}, {0., 0., 0.}, 1, 0);

  double cutoff = 1.;

  autopas::FlopCounterFunctor<Particle> functor(cutoff);

  // This is a basic check for the global calculations, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
  {
#if defined(AUTOPAS_OPENMP)
#pragma omp sections
#endif
    {
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      functor.AoSFunctor(p1, p2, newton3);
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      functor.AoSFunctor(p3, p4, newton3);
    }
  }
}

TEST_F(FlopCounterTest, testFlopCounterSoAOpenMP) {
  bool newton3 = true;
  Molecule p1({0., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p2({0.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p3({0., 1., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p4({0.1, 1.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p5({1., 1., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p6({1.1, 1.2, 0.3}, {0., 0., 0.}, 1, 0);

  Molecule p7({1., 0., 0.}, {0., 0., 0.}, 0, 0);
  Molecule p8({1.1, 0.2, 0.3}, {0., 0., 0.}, 1, 0);

  double cutoff = 1.;

  autopas::FlopCounterFunctor<Particle> functor(cutoff);

  autopas::FullParticleCell<Particle> cell1;
  cell1.addParticle(p1);
  cell1.addParticle(p2);

  autopas::FullParticleCell<Particle> cell2;
  cell2.addParticle(p3);
  cell2.addParticle(p4);

  autopas::FullParticleCell<Particle> cell3;
  cell3.addParticle(p5);
  cell3.addParticle(p6);

  autopas::FullParticleCell<Particle> cell4;
  cell3.addParticle(p7);
  cell3.addParticle(p8);

  functor.SoALoader(cell1, cell1._particleSoABuffer, 0);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0);
  functor.SoALoader(cell3, cell3._particleSoABuffer, 0);
  functor.SoALoader(cell4, cell4._particleSoABuffer, 0);

  // This is a basic check for the accumulated values, by checking the handling of two particle interactions in
  // parallel. If interactions are dangerous, archer will complain.

  // first functors on one cell
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
  {
#if defined(AUTOPAS_OPENMP)
#pragma omp sections
#endif
    {
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      functor.SoAFunctorSingle(cell1._particleSoABuffer, newton3);
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      functor.SoAFunctorSingle(cell2._particleSoABuffer, newton3);
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      functor.SoAFunctorSingle(cell3._particleSoABuffer, newton3);
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      functor.SoAFunctorSingle(cell4._particleSoABuffer, newton3);
    }
  }

  // functors on two cells
#if defined(AUTOPAS_OPENMP)
#pragma omp parallel
#endif
  {
#if defined(AUTOPAS_OPENMP)
#pragma omp sections
#endif
    {
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      functor.SoAFunctorPair(cell1._particleSoABuffer, cell2._particleSoABuffer, newton3);
#if defined(AUTOPAS_OPENMP)
#pragma omp section
#endif
      functor.SoAFunctorPair(cell3._particleSoABuffer, cell4._particleSoABuffer, newton3);
    }
  }
}
