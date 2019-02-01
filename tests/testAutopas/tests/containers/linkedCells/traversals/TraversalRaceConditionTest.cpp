/**
 * @file TraversalRaceConditionTest.cpp
 * @author seckler
 * @date 18.04.18
 */

#include "TraversalRaceConditionTest.h"
#include "testingHelpers/GridGenerator.h"

/**
 * Idea: create mesh of particles and iterate with the SimpleFunctor.
 * All non-border particles should have F=0 at the end.
 *
 * Failing this test means that the traversal is incomplete or a race condition
 * occurred. Passing this test does not guarantee that there is no race
 * condition. Multiple execution is advised until a deterministic test is
 * implemented.
 *
 * Attention: If the traversal traverses over no particles this test will pass.
 * TODO: when periodic boundaries are implemented also border particles will
 * have F=0
 */
TEST_F(TraversalRaceConditionTest, testRCNonDeterministic) {
  double cellLength = 1;
  std::array<size_t, 3> particlesPerDimension = {30, 30, 30};
  std::array<double, 3> boxMin = {0., 0., 0.};
  std::array<double, 3> boxMax = {(double)particlesPerDimension[0], (double)particlesPerDimension[1],
                                  (double)particlesPerDimension[2]};

#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(8);
#endif

  /// @todo: test all containers
  for (auto &traversalLC : autopas::LinkedCells<Particle, FPCell>::allLCApplicableTraversals()) {
    if (traversalLC == autopas::TraversalOption::c01) {
      // c01 traversal does not work with newton3.
      // Here only one traversal is tested.
      continue;
    }
    autopas::AutoPas<Particle, FPCell> autoPas;

    // generates one cell per particle + 1 halo layer
    auto containerList = {autopas::ContainerOption::linkedCells};
    auto traversalList = {traversalLC};
    autoPas.init();

    auto defaultParticle = Particle({0, 0, 0}, {0, 0, 0}, 0);
    GridGenerator::fillWithParticles(autoPas, particlesPerDimension, defaultParticle);

    SimpleFunctor functor;

    autoPas.iteratePairwise(&functor);

    for (auto particleIterator = autoPas.begin(); particleIterator.isValid(); ++particleIterator) {
      if (particleIterator->getR()[0] == .5 || particleIterator->getR()[0] == particlesPerDimension[0] - .5 ||
          particleIterator->getR()[1] == .5 || particleIterator->getR()[1] == particlesPerDimension[1] - .5 ||
          particleIterator->getR()[2] == .5 || particleIterator->getR()[2] == particlesPerDimension[2] - .5)
        continue;
      // for debugging:
      //    particleIterator->print();

      // although these are doubles this should be exactly zero
      ASSERT_EQ(particleIterator->getF()[0], 0) << "in traversal: " << traversalLC;
      ASSERT_EQ(particleIterator->getF()[1], 0) << "in traversal: " << traversalLC;
      ASSERT_EQ(particleIterator->getF()[2], 0) << "in traversal: " << traversalLC;
    }
  }
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}