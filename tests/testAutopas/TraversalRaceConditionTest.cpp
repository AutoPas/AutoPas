#include "TraversalRaceConditionTest.h"

void TraversalRaceConditionTest::fillWithParticles(
    AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> &autoPas,
    std::array<size_t, 3> particlesPerDim) {
  size_t id = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        auto p = PrintableMolecule({x + .5, y + .5, z + .5}, {0, 0, 0}, id++);
        autoPas.addParticle(p);
      }
    }
  }
}

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

//  for (auto &traversalLC : autopas::LinkedCells<PrintableMolecule, FullParticleCell<PrintableMolecule>>::allLCApplicableTraversals()) {
    AutoPas<PrintableMolecule, autopas::FullParticleCell<PrintableMolecule>> autoPas;

    // generates one cell per particle + 1 halo layer
//    autoPas.init(boxMin, boxMax, cellLength, {autopas::ContainerOptions::linkedCells}, {traversalLC});
    autoPas.init(boxMin, boxMax, cellLength, {autopas::ContainerOptions::linkedCells}, {autopas::TraversalOptions::c08});

    fillWithParticles(autoPas, particlesPerDimension);

    SimpleFunctor functor;

    /// @todo: test all traversals -> no autotuning
    autoPas.iteratePairwise(&functor, autopas::aos);

    for (auto particleIterator = autoPas.begin(); particleIterator.isValid(); ++particleIterator) {
      if (particleIterator->getR()[0] == .5 || particleIterator->getR()[0] == particlesPerDimension[0] - .5 ||
          particleIterator->getR()[1] == .5 || particleIterator->getR()[1] == particlesPerDimension[1] - .5 ||
          particleIterator->getR()[2] == .5 || particleIterator->getR()[2] == particlesPerDimension[2] - .5)
        continue;
      // for debugging:
      //    particleIterator->print();

      // although these are doubles this should be exactly zero
      ASSERT_EQ(particleIterator->getF()[0], 0);
      ASSERT_EQ(particleIterator->getF()[1], 0);
      ASSERT_EQ(particleIterator->getF()[2], 0);
  }
//  }
}