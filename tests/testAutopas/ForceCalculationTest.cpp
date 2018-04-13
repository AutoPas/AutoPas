//
// Created by ga68cat on 4/13/18.
//

#include "ForceCalculationTest.h"
#include <algorithm>
#include "../../examples/md/mdutils.h"

void ForceCalculationTest::fillWithParticles(
    AutoPas<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>>
        &autoPas,
    vector<size_t> particlesPerDim, double spacing) {
  size_t id = 0;
  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        auto p = autopas::MoleculeLJ(
            {(x + 1) * spacing, (y + 1) * spacing, (z + 1) * spacing},
            {0, 0, 0}, id++);
        autoPas.addParticle(p);
      }
    }
  }
}

TEST_F(ForceCalculationTest, testLJwithU0AoS) {
  AutoPas<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>>
      autoPas;

  double cutoff = 1.1;
  double epsilon = 1.;
  double sigma = 1.;
  array<double, 3> boxMin = {0., 0., 0.};
  array<double, 3> boxMax = {3., 3., 3.};

  double expectedForces[4][3] = {
      {-24, -24, 0}, {24, -24, 0}, {-24, 24, 0}, {24, 24, 0}};
  double tolerance = 1e-13;

  autoPas.init(boxMin, boxMax, cutoff, linkedCells);

  fillWithParticles(autoPas, {2, 2, 1}, 1);

  autopas::MoleculeLJ::setEpsilon(epsilon);
  autopas::MoleculeLJ::setSigma(sigma);

  LJFunctor<MoleculeLJ, FullParticleCell<MoleculeLJ>>::setGlobals(
      cutoff, epsilon, sigma, 0.0);
  LJFunctor<MoleculeLJ, FullParticleCell<MoleculeLJ>> functor;

  autoPas.iteratePairwise(&functor, autopas::aos);

  for (auto p = autoPas.begin(); p.isValid(); ++p) {
    auto id = p->getID();
    ASSERT_NEAR(expectedForces[id][0], p->getF()[0], tolerance);
    ASSERT_NEAR(expectedForces[id][1], p->getF()[1], tolerance);
    ASSERT_NEAR(expectedForces[id][2], p->getF()[2], tolerance);
  }
}