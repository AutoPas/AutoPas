/**
 * @file ForceCalculationTest.cpp
 * @author F. Gratl
 * @date 13.04.18
 */

#include "ForceCalculationTest.h"

#include "autopas/AutoPas.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopasTools/generators/GridGenerator.h"
#include "testingHelpers/commonTypedefs.h"

void ForceCalculationTest::testLJ(double particleSpacing, double cutoff, autopas::DataLayoutOption dataLayoutOption,
                                  std::array<std::array<double, 3>, 4> expectedForces, double tolerance) {
  autopas::AutoPas<Molecule, FMCell> autoPas;
  std::array<double, 3> boxMin = {0., 0., 0.};
  std::array<double, 3> boxMax = {3., 3., 3.};

  //@todo test this with all containers and traversals
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(cutoff);
  autoPas.setAllowedContainers({autopas::ContainerOption::linkedCells});
  autoPas.setAllowedTraversals({autopas::TraversalOption::c08});
  autoPas.setAllowedDataLayouts({dataLayoutOption});

  autoPas.init();
  Molecule defaultParticle;

  autopasTools::generators::GridGenerator::fillWithParticles(autoPas, {2, 2, 1}, defaultParticle,
                                                             {particleSpacing, particleSpacing, particleSpacing});

  autopas::LJFunctor<Molecule> functor(cutoff);
  functor.setParticleProperties(24, 1);

  autoPas.iteratePairwise(&functor);

  for (auto p = autoPas.begin(); p.isValid(); ++p) {
    auto id = p->getID();
    EXPECT_NEAR(expectedForces[id][0], p->getF()[0], tolerance) << "ParticleID: " << id;
    EXPECT_NEAR(expectedForces[id][1], p->getF()[1], tolerance) << "ParticleID: " << id;
    EXPECT_NEAR(expectedForces[id][2], p->getF()[2], tolerance) << "ParticleID: " << id;
  }
}

TEST_F(ForceCalculationTest, testLJwithU0AoS) {
  double spacing = 1;
  double cutoff = 1.1;

  std::array<std::array<double, 3>, 4> expectedForces = {{{-24, -24, 0}, {24, -24, 0}, {-24, 24, 0}, {24, 24, 0}}};
  double tolerance = 1e-13;

  testLJ(spacing, cutoff, autopas::DataLayoutOption::aos, expectedForces, tolerance);
}

TEST_F(ForceCalculationTest, testLJwithU0SoA) {
  double spacing = 1;
  double cutoff = 1.1;

  std::array<std::array<double, 3>, 4> expectedForces = {{{-24, -24, 0}, {24, -24, 0}, {-24, 24, 0}, {24, 24, 0}}};
  double tolerance = 1e-13;

  testLJ(spacing, cutoff, autopas::DataLayoutOption::soa, expectedForces, tolerance);
}

TEST_F(ForceCalculationTest, testLJwithF0AoS) {
  double spacing = std::pow(2, 1.0 / 6);
  double cutoff = 1.3;

  std::array<std::array<double, 3>, 4> expectedForces = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
  double tolerance = 1e-13;

  testLJ(spacing, cutoff, autopas::DataLayoutOption::aos, expectedForces, tolerance);
}

TEST_F(ForceCalculationTest, testLJwithF0SoA) {
  double spacing = std::pow(2, 1.0 / 6);
  double cutoff = 1.3;

  std::array<std::array<double, 3>, 4> expectedForces = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
  double tolerance = 1e-13;

  testLJ(spacing, cutoff, autopas::DataLayoutOption::soa, expectedForces, tolerance);
}