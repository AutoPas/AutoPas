/**
 * @file ForceCalculationTest.cpp
 * @author F. Gratl
 * @date 13.04.18
 */

#include "ForceCalculationTest.h"

void ForceCalculationTest::testLJ(double particleSpacing, double cutoff, autopas::DataLayoutOption dataLayoutOption,
                                  std::array<std::array<double, 3>, 4> expectedForces, double tolerance) {
  autopas::AutoPas<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>> autoPas;

  double epsilon = 1.;
  double sigma = 1.;
  std::array<double, 3> boxMin = {0., 0., 0.};
  std::array<double, 3> boxMax = {3., 3., 3.};

  //@todo test this with all containers and traversals
  autoPas.init(boxMin, boxMax, cutoff, 0, 1, {autopas::ContainerOptions::linkedCells},
               {autopas::TraversalOptions::c08});
  autopas::MoleculeLJ defaultParticle;

  GridGenerator::fillWithParticles(autoPas, {2, 2, 1}, defaultParticle,
                                   {particleSpacing, particleSpacing, particleSpacing});

  autopas::MoleculeLJ::setEpsilon(epsilon);
  autopas::MoleculeLJ::setSigma(sigma);

  autopas::LJFunctor<autopas::MoleculeLJ, autopas::FullParticleCell<autopas::MoleculeLJ>> functor(cutoff, epsilon,
                                                                                                  sigma, 0.0);

  autoPas.iteratePairwise(&functor, dataLayoutOption);

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

  testLJ(spacing, cutoff, autopas::aos, expectedForces, tolerance);
}

TEST_F(ForceCalculationTest, testLJwithU0SoA) {
  double spacing = 1;
  double cutoff = 1.1;

  std::array<std::array<double, 3>, 4> expectedForces = {{{-24, -24, 0}, {24, -24, 0}, {-24, 24, 0}, {24, 24, 0}}};
  double tolerance = 1e-13;

  testLJ(spacing, cutoff, autopas::soa, expectedForces, tolerance);
}

TEST_F(ForceCalculationTest, testLJwithF0AoS) {
  double spacing = std::pow(2, 1.0 / 6);
  double cutoff = 1.3;

  std::array<std::array<double, 3>, 4> expectedForces = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
  double tolerance = 1e-13;

  testLJ(spacing, cutoff, autopas::aos, expectedForces, tolerance);
}

TEST_F(ForceCalculationTest, testLJwithF0SoA) {
  double spacing = std::pow(2, 1.0 / 6);
  double cutoff = 1.3;

  std::array<std::array<double, 3>, 4> expectedForces = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
  double tolerance = 1e-13;

  testLJ(spacing, cutoff, autopas::soa, expectedForces, tolerance);
}