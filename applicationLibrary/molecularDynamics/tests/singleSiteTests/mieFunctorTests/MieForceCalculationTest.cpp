#include "MieForceCalculationTest.h"
/**
 * @file ForceCalculationTest.cpp
 * @author F. Gratl
 * @date 13.04.18
 */

#include "autopas/AutoPasDecl.h"
#include "autopasTools/generators/GridGenerator.h"
#include "molecularDynamicsLibrary/MieFunctor.h"
#include "autopas/AutoPasImpl.h"


extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::iteratePairwise(
    mdLib::MieFunctor<Molecule, /* shifting */ false, /*mixing*/ false, autopas::FunctorN3Modes::Both,
                     /*globals*/ false, /*relevantForTuning*/ true> *);

void MieForceCalculationTest::testMie(uint16_t n, uint16_t m, double particleSpacing, double cutoff, autopas::DataLayoutOption dataLayoutOption, std::array<std::array<double, 3>, 4> expectedForces, double tolerance) {
  autopas::AutoPas<Molecule> autoPas;
  std::array<double, 3> boxMin = {0., 0., 0.};
  std::array<double, 3> boxMax = {3., 3., 3.};

  //@todo test this with all containers and traversals
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(cutoff);
  autoPas.setAllowedContainers({autopas::ContainerOption::linkedCells});
  autoPas.setAllowedTraversals({autopas::TraversalOption::lc_c08});
  autoPas.setAllowedDataLayouts({dataLayoutOption});

  autoPas.init();
  Molecule defaultParticle;

  autopasTools::generators::GridGenerator::fillWithParticles(autoPas, {2, 2, 1}, defaultParticle,
                                                             {particleSpacing, particleSpacing, particleSpacing});
  mdLib::MieFunctor<Molecule> functor(cutoff,n,m);
  functor.setParticleProperties(1, 1);

  autoPas.iteratePairwise(&functor);

  for (auto p = autoPas.begin(); p.isValid(); ++p) {
    auto id = p->getID();
    EXPECT_NEAR(expectedForces[id][0], p->getF()[0], tolerance) << "ParticleID: " << id;
    EXPECT_NEAR(expectedForces[id][1], p->getF()[1], tolerance) << "ParticleID: " << id;
    EXPECT_NEAR(expectedForces[id][2], p->getF()[2], tolerance) << "ParticleID: " << id;
  }
}

TEST_F(MieForceCalculationTest, testMiewithU0AoS) {
  double spacing = 1;
  double cutoff = 1.1;

  std::array<std::array<double, 3>, 4> expectedForces = {{{-4, -4, 0}, {4, -4, 0}, {-4, 4, 0}, {4, 4, 0}}};
  double tolerance = 1e-13;

  testMie(2,1, spacing, cutoff, autopas::DataLayoutOption::aos, expectedForces, tolerance);
}

TEST_F(MieForceCalculationTest, testMiewithU0SoA) {
  double spacing = 1;
  double cutoff = 1.1;

  std::array<std::array<double, 3>, 4> expectedForces = {{{-4, -4, 0}, {4, -4, 0}, {-4, 4, 0}, {4, 4, 0}}};
  double tolerance = 1e-13;

  testMie(2,1, spacing, cutoff, autopas::DataLayoutOption::soa, expectedForces, tolerance);
}

TEST_F(MieForceCalculationTest, testMiewithF0AoS) {
  double spacing = 2;
  double cutoff = 2.5;

  std::array<std::array<double, 3>, 4> expectedForces = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
  double tolerance = 1e-13;

  testMie(2,1, spacing, cutoff, autopas::DataLayoutOption::aos, expectedForces, tolerance);
}

TEST_F(MieForceCalculationTest, testMiewithF0SoA) {
  double spacing = 2.0;
  double cutoff = 2.5;

  std::array<std::array<double, 3>, 4> expectedForces = {{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
  double tolerance = 1e-13;

  testMie(2,1, spacing, cutoff, autopas::DataLayoutOption::soa, expectedForces, tolerance);
}

