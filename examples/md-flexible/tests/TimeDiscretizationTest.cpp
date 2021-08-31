/**
 * @file TimeDiscretizationTest.cpp
 * @author N. Fottner
 * @date 05/22/19.
 */

#include "TimeDiscretizationTest.h"

#include <memory>

#include "autopas/utils/ArrayMath.h"
#include "autopasTools/generators/GridGenerator.h"
#include "src/TimeDiscretization.h"
#include "src/configuration/MDFlexConfig.h"

void TimeDiscretizationTest::fillWithParticlesAndInit(autopas::AutoPas<Molecule> &autopas) {
  autopas.setBoxMin({0., 0., 0.});
  autopas.setBoxMax({5., 5., 5.});
  autopas.init();
  Molecule dummy;
  dummy.setF({0., 0., 1.});
  dummy.setV({0., 0., 1.});
  autopasTools::generators::GridGenerator::fillWithParticles(autopas, {2, 2, 2}, dummy, {1, 1, 1}, {0., 0., 0.});
}

TEST_F(TimeDiscretizationTest, testCalculateVelocities) {
  auto autoPas = std::make_shared<autopas::AutoPas<Molecule>>();
  fillWithParticlesAndInit(*autoPas);

  TimeDiscretization::calculateVelocities(*autoPas, _particlePropertiesLibrary, 0.1);
  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    // only velocity in one direction is expected, as the force is initialized to point only in z-direction.
    EXPECT_EQ(iter->getV()[0], 0);
    EXPECT_EQ(iter->getV()[1], 0);
    // Störmer-Verlet: 1 + (0+1)/2 * 0.1 = 1.05
    EXPECT_NEAR(iter->getV()[2], 1.05, 1e-13);

    // set force for next iteration
    iter->setOldF(iter->getF());
    iter->setF({0, 0, 2});
  }

  TimeDiscretization::calculateVelocities(*autoPas, _particlePropertiesLibrary, 0.1);
  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    // only velocity in one direction is expected
    EXPECT_EQ(iter->getV()[0], 0);
    EXPECT_EQ(iter->getV()[1], 0);
    // Störmer-Verlet: 1.05 + (1+2)/2 * 0.1 = 1.2
    EXPECT_NEAR(iter->getV()[2], 1.2, 1e-13);
  }
}

TEST_F(TimeDiscretizationTest, testCalculatePositions) {
  auto autoPas = std::make_shared<autopas::AutoPas<Molecule>>();
  auto autoPasRef = std::make_shared<autopas::AutoPas<Molecule>>();
  fillWithParticlesAndInit(*autoPas);
  fillWithParticlesAndInit(*autoPasRef);

  TimeDiscretization::calculatePositions(*autoPas, _particlePropertiesLibrary, 0.1);
  for (auto iter = autoPas->begin(), iterRef = autoPasRef->begin(); iter.isValid(); ++iter, ++iterRef) {
    // only change in one direction is expected
    EXPECT_EQ(iter->getR()[0], iterRef->getR()[0]);
    EXPECT_EQ(iter->getR()[1], iterRef->getR()[1]);
    // Störmer-Verlet: 0.1 * 1 + 0.1^2 * (1 / 2) = 0.105
    EXPECT_NEAR(iter->getR()[2], iterRef->getR()[2] + 0.105, 1e-13);

    // expect force to be reset
    EXPECT_THAT(iter->getF(), ::testing::ElementsAreArray({0, 0, 0}));

    // set force and velocity for next iteration
    iter->setF({0, 0, 2});
    iter->setV({0, 0, .5});

    // update reference position
    iterRef->setR(iter->getR());
  }

  TimeDiscretization::calculatePositions(*autoPas, _particlePropertiesLibrary, 0.1);
  for (auto iter = autoPas->begin(), iterRef = autoPasRef->begin(); iter.isValid(); ++iter, ++iterRef) {
    // only velocity in one direction is expected
    EXPECT_EQ(iter->getR()[0], iterRef->getR()[0]);
    EXPECT_EQ(iter->getR()[1], iterRef->getR()[1]);
    // Störmer-Verlet: 0.1 * .5 + 0.1^2 * (2 / 2) = 0.06
    EXPECT_NEAR(iter->getR()[2], iterRef->getR()[2] + 0.06, 1e-13);
  }
}

// @todo: move tests to new class SimulationTest.cpp -> Issue #641
// https://github.com/AutoPas/AutoPas/issues/641

// TEST_F(TimeDiscretizationTest, calculatePairwiseForces) {
//   auto autoPas = std::make_shared<autopas::AutoPas<Molecule>>();
//   fillWithParticlesAndInit(*autoPas);

//   ParticlePropertiesLibraryType particlePropertiesLibrary = ParticlePropertiesLibraryType(3.0);
//   particlePropertiesLibrary.addType(0, 4.0, 4.0, 4.0);
//   particlePropertiesLibrary.calculateMixingCoefficients();

//   bool wasTuningIteration = false;

//   TimeDiscretization::calculatePairwiseForces(*autoPas, particlePropertiesLibrary, 0.1,
//                                               MDFlexConfig::FunctorOption::lj12_6, wasTuningIteration);

//   const std::vector<std::array<double, 3>> expectedForces = {
//       {-3.22083e+09, -3.22083e+09, -3.22083e+09}, {3.22083e+09, -3.22083e+09, -3.22083e+09},
//       {-3.22083e+09, 3.22083e+09, -3.22083e+09},  {3.22083e+09, 3.22083e+09, -3.22083e+09},
//       {-3.22083e+09, -3.22083e+09, 3.22083e+09},  {3.22083e+09, -3.22083e+09, 3.22083e+09},
//       {-3.22083e+09, 3.22083e+09, 3.22083e+09},   {3.22083e+09, 3.22083e+09, 3.22083e+09}};

//   int particleIndex = 0;
//   for (auto particle = autoPas->begin(); particle.isValid(); ++particle) {
//     const std::array<double, 3> force = particle->getF();

//     EXPECT_NEAR(force[0], expectedForces[particleIndex][0], 1e+9);
//     EXPECT_NEAR(force[1], expectedForces[particleIndex][1], 1e+9);
//     EXPECT_NEAR(force[2], expectedForces[particleIndex][2], 1e+9);

//     ++particleIndex;
//   }
// }

// TEST_F(TimeDiscretizationTest, calculateGlobalForces) {
//   auto autoPas = std::make_shared<autopas::AutoPas<Molecule>>();
//   fillWithParticlesAndInit(*autoPas);

//   const std::array<double, 3> globalForce = {0.0, -1.0, 0.0};

//   TimeDiscretization::calculateGlobalForces(*autoPas, globalForce);

//   const std::array<double, 3> expectedForce = {0, -1, 1};

//   for (auto particle = autoPas->begin(); particle.isValid(); ++particle) {
//     const std::array<double, 3> force = particle->getF();

//     EXPECT_EQ(force[0], expectedForce[0]);
//     EXPECT_EQ(force[1], expectedForce[1]);
//     EXPECT_EQ(force[2], expectedForce[2]);
//   }
// }
