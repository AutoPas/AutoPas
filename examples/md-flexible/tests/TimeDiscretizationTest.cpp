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

namespace {
template <class MoleculeType>
void fillWithParticlesAndInit(autopas::AutoPas<MoleculeType> &autopasContainer) {
  // Init autopas
  autopasContainer.setBoxMin({0., 0., 0.});
  autopasContainer.setBoxMax({5., 5., 5.});
  autopasContainer.init();

  // Define dummy particle
  MoleculeType dummy;
  dummy.setF({0., 0., 1.});
  dummy.setV({0., 0., 1.});

  // Use dummy to fill container
  autopasTools::generators::GridGenerator::fillWithParticles(autopasContainer, {2, 2, 2}, dummy, {1, 1, 1}, {0., 0., 0.});
}

template<> void fillWithParticlesAndInit<MultisiteMolecule>(autopas::AutoPas<MultisiteMolecule> &autopasContainer) {
  // Init autopas
  autopasContainer.setBoxMin({0., 0., 0.});
  autopasContainer.setBoxMax({5., 5., 5.});
  autopasContainer.init();

  // Define dummy particle
  MultisiteMolecule dummy;
  dummy.setF({0., 0., 1.});
  dummy.setV({0., 0., 1.});
  dummy.setTorque({1., 0., 0.});
  dummy.setAngularVel({1., 0., 0.});
  dummy.setQ({0.7071067811865475, 0.7071067811865475, 0., 0.});

  // Use dummy to fill container
  autopasTools::generators::GridGenerator::fillWithParticles(autopasContainer, {2, 2, 2}, dummy, {1, 1, 1}, {0., 0., 0.});
}

/**
 * Initialise particle properties library.
 * This function should have a valid molecule type.
 * @tparam MoleculeType
 * @param PPL
 */
template <class MoleculeType>
void initPPL(ParticlePropertiesLibrary<> &PPL) {
  autopas::utils::ExceptionHandler::exception("initPPL should not be called with this molecule type!");
}

template<> void initPPL<Molecule>(ParticlePropertiesLibrary<> &PPL) {
  PPL.addSimpleType(0, 1, 1, 1);
  PPL.calculateMixingCoefficients();
}

template<> void initPPL<MultisiteMolecule>(ParticlePropertiesLibrary<> &PPL) {
  PPL.addSiteType(0, 0.5, 1, 1);
  PPL.addMolType(0, {0, 0}, {{-0.05, 0, 0}, {0.05, 0, 0}}, {1., 1., 1.});
  PPL.calculateMixingCoefficients();
}

template<class MoleculeType> void testCalculateVelocitiesImpl() {
  auto autoPas = std::make_shared<autopas::AutoPas<MoleculeType>>();
  ParticlePropertiesLibrary<double, size_t> PPL();

  fillWithParticlesAndInit(*autoPas);
  initPPL<MoleculeType>(PPL);


}

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
  fillWithParticlesAndInit(*autoPas);

  // The reference positions are the positiosn of the particles in the AutoPas container before
  // calling calculatePositions.
  std::vector<std::array<double, 3>> referencePositions = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
                                                           {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

  size_t index = 0;
  TimeDiscretization::calculatePositions(*autoPas, _particlePropertiesLibrary, 0.1, {0., 0., 0.});
  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    // only change in one direction is expected
    EXPECT_EQ(iter->getR()[0], referencePositions[index][0]);
    EXPECT_EQ(iter->getR()[1], referencePositions[index][1]);
    // Störmer-Verlet: 0.1 * 1 + 0.1^2 * (1 / 2) = 0.105
    EXPECT_NEAR(iter->getR()[2], referencePositions[index][2] + 0.105, 1e-13);

    // expect force to be reset
    EXPECT_THAT(iter->getF(), ::testing::ElementsAreArray({0, 0, 0}));

    // set force and velocity for next iteration
    iter->setF({0, 0, 2});
    iter->setV({0, 0, .5});

    ++index;
  }

  // The reference positions are the positiosn of the particles in the AutoPas container before
  // calling calculatePositions.
  referencePositions = {{0, 0, 0.105}, {1, 0, 0.105}, {0, 1, 0.105}, {1, 1, 0.105},
                        {0, 0, 1.105}, {1, 0, 1.105}, {0, 1, 1.105}, {1, 1, 1.105}};

  TimeDiscretization::calculatePositions(*autoPas, _particlePropertiesLibrary, 0.1, {0., 0., 0.});
  index = 0;

  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(iter->getR()[0], referencePositions[index][0]);
    EXPECT_EQ(iter->getR()[1], referencePositions[index][1]);
    // Störmer-Verlet: 0.1 * .5 + 0.1^2 * (2 / 2) = 0.06
    EXPECT_NEAR(iter->getR()[2], referencePositions[index][2] + 0.06, 1e-13);
    ++index;
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
