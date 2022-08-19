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
  auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);

  fillWithParticlesAndInit(*autoPas);
  initPPL<MoleculeType>(*PPL);

  // First timestep
  TimeDiscretization::calculateVelocities(*autoPas, *PPL, 0.1);
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

  // Second timestep
  TimeDiscretization::calculateVelocities(*autoPas, *PPL, 0.1);
  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    // only velocity in one direction is expected
    EXPECT_EQ(iter->getV()[0], 0);
    EXPECT_EQ(iter->getV()[1], 0);
    // Störmer-Verlet: 1.05 + (1+2)/2 * 0.1 = 1.2
    EXPECT_NEAR(iter->getV()[2], 1.2, 1e-13);
  }
}
}

// todo this really should be extended to include non-zero global force calculations
template<class MoleculeType> void testCalculatePositionsImpl() {
  auto autoPas = std::make_shared<autopas::AutoPas<Molecule>>();
  auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);

  fillWithParticlesAndInit(*autoPas);
  initPPL<MoleculeType>(*PPL);

  // The reference positions are the position of the particles in the AutoPas container before
  // calling calculatePositions.
  const std::vector<std::array<double, 3>> referencePositions1 = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
                                                                {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

  size_t index = 0;
  TimeDiscretization::calculatePositions(*autoPas, *PPL, 0.1, {0., 0., 0.});
  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    // only change in one direction is expected
    EXPECT_EQ(iter->getR()[0], referencePositions1[index][0]);
    EXPECT_EQ(iter->getR()[1], referencePositions1[index][1]);
    // Störmer-Verlet: 0.1 * 1 + 0.1^2 * (1 / 2) = 0.105
    EXPECT_NEAR(iter->getR()[2], referencePositions1[index][2] + 0.105, 1e-13);

    // expect force to be reset
    const std::array<double, 3> expectedF = {0., 0., 0.};
    EXPECT_EQ(iter->getF()[0], expectedF[0]);
    EXPECT_EQ(iter->getF()[1], expectedF[1]);
    EXPECT_EQ(iter->getF()[2], expectedF[2]);

    // set force and velocity for next iteration
    iter->setF({0, 0, 2});
    iter->setV({0, 0, .5});

    ++index;
  }

  // The reference positions are the position of the particles in the AutoPas container before
  // calling calculatePositions.
  const std::vector<std::array<double, 3>> referencePositions2 = {{0, 0, 0.105}, {1, 0, 0.105}, {0, 1, 0.105}, {1, 1, 0.105},
                        {0, 0, 1.105}, {1, 0, 1.105}, {0, 1, 1.105}, {1, 1, 1.105}};

  TimeDiscretization::calculatePositions(*autoPas, *PPL, 0.1, {0., 0., 0.});
  index = 0;

  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(iter->getR()[0], referencePositions2[index][0]);
    EXPECT_EQ(iter->getR()[1], referencePositions2[index][1]);
    // Störmer-Verlet: 0.1 * .5 + 0.1^2 * (2 / 2) = 0.06
    EXPECT_NEAR(iter->getR()[2], referencePositions2[index][2] + 0.06, 1e-13);
    ++index;
  }
}

template<class MoleculeType> void testCalculateQuaternionsImpl() {
  // todo
}

// todo extend to include non-zero global force calculations once handling for these have been added
template<> void testCalculateQuaternionsImpl<MultisiteMolecule>() {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayMath::mul;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::cross;
  using autopas::utils::quaternion::qMul;
  using autopas::utils::quaternion::qConjugate;
  using autopas::utils::quaternion::convertQuaternionTo3DVec;

  auto autopasContainer = std::make_shared<autopas::AutoPas<MultisiteMolecule>>();
  auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);

  const double deltaT = 0.1;

  // Init autopas
  autopasContainer->setBoxMin({0., 0., 0.});
  autopasContainer->setBoxMax({4., 4., 4.});
  autopasContainer->init();

  // Init PPL
  PPL->addSiteType(0, 1, 1, 0.5);
  PPL->addMolType(0, {0, 0, 0}, {{0.74349607, 1.20300191, 0.}, {0.3249197, -1.37638192, 0.},
                                 {-1.37638192, -0.3249197, 0.}}, {5.23606798, 0.76393202, 6.});
  // comment on seemingly random site positions + MoI:
  // this molecule looks like
  //
  //             x
  //             |
  //     sqrt(2) |
  //            CoM
  //   sqrt(2)/     \ sqrt(2)
  //        /         \
  //      x             x
  //
  // Site positions have been chosen such the momentOfInertia is diagonal (and thus represented only by 3 doubles)
  PPL->calculateMixingCoefficients();

  // add particle
  MultisiteMolecule mol;
  mol.setR({2., 2., 2.});
  mol.setQ({0.7071067811865475, 0.7071067811865475, 0., 0.});
  mol.setF({0., 0., 1.});
  mol.setV({0., 0., 1.});
  mol.setTorque({1., 0., 0.});
  mol.setAngularVel({1., 0., 0.});
  mol.setTypeId(0);
  autopasContainer->addParticle(mol);

  // derive expected new quaternion by working through algorithm from Rozmanov, 2010, Robust rotational-velocity-Verlet
  // integration methods (method A) step-by-step. To ensure no mistakes, we strictly follow algorithm in paper (as
  // opposed to variation in function)

  const auto momentOfInertiaM = PPL->getMomentOfInertia(mol.getTypeId());

  const auto angularMomentumW0 = mul(mol.getAngularVel(), momentOfInertiaM);

  // (17)
  const auto angularMomentumM0 = convertQuaternionTo3DVec(qMul(qConjugate(mol.getQ()),qMul(angularMomentumW0, mol.getQ())));

  // (18)
  const auto torqueM0 = convertQuaternionTo3DVec(qMul(qConjugate(mol.getQ()),qMul(mol.getTorque(), mol.getQ())));

  // (19)
  const auto angularVelM0 = div(angularMomentumM0, momentOfInertiaM);

  // (20)
  const auto derivAngularMomentumM0 = sub(torqueM0, cross(angularVelM0, angularMomentumM0));

  // (21)
  const auto angularMomentumMHalf = add(angularMomentumM0, mulScalar(derivAngularMomentumM0, 0.5 * deltaT));

  // (22)
  const auto derivQHalf0 = mulScalar(qMul(mol.getQ(), div(angularMomentumMHalf, momentOfInertiaM) ), 0.5);

  // (23)
  const auto qHalf0 = add(mol.getQ(), mulScalar(derivQHalf0, 0.5 * deltaT));

  // (24)
  const auto angularMomentumWHalf = add(angularMomentumW0, mulScalar(mol.getTorque(), 0.5 * deltaT));

  // (25)
  const auto
}


TEST_F(TimeDiscretizationTest, testCalculateVelocitiesSimple) {
    testCalculateVelocitiesImpl<Molecule>();
}

TEST_F(TimeDiscretizationTest, testCalculateVelocitiesMultisite) {
  testCalculateVelocitiesImpl<MultisiteMolecule>();
}

TEST_F(TimeDiscretizationTest, testCalculatePositionsSimple) {
  testCalculatePositionsImpl<Molecule>();
}

TEST_F(TimeDiscretizationTest, testCalculatePositionsMultisite) {
  testCalculatePositionsImpl<MultisiteMolecule>();
}

// @todo: move tests to new class SimulationTest.cpp -> Issue #641
// https://github.com/AutoPas/AutoPas/issues/641
// @note: since this issue was made, these tests have been converted to templates for either Molecule or MultisiteMolecule

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
