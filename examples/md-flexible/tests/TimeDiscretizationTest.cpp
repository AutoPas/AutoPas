/**
 * @file TimeDiscretizationTest.cpp
 * @author N. Fottner
 * @date 05/22/19.
 */

#include "TimeDiscretizationTest.h"

#include <memory>

namespace {
void fillWithParticlesAndInit(autopas::AutoPas<ParticleType> &autopasContainer) {
  // Init autopas
  autopasContainer.setBoxMin({0., 0., 0.});
  autopasContainer.setBoxMax({5., 5., 5.});
  autopasContainer.init();

  // Define dummy particle
  ParticleType dummy;
  dummy.setF({0., 0., 1.});
  dummy.setV({0., 0., 1.});
#ifdef MD_FLEXIBLE_USE_MULTI_SITE
  dummy.setTorque({1., 0., 0.});
  dummy.setAngularVel({1., 0., 0.});
  dummy.setQ({0.7071067811865475, 0.7071067811865475, 0., 0.});
#endif

  // Use dummy to fill container
  autopasTools::generators::GridGenerator::fillWithParticles(autopasContainer, {2, 2, 2}, dummy, {1, 1, 1}, {0., 0., 0.});
}

void initPPL(ParticlePropertiesLibrary<> &PPL) {
  PPL.addSiteType(0, 1., 1., 1.);
#ifdef MD_FLEXIBLE_USE_MULTI_SITE
  PPL.addMolType(0, {0, 0}, {{-0.05, 0, 0}, {0.05, 0, 0}}, {1., 1., 1.});
#endif
  PPL.calculateMixingCoefficients();
}
}

template<class MoleculeType> void testCalculateAngularVelocitiesImpl() {
  // todo
}

template<> void testCalculateAngularVelocitiesImpl<MultisiteMolecule>() {

}


TEST_F(TimeDiscretizationTest, testCalculateVelocities) {
  auto autoPas = std::make_shared<autopas::AutoPas<ParticleType>>();
  auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);

  fillWithParticlesAndInit(*autoPas);
  initPPL(*PPL);

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

// todo this really should be extended to include non-zero global force calculations
TEST_F(TimeDiscretizationTest, testCalculatePositions) {
  auto autoPas = std::make_shared<autopas::AutoPas<ParticleType>>();
  auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);

  fillWithParticlesAndInit(*autoPas);
  initPPL(*PPL);

  // The reference positions are the position of the particles in the AutoPas container before
  // calling calculatePositions.
  const std::vector<std::array<double, 3>> referencePositions1 = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
                                                                  {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

  // Create vector of the forces the particles were initialized with
  std::vector<std::array<double, 3>> referenceForces{};
  referenceForces.reserve(8);
  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    referenceForces.push_back(iter->getF());
  }

  size_t index = 0;
  TimeDiscretization::calculatePositionsAndUpdateForces(*autoPas, *PPL, 0.1, {0., 0., 0.}, false);
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

    // expect old force to be the force the particles were initialized with
    EXPECT_DOUBLE_EQ(iter->getOldF()[0], referenceForces[index][0]);
    EXPECT_DOUBLE_EQ(iter->getOldF()[1], referenceForces[index][1]);
    EXPECT_DOUBLE_EQ(iter->getOldF()[2], referenceForces[index][2]);

    // set force and velocity for next iteration
    iter->setF({0, 0, 2});
    iter->setV({0, 0, .5});

    ++index;
  }

  // The reference positions are the position of the particles in the AutoPas container before
  // calling calculatePositions.
  const std::vector<std::array<double, 3>> referencePositions2 = {{0, 0, 0.105}, {1, 0, 0.105}, {0, 1, 0.105}, {1, 1, 0.105},
                                                                  {0, 0, 1.105}, {1, 0, 1.105}, {0, 1, 1.105}, {1, 1, 1.105}};

  TimeDiscretization::calculatePositionsAndUpdateForces(*autoPas, *PPL, 0.1, {0., 0., 0.}, false);
  index = 0;

  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(iter->getR()[0], referencePositions2[index][0]);
    EXPECT_EQ(iter->getR()[1], referencePositions2[index][1]);
    // Störmer-Verlet: 0.1 * .5 + 0.1^2 * (2 / 2) = 0.06
    EXPECT_NEAR(iter->getR()[2], referencePositions2[index][2] + 0.06, 1e-13);
    ++index;
  }

  // Check that force is reset correctly to some non-zero global force.
  TimeDiscretization::calculatePositionsAndUpdateForces(*autoPas, *PPL, 0.1, {-4.5, 2.3, 0.01}, false);
  for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    EXPECT_DOUBLE_EQ(iter->getF()[0], -4.5);
    EXPECT_DOUBLE_EQ(iter->getF()[0], 2.3);
    EXPECT_DOUBLE_EQ(iter->getF()[0], 0.01);
  }
}

TEST_F(TimeDiscretizationTest, testCalculateAngularVelocities) {
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mul;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::quaternion::rotatePosition;
  using autopas::utils::quaternion::rotatePositionBackwards;
  using autopas::utils::quaternion::qMul;
  using autopas::utils::quaternion::qConjugate;
  using autopas::utils::quaternion::convertQuaternionTo3DVec;

  auto autopasContainer = std::make_shared<autopas::AutoPas<ParticleType>>();
  auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);

  const double deltaT = 0.1;

  // Init autopas
  autopasContainer->setBoxMin({0., 0., 0.});
  autopasContainer->setBoxMax({4., 4., 4.});
  autopasContainer->init();

  // Init PPL
  PPL->addSiteType(0, 1, 1, 0.5);
#ifdef MD_FLEXIBLE_USE_MULTI_SITE
  const std::array<double, 3> momentOfInertiaM = {5.23606798, 0.76393202, 6.};
  PPL->addMolType(0, {0, 0, 0},
                  {{0.74349607, 1.20300191, 0.}, {0.3249197, -1.37638192, 0.}, {-1.37638192, -0.3249197, 0.}},
                  momentOfInertiaM);
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

  // add particles
  ParticleType mol1;
  mol1.setR({2., 2., 2.});
  mol1.setQ({0.7071067811865475, 0.7071067811865475, 0., 0.});
  mol1.setF({0., 0., 1.});
  mol1.setV({0., 0., 1.});
  mol1.setTorque({1., 0., 0.});
  mol1.setAngularVel({1., 0., 0.});
  mol1.setTypeId(0);
  mol1.setID(0);
  autopasContainer->addParticle(mol1);

  ParticleType mol2;
  mol2.setR({3., 3., 3.});
  mol2.setQ({-0.7071067811865475, 0.7071067811865475, 0., 0.});
  mol2.setF({0., 0., 1.});
  mol2.setV({0., 0., 1.});
  mol2.setTorque({0., 1., 0.});
  mol2.setAngularVel({1., 0., 0.});
  mol2.setTypeId(0);
  mol2.setID(1);
  autopasContainer->addParticle(mol2);

  // Derive expected Angular Velocities
  // mol 1
  // convert angular velocity to angular momentum (requiring some awkward rotations)
  const auto angVelWHalf1 = mol1.getAngularVel();
  const auto angVelMHalf1 = rotatePositionBackwards(mol1.getQ(), angVelWHalf1);
  const auto angMomMHalf1 = mul(angVelMHalf1, momentOfInertiaM);
  const auto angMomWHalf1 = rotatePosition(mol1.getQ(), angMomMHalf1);

  // half step with angular momentum
  const auto angMomWFullStep1 = add(angMomWHalf1, mulScalar(mol1.getTorque(), 0.5 * deltaT));

  // convert angular momentum back to angular velocity
  const auto angMomMFullStep1 = rotatePositionBackwards(mol1.getQ(), angMomWFullStep1);
  const auto angVelMFullStep1 = div(angMomMFullStep1, momentOfInertiaM);
  const auto angVelWFullStep1 = rotatePosition(mol1.getQ(), angVelMFullStep1);

  // mol 2
  // convert angular velocity to angular momentum (requiring some awkward rotations)
  const auto angVelWHalf2 = mol2.getAngularVel();
  const auto angVelMHalf2 = rotatePositionBackwards(mol2.getQ(), angVelWHalf2);
  const auto angMomMHalf2 = mul(angVelMHalf2, momentOfInertiaM);
  const auto angMomWHalf2 = rotatePosition(mol2.getQ(), angMomMHalf2);

  // half step with angular momentum
  const auto angMomWFullStep2 = add(angMomWHalf2, mulScalar(mol2.getTorque(), 0.5 * deltaT));

  // convert angular momentum back to angular velocity
  const auto angMomMFullStep2 = rotatePositionBackwards(mol2.getQ(), angMomWFullStep2);
  const auto angVelMFullStep2 = div(angMomMFullStep2, momentOfInertiaM);
  const auto angVelWFullStep2 = rotatePosition(mol2.getQ(), angVelMFullStep2);

  // obtain angular velocities as determined by TimeDiscretization::calculateAngularVelocities
  TimeDiscretization::calculateAngularVelocities(*autopasContainer, *PPL, deltaT);

  // compare
  auto mol = autopasContainer->begin(autopas::IteratorBehavior::owned);

  ASSERT_NEAR(mol->getAngularVel()[0], angVelWFullStep1[0], 1e-13);
  ASSERT_NEAR(mol->getAngularVel()[1], angVelWFullStep1[1], 1e-13);
  ASSERT_NEAR(mol->getAngularVel()[2], angVelWFullStep1[2], 1e-13);

  ++mol;

  ASSERT_NEAR(mol->getAngularVel()[0], angVelWFullStep2[0], 1e-13);
  ASSERT_NEAR(mol->getAngularVel()[1], angVelWFullStep2[1], 1e-13);
  ASSERT_NEAR(mol->getAngularVel()[2], angVelWFullStep2[2], 1e-13);

  // Check no additional molecules were created
  ++mol;
  ASSERT_FALSE(mol.isValid());

#else
  PPL->calculateMixingCoefficients();

  // add particles
  ParticleType mol1;
  mol1.setR({2., 2., 2.});
  mol1.setF({0., 0., 1.});
  mol1.setV({0., 0., 1.});
  mol1.setTypeId(0);
  mol1.setID(0);
  autopasContainer->addParticle(mol1);

  // try to calculate angular velocity
  EXPECT_ANY_THROW(TimeDiscretization::calculateAngularVelocities(*autopasContainer, *PPL, deltaT););
#endif
}

TEST_F(TimeDiscretizationTest, testCalculateQuaternion) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayMath::mul;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::cross;
  using autopas::utils::ArrayMath::L2Norm;
  using autopas::utils::ArrayMath::normalize;
  using autopas::utils::quaternion::qMul;
  using autopas::utils::quaternion::qConjugate;
  using autopas::utils::quaternion::convertQuaternionTo3DVec;

  auto autopasContainer = std::make_shared<autopas::AutoPas<ParticleType>>();
  auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);

  const double deltaT = 0.1;

#ifdef MD_FLEXIBLE_USE_MULTI_SITE

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
  ParticleType mol;
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

  // (17)
  const auto angularVelocityM0 = convertQuaternionTo3DVec(qMul(qConjugate(mol.getQ()),qMul(mol.getAngularVel(), mol.getQ())));
  const auto angularMomentumM0 = mul(angularVelocityM0, momentOfInertiaM);

  const auto angularMomentumW0 = convertQuaternionTo3DVec(qMul(mol.getQ(),qMul(angularMomentumM0, qConjugate(mol.getQ())))); // this is used later

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
  auto qHalfK = qHalf0;
  auto qHalfKp1 = qHalf0;
  auto derivQHalfKp1 = derivQHalf0;
  qHalfK[0] = 1e10 * qHalfKp1[0]; // ensuring while loop runs at least once
  while (L2Norm(sub(qHalfKp1, qHalfK)) > 1e-13) {
    qHalfK = qHalfKp1;
    const auto angularMomentumMHalfKp1 = convertQuaternionTo3DVec(qMul(qConjugate(qHalfK), qMul(angularMomentumWHalf, qHalfK)));
    const auto angularVelocityHalfKp1 = div(angularMomentumMHalfKp1,momentOfInertiaM);
    derivQHalfKp1 = mulScalar(qMul(qHalfK, angularVelocityHalfKp1), 0.5);
    qHalfKp1 = normalize(add(mol.getQ(), mulScalar(derivQHalfKp1, 0.5 * deltaT)));
  }

  // (26)
  const auto qExpected = normalize(add(mol.getQ(), mulScalar(derivQHalfKp1, deltaT)));

  // Obtaining angularVelocityWHalf (Not part of original algorithm but needed for implementation in md-flexible)
  const auto angularVelocityMHalf = div(angularMomentumMHalf, momentOfInertiaM);
  const auto angularVelocityWHalfExpected = convertQuaternionTo3DVec(qMul(mol.getQ(),qMul(angularVelocityMHalf, qConjugate(mol.getQ()))));

  // Run TimeDiscretization::calculateQuaternions
  TimeDiscretization::calculateQuaternions(*autopasContainer, *PPL, deltaT, {0, 0, 0});

  auto resultantMol = autopasContainer->begin(autopas::IteratorBehavior::owned);

  // Compare values
  ASSERT_NEAR(qExpected[0], resultantMol->getQ()[0], 1e-13);
  ASSERT_NEAR(qExpected[1], resultantMol->getQ()[1], 1e-13);
  ASSERT_NEAR(qExpected[2], resultantMol->getQ()[2], 1e-13);
  ASSERT_NEAR(qExpected[3], resultantMol->getQ()[3], 1e-13);

  ASSERT_NEAR(angularVelocityWHalfExpected[0], resultantMol->getAngularVel()[0], 1e-13);
  ASSERT_NEAR(angularVelocityWHalfExpected[1], resultantMol->getAngularVel()[1], 1e-13);
  ASSERT_NEAR(angularVelocityWHalfExpected[2], resultantMol->getAngularVel()[2], 1e-13);

  // Confirm no extra molecules were created
  ++resultantMol;
  ASSERT_FALSE(resultantMol.isValid());

#else
  // Init autopas
  autopasContainer->setBoxMin({0., 0., 0.});
  autopasContainer->setBoxMax({4., 4., 4.});
  autopasContainer->init();

  // Init PPL
  PPL->addSiteType(0, 1, 1, 0.5);
  PPL->calculateMixingCoefficients();

  // add particle
  ParticleType mol;
  mol.setR({2., 2., 2.});
  mol.setF({0., 0., 1.});
  mol.setV({0., 0., 1.});
  mol.setTypeId(0);
  autopasContainer->addParticle(mol);

  // Try to calculate the quaternion
  EXPECT_ANY_THROW(TimeDiscretization::calculateQuaternions(*autopasContainer, *PPL, deltaT, {0, 0, 0}));
#endif
}

/**
 * Test the mechanism that throws an exception when particles travel faster than skin/2/rebuildFreq
 */
TEST_F(TimeDiscretizationTest, testFastParticlesCheck) {
  auto autoPas = std::make_shared<autopas::AutoPas<ParticleType>>();
  auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);

  autoPas->setBoxMin({0., 0., 0.});
  autoPas->setBoxMax({10., 10., 10.});
  autoPas->setVerletSkinPerTimestep(.02);
  autoPas->setVerletRebuildFrequency(10);
  autoPas->init();

  initPPL(*PPL);

  // slow particle -> no exception
  autoPas->addParticle(ParticleType({0., 0., 0.}, {0.1, 0., 0.}, 0));
  EXPECT_NO_THROW(
      TimeDiscretization::calculatePositionsAndUpdateForces(*autoPas, PPL, 0.1, {0., 0., 0.}, true));
  // fast particle -> exception
  autoPas->begin()->setV({1., 0., 0.});
  EXPECT_THROW(TimeDiscretization::calculatePositionsAndUpdateForces(*autoPas, PPL, 0.1, {0., 0., 0.}, true),
               std::runtime_error);
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
