//
// Created by johnny on 01.12.23.
//
#include "temporaryTests.h"

#if defined(I_JUST_WANNA_RUN_UNIT_TESTS)
namespace temporaryTests {
void doubleEqual(const double lhs, const double rhs, const double deviation = 0.001){
  if(std::abs(lhs -rhs) > deviation){
    std::cout << "lhs " << lhs << " unequal to rhs " << rhs << std::endl;
    exit(1);
  }
}


void fillWithParticlesAndInit(autopas::AutoPas<ParticleType> &autopasContainer, MoleculeContainer &moleculeContainer,
                              ParticlePropertiesLibrary<> &PPL) {
  autopasContainer.setBoxMin({0., 0., 0.});
  autopasContainer.setBoxMax({5., 5., 5.});
  autopasContainer.init();

  MoleculeType dummyMol;
  dummyMol.setF({0., 0., 1.});
  dummyMol.setTypeId(0);
  dummyMol.setV({0., 0., 1.});

  ParticleType dummySite;
  dummySite.setF({0., 0., 0.});
  dummySite.setTypeId(0);

  // strong similarity with fillWithParticles but different enough using the helper function would not make things significantly easier
  const std::array<int, 3> particlesPerDim{2, 2, 2};
  const std::array<double, 3> spacing{1., 1., 1.};
  int molId = 0;
  int siteId = 0;

  for (unsigned int z = 0; z < particlesPerDim[2]; ++z) {
    for (unsigned int y = 0; y < particlesPerDim[1]; ++y) {
      for (unsigned int x = 0; x < particlesPerDim[0]; ++x) {
        dummyMol.setR({x * spacing[0], y * spacing[1], z * spacing[2]});
        dummyMol.setID(molId++);
        MoleculeType dummyMolCopy =
            dummyMol;  // can be removed if in the future moleculeContainer gets more than the std::move - push_back
        moleculeContainer.push_back(std::move(dummyMolCopy));
        const std::vector<std::array<double, 3>> unrotatedSitePositions = PPL.getSitePositions(dummyMol.getTypeId());
        const std::vector<std::array<double, 3>> rotatedSitePositions =
            autopas::utils::quaternion::rotateVectorOfPositions(dummyMol.getQuaternion(), unrotatedSitePositions);

        for (unsigned int siteIndex = 0; siteIndex < rotatedSitePositions.size(); siteIndex++) {
          const std::array<double, 3> rotatedSitePosition =
              autopas::utils::quaternion::rotatePosition(dummyMol.getQuaternion(), unrotatedSitePositions[siteIndex]);
          dummySite.setR(autopas::utils::ArrayMath::add(dummyMol.getR(), rotatedSitePosition));
          dummySite.setTypeId(PPL.getSiteTypes(dummyMol.getTypeId())[siteIndex]);
          dummySite.setID(siteId++);
          dummySite.setIndexInsideMolecule(siteId);
          dummySite.setMoleculeId(dummyMol.getID());
          autopasContainer.addParticle(dummySite);
        }
      }
    }
  }
}

void initPPL(ParticlePropertiesLibrary<> &PPL) {
#if MD_FLEXIBLE_MODE == MULTISITE
  PPL.addSiteType(0, 1., 1., 0.5);
  PPL.addMolType(0, {0, 0}, {{-0.05, 0, 0}, {0.05, 0, 0}}, {1., 1., 1.});
#else
  PPL.addSiteType(0, 1., 1., 1.);
#endif
  PPL.calculateMixingCoefficients();
}

void testCalculateVelocities() {
  std::cout << "testing calculateVelocities"<< std::endl;

  auto autoPas = std::make_shared<autopas::AutoPas<ParticleType>>();
  auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);
  auto moleculeContainer = MoleculeContainer();

  initPPL(*PPL);
#if MD_FLEXIBLE_MODE != MULTISITE or not defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH)
  fillWithParticlesAndInit(*autoPas);

#else
  fillWithParticlesAndInit(*autoPas, moleculeContainer, *PPL);
#endif

  std::cout << "first time calling calculateVelocities"<< std::endl;
  // First timestep
  TimeDiscretization::calculateVelocities(*autoPas, moleculeContainer, *PPL, 0.1);

  std::cout << "moleculeContainer stuff"<< std::endl;
  for(size_t i = 0; i < moleculeContainer.size(); i++){
    auto& molecule = moleculeContainer.get(i);
    doubleEqual(molecule.getV()[0], 0., 0.);
    doubleEqual(molecule.getV()[1], 0., 0.);
    doubleEqual(molecule.getV()[2], 1.05, 0.0001);
    molecule.setOldF(molecule.getF());
    //molecule.setF({0, 0, 2});
    molecule.setF({0, 0, 2});
  }

  //for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
  //  iter->setF({0,0,1});
  //}

  std::cout << "second time calling calculateVelocities"<< std::endl;
  TimeDiscretization::calculateVelocities(*autoPas, moleculeContainer, *PPL, 0.1);
  for(size_t i = 0; i < moleculeContainer.size(); i++){
    auto& molecule = moleculeContainer.get(i);
    doubleEqual(molecule.getV()[0], 0., 0.);
    doubleEqual(molecule.getV()[1], 0., 0.);
    doubleEqual(molecule.getV()[2], 1.2, 0.0001);
  }
  std::cout << "done testing calculateVelocities"<< std::endl;
}

void testCalculateQuaternion() {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::cross;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::L2Norm;
  using autopas::utils::ArrayMath::mul;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::normalize;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::quaternion::convertQuaternionTo3DVec;
  using autopas::utils::quaternion::qConjugate;
  using autopas::utils::quaternion::qMul;
  using autopas::utils::quaternion::rotateVectorOfPositions;

  auto autopasContainer = std::make_shared<autopas::AutoPas<ParticleType>>();
  auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);
  auto moleculeContainer = MoleculeContainer();

  const double deltaT = 0.1;

#if MD_FLEXIBLE_MODE == MULTISITE

  // Init autopas
  autopasContainer->setBoxMin({0., 0., 0.});
  autopasContainer->setBoxMax({4., 4., 4.});
  autopasContainer->init();

  // Init PPL
  PPL->addSiteType(0, 1, 1, 0.5);
  PPL->addMolType(0, {0, 0, 0},
                  {{0.74349607, 1.20300191, 0.}, {0.3249197, -1.37638192, 0.}, {-1.37638192, -0.3249197, 0.}},
                  {5.23606798, 0.76393202, 6.});
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

  const auto molId = 0;

  MoleculeType mol;
  mol.setR({2., 2., 2.});
  mol.setQuaternion({0., 0., 0., 1.});
  mol.setF({0.,0.,0.}); //mol.setF({0., 0., 1.});
  mol.setV({0., 0., 1.});
  mol.setTorque({0.,0.,0.});
  mol.setAngularVel({1., 0., 0.});
  mol.setTypeId(0);
  mol.setID(molId);

  auto molCopy = mol;
  moleculeContainer.push_back(std::move(molCopy));

  {
    // initialize sites with site forces necessary to apply a torque of {1., 0.,0.} on the molecule
    const auto siteTypes = PPL->getSiteTypes(mol.getTypeId());
    const auto unrotatedSitePositions = PPL->getSitePositions(mol.getTypeId());
    const auto rotatedSitePositions =
        autopas::utils::quaternion::rotateVectorOfPositions(mol.getQuaternion(), unrotatedSitePositions);
    auto id = 0;

    ParticleType site;
    site.setR(autopas::utils::ArrayMath::add(mol.getR(), rotatedSitePositions[0]));
    site.setID(id++);
    site.setTypeId(siteTypes[0]);
    site.setMoleculeId(molId);
    site.setIndexInsideMolecule(0);
    std::array<double, 3> unrotatedFSite0 = {0.,
                                             1/0.74349607,
                                             0.}; // cross(unrotatedSitePosition, F) = {0., 0., 1.};
    std::array<double, 3> rotatedFSite0 = autopas::utils::quaternion::rotatePosition({0.,0.,0.,1}, unrotatedFSite0);
    site.setF(rotatedFSite0);
    autopasContainer->addParticle(site);

    site.setR(autopas::utils::ArrayMath::add(mol.getR(), rotatedSitePositions[1]));
    site.setID(id++);
    site.setMoleculeId(molId);
    site.setTypeId(siteTypes[1]);
    site.setIndexInsideMolecule(1);
    site.setF({0., 0., 0.});
    autopasContainer->addParticle(site);

    site.setR(autopas::utils::ArrayMath::add(mol.getR(), rotatedSitePositions[2]));
    site.setID(id++);
    site.setMoleculeId(molId);
    site.setTypeId(siteTypes[2]);
    site.setIndexInsideMolecule(2);
    site.setF({0., 0., 0.});
    autopasContainer->addParticle(site);
  }

  // derive expected new quaternion by working through algorithm from Rozmanov, 2010, Robust rotational-velocity-Verlet
  // integration methods (method A) step-by-step. To ensure no mistakes, we strictly follow algorithm in paper (as
  // opposed to variation in function)

  const auto momentOfInertiaM = PPL->getMomentOfInertia(mol.getTypeId());

  // (17)
  const auto angularVelocityM0 =
      convertQuaternionTo3DVec(qMul(qConjugate(mol.getQuaternion()), qMul(mol.getAngularVel(), mol.getQuaternion())));
  const auto angularMomentumM0 = mul(angularVelocityM0, momentOfInertiaM);

  const auto angularMomentumW0 = convertQuaternionTo3DVec(
      qMul(mol.getQuaternion(), qMul(angularMomentumM0, qConjugate(mol.getQuaternion()))));  // this is used later

  // (18)
  const std::array<double, 3> expectedTorqueOnMol = {0.,0.,1.};
  const auto torqueM0 =
      convertQuaternionTo3DVec(qMul(qConjugate(mol.getQuaternion()), qMul(expectedTorqueOnMol, mol.getQuaternion())));

  // (19)
  const auto angularVelM0 = div(angularMomentumM0, momentOfInertiaM);

  // (20)
  const auto derivAngularMomentumM0 = sub(torqueM0, cross(angularVelM0, angularMomentumM0));

  // (21)
  const auto angularMomentumMHalf = add(angularMomentumM0, mulScalar(derivAngularMomentumM0, 0.5 * deltaT));

  // (22)
  const auto derivQHalf0 = mulScalar(qMul(mol.getQuaternion(), div(angularMomentumMHalf, momentOfInertiaM)), 0.5);

  // (23)
  const auto qHalf0 = normalize(add(mol.getQuaternion(), mulScalar(derivQHalf0, 0.5 * deltaT)));

  // (24)
  const auto angularMomentumWHalf = add(angularMomentumW0, mulScalar(mol.getTorque(), 0.5 * deltaT));

  // (25)
  auto qHalfK = qHalf0;
  auto qHalfKp1 = qHalf0;
  auto derivQHalfKp1 = derivQHalf0;
  qHalfK[0] += 2e-13;  // ensuring while loop runs at least once
  while (L2Norm(sub(qHalfKp1, qHalfK)) > 1e-13) {
    qHalfK = qHalfKp1;
    const auto angularMomentumMHalfKp1 =
        convertQuaternionTo3DVec(qMul(qConjugate(qHalfK), qMul(angularMomentumWHalf, qHalfK)));
    const auto angularVelocityHalfKp1 = div(angularMomentumMHalfKp1, momentOfInertiaM);
    derivQHalfKp1 = mulScalar(qMul(qHalfK, angularVelocityHalfKp1), 0.5);
    qHalfKp1 = normalize(add(mol.getQuaternion(), mulScalar(derivQHalfKp1, 0.5 * deltaT)));
  }

  // (26)
  const auto qExpected = normalize(add(mol.getQuaternion(), mulScalar(derivQHalfKp1, deltaT)));

  // Obtaining angularVelocityWHalf (Not part of original algorithm but needed for implementation in md-flexible)
  const auto angularVelocityMHalf = div(angularMomentumMHalf, momentOfInertiaM);
  const auto angularVelocityWHalfExpected =
      convertQuaternionTo3DVec(qMul(mol.getQuaternion(), qMul(angularVelocityMHalf, qConjugate(mol.getQuaternion()))));

  // Run TimeDiscretization::calculateQuaternions
  TimeDiscretization::gatherTorquesFromForces(*autopasContainer, moleculeContainer, *PPL);
  TimeDiscretization::calculateQuaternionsAndResetTorques(*autopasContainer, moleculeContainer, *PPL, deltaT, {0, 0, 0});

  //auto resultantMol = autopasContainer->begin(autopas::IteratorBehavior::owned);
  auto resultantMol = moleculeContainer.get(0);

  // Compare values
  doubleEqual(qExpected[0], resultantMol.getQuaternion()[0]);
  doubleEqual(qExpected[1], resultantMol.getQuaternion()[1]);
  doubleEqual(qExpected[2], resultantMol.getQuaternion()[2]);
  doubleEqual(qExpected[3], resultantMol.getQuaternion()[3]);


  doubleEqual(angularVelocityWHalfExpected[0], resultantMol.getAngularVel()[0]);
  doubleEqual(angularVelocityWHalfExpected[1], resultantMol.getAngularVel()[1]);
  doubleEqual(angularVelocityWHalfExpected[2], resultantMol.getAngularVel()[2]);

  doubleEqual(0., resultantMol.getTorque()[0]);
  doubleEqual(0., resultantMol.getTorque()[1]);
  doubleEqual(0., resultantMol.getTorque()[2]);

  // Reset torques using a non-zero global force
  TimeDiscretization::calculateQuaternionsAndResetTorques(*autopasContainer, moleculeContainer, *PPL, deltaT, {0.1, -10, 1.});

  //auto resultantMolWithGlobalForce = autopasContainer->begin(autopas::IteratorBehavior::owned);
  auto resultantMolWithGlobalForce = moleculeContainer.get(0);

  // From the current quaternion, calculate the expected torque
  std::array<double, 3> expectedTorque{0., 0., 0.};
  const auto unrotatedSitePositions = PPL->getSitePositions(0);
  const auto rotatedSitePositions =
      rotateVectorOfPositions(resultantMolWithGlobalForce.getQuaternion(), unrotatedSitePositions);
  for (size_t site = 0; site < PPL->getNumSites(0); site++) {
    expectedTorque = add(expectedTorque, cross(rotatedSitePositions[site], {0.1, -10., 1.}));
  }

  doubleEqual(expectedTorque[0], resultantMolWithGlobalForce.getTorque()[0]);
  doubleEqual(expectedTorque[1], resultantMolWithGlobalForce.getTorque()[1]);
  doubleEqual(expectedTorque[2], resultantMolWithGlobalForce.getTorque()[2]);


  //check that site positions are consistent with current position and rotation of molecule
  for (auto iter = autopasContainer->begin(); iter.isValid(); ++iter) {
    std::array<double, 3> expectedSitePosition{autopas::utils::ArrayMath::add(resultantMolWithGlobalForce.getR(), rotatedSitePositions[iter->getIndexInsideMolecule()])};
    doubleEqual(expectedSitePosition[0], iter->getR()[0]);
    doubleEqual(expectedSitePosition[1], iter->getR()[1]);
    doubleEqual(expectedSitePosition[2], iter->getR()[2]);

  }

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
  EXPECT_ANY_THROW(TimeDiscretization::calculateQuaternionsAndResetTorques(*autopasContainer, *PPL, deltaT, {0, 0, 0}));
#endif
}

void testCalculateAngularVelocities() {

    using autopas::utils::ArrayMath::add;
    using autopas::utils::ArrayMath::div;
    using autopas::utils::ArrayMath::mul;
    using autopas::utils::ArrayMath::mulScalar;
    using autopas::utils::quaternion::convertQuaternionTo3DVec;
    using autopas::utils::quaternion::qConjugate;
    using autopas::utils::quaternion::qMul;
    using autopas::utils::quaternion::rotatePosition;
    using autopas::utils::quaternion::rotatePositionBackwards;

    auto autopasContainer = std::make_shared<autopas::AutoPas<ParticleType>>();
    auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);
    auto moleculeContainer = MoleculeContainer();

    const double deltaT = 0.1;

    // Init autopas
    autopasContainer->setBoxMin({0., 0., 0.});
    autopasContainer->setBoxMax({4., 4., 4.});
    autopasContainer->init();

    // Init PPL
    PPL->addSiteType(0, 1, 1, 0.5);
#if MD_FLEXIBLE_MODE == MULTISITE
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

    // add Molecules;
    // Since calculateAngularVelocities ASSUMES THAT THE TORQUE ACTING ON A MOLECULE HAS ALREADY BEEN COMPUTED
    // based on the forces acting on the molecules the sites don't even have to be initialized
    MoleculeType mol1;
    mol1.setR({2., 2., 2.});
    mol1.setQuaternion({0.7071067811865475, 0.7071067811865475, 0., 0.});
    mol1.setF({0., 0., 1.});
    mol1.setV({0., 0., 1.});
    mol1.setTorque({1., 0., 0.});
    mol1.setAngularVel({1., 0., 0.});
    mol1.setTypeId(0);
    mol1.setID(0);
    {
    auto mol1Copy = mol1;
    moleculeContainer.push_back(std::move(mol1Copy));
    }

    MoleculeType mol2;
    mol2.setR({3., 3., 3.});
    mol2.setQuaternion({-0.7071067811865475, 0.7071067811865475, 0., 0.});
    mol2.setF({0., 0., 1.});
    mol2.setV({0., 0., 1.});
    mol2.setTorque({0., 1., 0.});
    mol2.setAngularVel({1., 0., 0.});
    mol2.setTypeId(0);
    mol2.setID(1);
    {
    auto mol2Copy = mol2;
    moleculeContainer.push_back(std::move(mol2Copy));
    }

    // Derive expected Angular Velocities
    // mol 1
    // convert angular velocity to angular momentum (requiring some awkward rotations)
    const auto angVelWHalf1 = mol1.getAngularVel();
    const auto angVelMHalf1 = rotatePositionBackwards(mol1.getQuaternion(), angVelWHalf1);
    const auto angMomMHalf1 = mul(angVelMHalf1, momentOfInertiaM);
    const auto angMomWHalf1 = rotatePosition(mol1.getQuaternion(), angMomMHalf1);

    // half step with angular momentum
    const auto angMomWFullStep1 = add(angMomWHalf1, mulScalar(mol1.getTorque(), 0.5 * deltaT));

    // convert angular momentum back to angular velocity
    const auto angMomMFullStep1 = rotatePositionBackwards(mol1.getQuaternion(), angMomWFullStep1);
    const auto angVelMFullStep1 = div(angMomMFullStep1, momentOfInertiaM);
    const auto angVelWFullStep1 = rotatePosition(mol1.getQuaternion(), angVelMFullStep1);

    // mol 2
    // convert angular velocity to angular momentum (requiring some awkward rotations)
    const auto angVelWHalf2 = mol2.getAngularVel();
    const auto angVelMHalf2 = rotatePositionBackwards(mol2.getQuaternion(), angVelWHalf2);
    const auto angMomMHalf2 = mul(angVelMHalf2, momentOfInertiaM);
    const auto angMomWHalf2 = rotatePosition(mol2.getQuaternion(), angMomMHalf2);

    // half step with angular momentum
    const auto angMomWFullStep2 = add(angMomWHalf2, mulScalar(mol2.getTorque(), 0.5 * deltaT));

    // convert angular momentum back to angular velocity
    const auto angMomMFullStep2 = rotatePositionBackwards(mol2.getQuaternion(), angMomWFullStep2);
    const auto angVelMFullStep2 = div(angMomMFullStep2, momentOfInertiaM);
    const auto angVelWFullStep2 = rotatePosition(mol2.getQuaternion(), angVelMFullStep2);

    // obtain angular velocities as determined by TimeDiscretization::calculateAngularVelocities
    TimeDiscretization::calculateAngularVelocities(*autopasContainer, moleculeContainer, *PPL, deltaT);

    // compare
    auto computedMol1 = moleculeContainer.get(0);

    doubleEqual(computedMol1.getAngularVel()[0], angVelWFullStep1[0], 1e-15);
    doubleEqual(computedMol1.getAngularVel()[1], angVelWFullStep1[1], 1e-15);
    doubleEqual(computedMol1.getAngularVel()[2], angVelWFullStep1[2], 1e-15);

    auto computedMol2 = moleculeContainer.get(1);

    doubleEqual(computedMol2.getAngularVel()[0], angVelWFullStep2[0], 1e-15);
    doubleEqual(computedMol2.getAngularVel()[1], angVelWFullStep2[1], 1e-15);
    doubleEqual(computedMol2.getAngularVel()[2], angVelWFullStep2[2], 1e-15);

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

/*
void testCalculatePositions() {
    auto autoPas = std::make_shared<autopas::AutoPas<ParticleType>>();
    auto PPL = std::make_shared<ParticlePropertiesLibrary<>>(1.0);
    auto moleculeContainer = MoleculeContainer();

    initPPL(*PPL);
    fillWithParticlesAndInit(*autoPas, moleculeContainer, *PPL);

    // Set verlet skin per timestep to something large so no error messages are displayed
    autoPas->setVerletSkinPerTimestep(1.);

    // The reference positions are the position of the particles in the AutoPas container before
    // calling calculatePositions.
    const std::vector<std::array<double, 3>> referencePositions1 = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
                                                                    {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};

    // Create vector of the forces the particles were initialized with
    std::vector<std::array<double, 3>> referenceForces{};
    referenceForces.reserve(8);
    for(size_t i = 0; i < moleculeContainer.size(); i++) {
      referenceForces.push_back(moleculeContainer.get(i).getF());
    }
    //for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    //referenceForces.push_back(iter->getF());
    //}

    size_t index = 0;
    TimeDiscretization::calculatePositionsAndResetForces(*autoPas, moleculeContainer, *PPL, 0.1, {0., 0., 0.}, false);
    for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    // only change in one direction is expected
    EXPECT_EQ(iter->getR()[0], referencePositions1[index][0]);
    EXPECT_EQ(iter->getR()[1], referencePositions1[index][1]);
    // Störmer-Verlet: 0.1 * 1 + 0.1^2 * (1 / 2) = 0.105
    EXPECT_DOUBLE_EQ(iter->getR()[2], referencePositions1[index][2] + 0.105);

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
    const std::vector<std::array<double, 3>> referencePositions2 = {{0, 0, 0.105}, {1, 0, 0.105}, {0, 1, 0.105},
                                                                    {1, 1, 0.105}, {0, 0, 1.105}, {1, 0, 1.105},
                                                                    {0, 1, 1.105}, {1, 1, 1.105}};

    TimeDiscretization::calculatePositionsAndResetForces(*autoPas, *PPL, 0.1, {0., 0., 0.}, false);
    index = 0;

    for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    EXPECT_EQ(iter->getR()[0], referencePositions2[index][0]);
    EXPECT_EQ(iter->getR()[1], referencePositions2[index][1]);
    // Störmer-Verlet: 0.1 * .5 + 0.1^2 * (2 / 2) = 0.06
    EXPECT_DOUBLE_EQ(iter->getR()[2], referencePositions2[index][2] + 0.06);
    ++index;
    }

    // Check that force is reset correctly to some non-zero global force.
    TimeDiscretization::calculatePositionsAndResetForces(*autoPas, *PPL, 0.1, {-4.5, 2.3, 0.01}, false);
    for (auto iter = autoPas->begin(); iter.isValid(); ++iter) {
    EXPECT_DOUBLE_EQ(iter->getF()[0], -4.5);
    EXPECT_DOUBLE_EQ(iter->getF()[1], 2.3);
    EXPECT_DOUBLE_EQ(iter->getF()[2], 0.01);
    }
}
*/

} //namespace

#endif