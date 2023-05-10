/**
 * @file ReflectiveBoundaryConditionTest.cpp
 * @author F. Gratl
 * @date 26.04.22
 */

#include "ReflectiveBoundaryConditionTest.h"

#include <gmock/gmock-matchers.h>

#include "autopas/AutoPasDecl.h"
#include "autopas/utils/ArrayMath.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"
#include "src/options/BoundaryTypeOption.h"

extern template class autopas::AutoPas<ParticleType>;

/**
 * Very simple test of reflective boundaries in all 3 dimension. Places identical particles on every face and tests that
 * the particle receives the correct force. Get's input parameters from test suite below.
 *
 * Note, this test is not designed to deal with multiple particle types - see reflectiveZoningTest for this.
 *
 * The multi-site force calculation works a bit differently to in RegularGridDecomposition. Here, we mirror the molecule
 * by
 * a) reflecting the unrotated rigid bodies in the dimension given by the dimension normal to the reflective boundary. (So that
 * the unrotated body is a mirror image of the unrotated mirror body)
 * b) reflecting the quaternion with qMirror so that the rotation is mirrored.
 *
 * Part (a) involves created a new molecule type for reflections in all 3 dimensions, so is not suitable for real simulations,
 * but we use it here so that we can compare results of the actually used algorithm with a fairly different algorithm.
 */
TEST_P(ReflectiveBoundaryConditionTest, simpleReflectionTest) {
  using namespace autopas::utils::ArrayMath::literals;

  // initialise AutoPas container & domainDecomposition
  MDFlexConfig config(0, nullptr);
  config.epsilonMap.value.clear();
  config.sigmaMap.value.clear();
  config.massMap.value.clear();

  const std::array<double, 3> boxMin = {0., 0., 0.};
  const std::array<double, 3> boxMax = {5., 5., 5.};

  config.boxMin.value = boxMin;
  config.boxMax.value = boxMax;
  const std::array<double, 3> boxLength = boxMax - boxMin;
  config.subdivideDimension.value = {true, true, true};
  const double cutoff = 0.3;
  config.cutoff.value = cutoff;
  config.verletSkinRadiusPerTimestep.value = 0.02;
  config.verletRebuildFrequency.value = 10;
  const double sigma = 1.;
  config.addSiteType(0, 1., sigma, 1.);
  config.boundaryOption.value = {options::BoundaryTypeOption::reflective, options::BoundaryTypeOption::reflective,
                                 options::BoundaryTypeOption::reflective};

  RegularGridDecomposition domainDecomposition(config);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);
  auto particlePropertiesLibrary = std::make_shared<ParticlePropertiesLibraryType>(cutoff);

  autoPasContainer->setBoxMin(boxMin);
  autoPasContainer->setBoxMax(boxMax);
  autoPasContainer->setCutoff(cutoff);
  autoPasContainer->setVerletSkinPerTimestep(config.verletSkinRadiusPerTimestep.value);
  autoPasContainer->setVerletRebuildFrequency(config.verletRebuildFrequency.value);
  autoPasContainer->init();

  particlePropertiesLibrary->addSiteType(0, 1., sigma, 1.);
#if MD_FLEXIBLE_MODE==MULTISITE
  // Correct Moment of Inertia is irrelevant to test, so accept that they are wrong
  particlePropertiesLibrary->addMolType(0, {0, 0, 0},
                                        {{0.074349607, 0.120300191, 0.}, {0.03249197, -0.137638192, 0.}, {-0.137638192, -0.03249197, 0.}}, {5.23606798, 0.76393202, 6.});
  particlePropertiesLibrary->addMolType(1, {0, 0, 0},
                                        {{-0.074349607, 0.120300191, 0.}, {-0.03249197, -0.137638192, 0.}, {0.137638192, -0.03249197, 0.}}, {5.23606798, 0.76393202, 6.});
  particlePropertiesLibrary->addMolType(2, {0, 0, 0},
                                        {{0.074349607, -0.120300191, 0.}, {0.03249197, 0.137638192, 0.}, {-0.137638192, 0.03249197, 0.}}, {5.23606798, 0.76393202, 6.});
  particlePropertiesLibrary->addMolType(3, {0, 0, 0},
                                        {{0.074349607, 0.120300191, -0.}, {0.03249197, -0.137638192, -0.}, {-0.137638192, -0.03249197, -0.}}, {5.23606798, 0.76393202, 6.});
#endif
  particlePropertiesLibrary->calculateMixingCoefficients();

  // get particle properties
  const std::array<double, 3> particlePosition = std::get<0>(GetParam());
#if MD_FLEXIBLE_MODE==MULTISITE
  const std::array<double, 4> particleQuaternion = std::get<1>(GetParam());
#endif

  const auto expectedPosition = particlePosition;
  std::array<double, 3> expectedForce{0., 0., 0.};
#if MD_FLEXIBLE_MODE==MULTISITE
  std::array<double, 3> expectedTorque{0., 0., 0.};
#endif

  // derive expected position
  auto addForceFromReflection = [&](const int dimensionOfBoundary, const bool isUpper) {
#if MD_FLEXIBLE_MODE==MULTISITE
    // get properties of mirror particle
    const auto mirroredPosition = [&] () {
      const auto distanceCenterOfMassToBoundary = isUpper ? boxMax[dimensionOfBoundary] - particlePosition[dimensionOfBoundary]
                                                          : particlePosition[dimensionOfBoundary] - boxMin[dimensionOfBoundary];
      auto mirroredPositionTmp = particlePosition;
      mirroredPositionTmp[dimensionOfBoundary] = isUpper ? boxMax[dimensionOfBoundary] + distanceCenterOfMassToBoundary
                                                         : boxMin[dimensionOfBoundary] - distanceCenterOfMassToBoundary;
      return mirroredPositionTmp;
    } ();
    const auto mirroredQuaternion = autopas::utils::quaternion::qMirror(particleQuaternion, dimensionOfBoundary);

    // get rotated site positions
    const auto unrotatedUntranslatedSitePositions = particlePropertiesLibrary->getSitePositions(0);
    const auto rotatedUntranslatedSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(particleQuaternion, unrotatedUntranslatedSitePositions);
    const auto unrotatedUntranslatedMirrorSitePositions = particlePropertiesLibrary->getSitePositions(1 + dimensionOfBoundary);
    const auto rotatedUntranslatedMirroredSitePositions = autopas::utils::quaternion::rotateVectorOfPositions(mirroredQuaternion, unrotatedUntranslatedMirrorSitePositions);

    // get expected force
    for (int siteOriginal = 0; siteOriginal < particlePropertiesLibrary->getNumSites(0); siteOriginal++) {
      for (int siteMirror = 0; siteMirror < particlePropertiesLibrary->getNumSites(0); siteMirror++) {
        const auto exactSitePosition =
            autopas::utils::ArrayMath::add(particlePosition, rotatedUntranslatedSitePositions[siteOriginal]);
        const auto exactMirrorSitePosition =
            autopas::utils::ArrayMath::add(mirroredPosition, rotatedUntranslatedMirroredSitePositions[siteMirror]);

        const auto displacementToMirrorParticle = autopas::utils::ArrayMath::sub(exactSitePosition, exactMirrorSitePosition);
        const auto distanceSquared = autopas::utils::ArrayMath::dot(displacementToMirrorParticle, displacementToMirrorParticle);

        const auto inverseDistanceSquared = 1. / distanceSquared;
        const auto lj2 = sigma * sigma * inverseDistanceSquared;
        const auto lj6 = lj2 * lj2 * lj2;
        const auto lj12 = lj6 * lj6;
        const auto lj12m6 = lj12 - lj6;
        const auto ljForceFactor = 24 * (lj12 + lj12m6) * inverseDistanceSquared;
        const auto forceContribution = autopas::utils::ArrayMath::mulScalar(displacementToMirrorParticle, ljForceFactor);
        const auto torqueContribution = autopas::utils::ArrayMath::cross(rotatedUntranslatedSitePositions[siteOriginal], forceContribution);
        expectedForce = autopas::utils::ArrayMath::add(expectedForce, forceContribution);
        expectedTorque = autopas::utils::ArrayMath::add(expectedTorque, torqueContribution);
      }
    }

    // if expectedForce is attractive towards the wall, reset everything
    if (expectedForce[dimensionOfBoundary] * (isUpper? 1 : -1) > 0) {
      expectedForce = {0., 0., 0.};
      expectedTorque = {0., 0., 0.};
    }
#else
    const auto distanceToBoundary = isUpper ? boxMax[dimensionOfBoundary] - particlePosition[dimensionOfBoundary]
                                            : particlePosition[dimensionOfBoundary] - boxMin[dimensionOfBoundary];
    const auto distanceToMirrorParticle = distanceToBoundary * 2.;
    const auto distanceSquared = distanceToMirrorParticle * distanceToMirrorParticle;

    const auto inverseDistanceSquared = 1. / distanceSquared;
    const auto lj2 = sigma * sigma * inverseDistanceSquared;
    const auto lj6 = lj2 * lj2 * lj2;
    const auto lj12 = lj6 * lj6;
    const auto lj12m6 = lj12 - lj6;
    const auto ljForceFactor = 24 * (lj12 + lj12m6) * inverseDistanceSquared;
    const auto force = ljForceFactor * distanceToMirrorParticle * (isUpper ? -1. : 1.);

    expectedForce[dimensionOfBoundary] += force;
#endif
  };

  for (int dimension = 0; dimension < 3; ++dimension) {
    if (particlePosition[dimension] < boxMin[dimension] + sixthRootOfTwo * sigma) {
      addForceFromReflection(dimension, false);
    } else if (particlePosition[dimension] > boxMax[dimension] - sixthRootOfTwo * sigma) {
      addForceFromReflection(dimension, true);
    }
  }

  // create particle and add to container
  // in the MPI case this is expected to be wrong for all but one rank
  if (domainDecomposition.isInsideLocalDomain(particlePosition)) {
    ParticleType particle;
    particle.setID(0);
    particle.setR(particlePosition);
    particle.setV({0., 0., 0.});
    particle.setF({0., 0., 0.});
#if MD_FLEXIBLE_MODE==MULTISITE
    particle.setQ(particleQuaternion);
    particle.setAngularVel({0., 0., 0.});
    particle.setTorque({0., 0., 0.});
#endif
    autoPasContainer->addParticle(particle);
#if not defined(AUTOPAS_INCLUDE_MPI)
  } else {
    GTEST_FAIL() << "Test Particle is not in the box -> setup is wrong!";
#endif
  }

  auto emigrants = autoPasContainer->updateContainer();

  // apply BCs + domain exchange
  domainDecomposition.exchangeMigratingParticles(*autoPasContainer, emigrants);
  domainDecomposition.reflectParticlesAtBoundaries(*autoPasContainer, *particlePropertiesLibrary);
  domainDecomposition.exchangeHaloParticles(*autoPasContainer);

  if (domainDecomposition.isInsideLocalDomain(particlePosition)) {
    EXPECT_EQ(1, autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned));
    // check particles have been successfully reflected (and not translated)
    const auto &reflectedParticle = autoPasContainer->begin(autopas::IteratorBehavior::owned);
    const auto &reflectedPosition = reflectedParticle->getR();
    const auto &reflectedVelocity = reflectedParticle->getV();
    const auto &reflectedForce = reflectedParticle->getF();
#if MD_FLEXIBLE_MODE==MULTISITE
    const auto &reflectedTorque = reflectedParticle->getTorque();
#endif
    for (size_t i = 0; i < 3; ++i) {
      EXPECT_NEAR(reflectedPosition[i], expectedPosition[i], 1e-13) << "Unexpected position[" << i << "]";
      EXPECT_NEAR(reflectedForce[i], expectedForce[i], 1e-13) << "Unexpected force[" << i << "]";
#if MD_FLEXIBLE_MODE==MULTISITE
      EXPECT_NEAR(reflectedTorque[i], expectedTorque[i], 1e-13) << "Unexpected torque[" << i << "]";
#endif
    }
#if not defined(AUTOPAS_INCLUDE_MPI)
  } else {
    GTEST_FAIL() << "Expected position is not in the box -> setup is wrong!";
  }
  // check that there are no halo particles. We can only guarantee this in the non-MPI case
  EXPECT_EQ(0, autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo));
#else
  }
#endif
}

/**
 * Tests reflective boundaries in all three dimensions for both single and multi site molecules.
 *
 * Expected results:
 * Particle 0: Repulsed in x (lower)
 * Particle 1: Repulsed in x (upper)
 * Particle 2,3: Same but y
 * Particle 4,5: Same but z
 * Particles 6-11: Same as 0-5 but with a non identity quaternion. Does not run with multi-site compilation.
 */
INSTANTIATE_TEST_SUITE_P(
    TestSimpleReflections, ReflectiveBoundaryConditionTest,
    testing::Values(/*position*/ /*velocity*/ /*quaternion*/
#if MD_FLEXIBLE_MODE==MULTISITE
                    std::make_tuple(std::array<double, 3>{0.5, 2.50, 2.50}, std::array<double, 4>{1., 0., 0., 0.}),
                    std::make_tuple(std::array<double, 3>{4.5, 2.50, 2.50}, std::array<double, 4>{1., 0., 0., 0.}),

                    std::make_tuple(std::array<double, 3>{2.50, 0.5, 2.50}, std::array<double, 4>{1., 0., 0., 0.}),
                    std::make_tuple(std::array<double, 3>{2.50, 4.5, 2.50}, std::array<double, 4>{1., 0., 0., 0.}),

                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 0.5}, std::array<double, 4>{1., 0., 0., 0.}),
                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 4.5}, std::array<double, 4>{1., 0., 0., 0.}),


                    std::make_tuple(std::array<double, 3>{0.5, 2.50, 2.50}, autopas::utils::ArrayMath::normalize(std::array<double, 4>{1., 0.5, 0.25, 0.125})),
                    std::make_tuple(std::array<double, 3>{4.5, 2.50, 2.50}, autopas::utils::ArrayMath::normalize(std::array<double, 4>{1., 0.5, 0.25, 0.125})),

                    std::make_tuple(std::array<double, 3>{2.50, 0.5, 2.50}, autopas::utils::ArrayMath::normalize(std::array<double, 4>{1., 0.5, 0.25, 0.125})),
                    std::make_tuple(std::array<double, 3>{2.50, 4.5, 2.50}, autopas::utils::ArrayMath::normalize(std::array<double, 4>{1., 0.5, 0.25, 0.125})),

                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 0.5}, autopas::utils::ArrayMath::normalize(std::array<double, 4>{1., 0.5, 0.25, 0.125})),
                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 4.5}, autopas::utils::ArrayMath::normalize(std::array<double, 4>{1., 0.5, 0.25, 0.125})))
#else
                    std::make_tuple(std::array<double, 3>{0.5, 2.50, 2.50}, std::array<double, 4>{1., 0., 0., 0.}),
                    std::make_tuple(std::array<double, 3>{4.5, 2.50, 2.50}, std::array<double, 4>{1., 0., 0., 0.}),

                    std::make_tuple(std::array<double, 3>{2.50, 0.5, 2.50}, std::array<double, 4>{1., 0., 0., 0.}),
                    std::make_tuple(std::array<double, 3>{2.50, 4.5, 2.50}, std::array<double, 4>{1., 0., 0., 0.}),

                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 0.5}, std::array<double, 4>{1., 0., 0., 0.}),
                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 4.5}, std::array<double, 4>{1., 0., 0., 0.}))
#endif
    //    ,ReflectiveBoundaryConditionTest::PrintToStringParamName());
);

/**
 * Implements the reflective boundary zoning test.
 * @param particlePosition
 * @param particleType Must be either 0 or 1. 0 corresponds to a sigma of 0.1, 1 corresponds to a sigma of 0.2.
 */
void testReflectiveBoundaryZoning(const std::array<double, 3> &particlePosition, int particleTypeID) {
  using namespace autopas::utils::ArrayMath::literals;
  if (particleTypeID != 0 and particleTypeID != 1) {
    std::cerr << "testReflectiveBoundaryZoning only takes particle types of 0 or 1 only!";
  }
  MDFlexConfig config(0, nullptr);
  config.epsilonMap.value.clear();
  config.sigmaMap.value.clear();
  config.massMap.value.clear();
#if MD_FLEXIBLE_MODE==MULTISITE
  config.molToSiteIdMap.clear();
  config.molToSitePosMap.clear();
  config.momentOfInertiaMap.clear();
#endif

  const std::array<double, 3> boxMin = {0., 0., 0.};
  const std::array<double, 3> boxMax = {5., 5., 5.};
  const std::array<double, 2> sigmas = {0.1, 0.2};
  const double cutoff = 0.3;

  config.boxMin.value = boxMin;
  config.boxMax.value = boxMax;
  const std::array<double, 3> boxLength = boxMax - boxMin;
  config.subdivideDimension.value = {true, true, true};
  config.cutoff.value = cutoff;
  config.verletSkinRadiusPerTimestep.value = 0.01;
  config.verletRebuildFrequency.value = 10;
  config.addSiteType(0, 1., sigmas[0], 1.);
  config.addSiteType(1, 1., sigmas[1], 1.);
  config.boundaryOption.value = {options::BoundaryTypeOption::reflective, options::BoundaryTypeOption::reflective,
                                 options::BoundaryTypeOption::reflective};

  RegularGridDecomposition domainDecomposition(config);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);
  auto particlePropertiesLibrary = std::make_shared<ParticlePropertiesLibraryType>(cutoff);

  autoPasContainer->setBoxMin(boxMin);
  autoPasContainer->setBoxMax(boxMax);
  autoPasContainer->setCutoff(cutoff);
  autoPasContainer->setVerletSkinPerTimestep(config.verletSkinRadiusPerTimestep.value);
  autoPasContainer->setVerletRebuildFrequency(config.verletRebuildFrequency.value);
  autoPasContainer->init();

  particlePropertiesLibrary->addSiteType(0, 1., sigmas[0], 1.);
  particlePropertiesLibrary->addSiteType(1, 1., sigmas[1], 1.);
#if MD_FLEXIBLE_MODE==MULTISITE
  particlePropertiesLibrary->addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.});
  particlePropertiesLibrary->addMolType(1, {1}, {{0., 0., 0.}}, {1., 1., 1.});
#endif
  particlePropertiesLibrary->calculateMixingCoefficients();

  std::array<bool, 3> expectReflection = {false, false, false};
  for (int dim = 0; dim < 3; dim++) {
    if (particlePosition[dim] < boxMin[dim] + sixthRootOfTwo * sigmas[particleTypeID] / 2.) {
      expectReflection[dim] = true;
    } else if (particlePosition[dim] > boxMax[dim] - sixthRootOfTwo * sigmas[particleTypeID] / 2.) {
      expectReflection[dim] = true;
    }
  }

  // create particle and add to container
  // in the MPI case this is expected to be wrong for all but one rank
  if (domainDecomposition.isInsideLocalDomain(particlePosition)) {
    ParticleType particle;
    particle.setID(0);
    particle.setR(particlePosition);
    particle.setF({0., 0., 0.});
#if MD_FLEXIBLE_MODE==MULTISITE
    particle.setQ({0., 0., 0., 1.});
    particle.setTorque({0., 0., 0.});
#endif
    particle.setTypeId(particleTypeID);
    autoPasContainer->addParticle(particle);
#if not defined(AUTOPAS_INCLUDE_MPI)
  } else {
    GTEST_FAIL() << "Test Particle is not in the box -> setup is wrong!";
#endif
  }

  domainDecomposition.reflectParticlesAtBoundaries(*autoPasContainer, *particlePropertiesLibrary);

  if (domainDecomposition.isInsideLocalDomain(particlePosition)) {
    auto returnedParticle = autoPasContainer->begin();

    const auto reflectedForce = returnedParticle->getF();
    // torque is irrelevant for single-site molecules.

    for (int dim = 0; dim < 3; dim++) {
      if (expectReflection[dim]) {
        EXPECT_NE(reflectedForce[dim], 0.)
            << "Particle does not experience reflective force in the " << dim << " dimension when it should.\n"
            << "Position = " << autopas::utils::ArrayUtils::to_string(particlePosition)
            << "; Actual Force = " << autopas::utils::ArrayUtils::to_string(reflectedForce) << ";";
      } else {
        EXPECT_DOUBLE_EQ(reflectedForce[dim], 0.)
            << "Particle experiences reflective force in the " << dim << " dimension when it shouldn't.\n"
            << "Position = " << autopas::utils::ArrayUtils::to_string(particlePosition)
            << "; Actual Force = " << autopas::utils::ArrayUtils::to_string(reflectedForce) << ";";
      }
    }
  }
}

/**
 * Tests that reflective boundaries are applied only to single-site particles near enough to the boundary that a repulsive
 * force is applied. For multi-site compilations, molecules consist solely of single sites.
 *
 * Four particles are created, with two different sigmas, as follows
 * * One, with a small sigma, is placed close to the boundary and is expected to be reflected.
 * * One, with a large sigma, is placed further from the boundary, but such that it is still expected to be reflected.
 * * One, with a small sigma, is placed at the same distance as the previous, but is not expected to be reflected.
 * * One, with a large sigma, is placed such that it is not expected to be reflected.
 *
 * This is repeated for all boundaries.
 */
TEST_F(ReflectiveBoundaryConditionTest, reflectiveSingleSiteZoningTest) {
  testReflectiveBoundaryZoning({0.05, 2.5, 2.5}, 0);
  testReflectiveBoundaryZoning({2.5, 0.05, 2.5}, 0);
  testReflectiveBoundaryZoning({2.5, 2.5, 0.05}, 0);
  testReflectiveBoundaryZoning({4.95, 2.5, 2.5}, 0);
  testReflectiveBoundaryZoning({2.5, 4.95, 2.5}, 0);
  testReflectiveBoundaryZoning({2.5, 2.5, 4.95}, 0);

  testReflectiveBoundaryZoning({0.1, 2.5, 2.5}, 1);
  testReflectiveBoundaryZoning({2.5, 0.1, 2.5}, 1);
  testReflectiveBoundaryZoning({2.5, 2.5, 0.1}, 1);
  testReflectiveBoundaryZoning({4.9, 2.5, 2.5}, 1);
  testReflectiveBoundaryZoning({2.5, 4.9, 2.5}, 1);
  testReflectiveBoundaryZoning({2.5, 2.5, 4.9}, 1);

  testReflectiveBoundaryZoning({0.1, 2.5, 2.5}, 0);
  testReflectiveBoundaryZoning({2.5, 0.1, 2.5}, 0);
  testReflectiveBoundaryZoning({2.5, 2.5, 0.1}, 0);
  testReflectiveBoundaryZoning({4.9, 2.5, 2.5}, 0);
  testReflectiveBoundaryZoning({2.5, 4.9, 2.5}, 0);
  testReflectiveBoundaryZoning({2.5, 2.5, 4.9}, 0);

  testReflectiveBoundaryZoning({0.15, 2.5, 2.5}, 1);
  testReflectiveBoundaryZoning({2.5, 0.15, 2.5}, 1);
  testReflectiveBoundaryZoning({2.5, 2.5, 0.15}, 1);
  testReflectiveBoundaryZoning({4.85, 2.5, 2.5}, 1);
  testReflectiveBoundaryZoning({2.5, 4.85, 2.5}, 1);
  testReflectiveBoundaryZoning({2.5, 2.5, 4.85}, 1);
}

/**
 * Places 4 molecules near a boundary:
 * * One where all sites are near enough that to experience repulsion. Expects repulsion.
 * * One where a site will be attracted to the boundary, but that the molecule overall experiences repulsion. Expects repulsion.
 * * One where a site will be repulsed from the boundary, but that the molecule overall experiences attraction. Expects no change in force on molecule.
 * * One with all sites far enough that they would all experience attraction. Expects no change in force on molecule.
 */
TEST_F(ReflectiveBoundaryConditionTest, reflectiveMultiSiteZoningTest) {
#ifndef MD_FLEXIBLE_USE_MULTI_SITE
  GTEST_SKIP() << "reflectiveMultiSiteZoningTest: Skipping as multi-site not compiled";
#else
  MDFlexConfig config(0, nullptr);
  config.epsilonMap.value.clear();
  config.sigmaMap.value.clear();
  config.massMap.value.clear();
  config.molToSiteIdMap.clear();
  config.molToSitePosMap.clear();
  config.momentOfInertiaMap.clear();

  const std::array<double, 3> boxMin = {0., 0., 0.};
  const std::array<double, 3> boxMax = {5., 5., 5.};
  const double cutoff = 0.3;

  config.boxMin.value = boxMin;
  config.boxMax.value = boxMax;
  const std::array<double, 3> boxLength = autopas::utils::ArrayMath::sub(boxMax, boxMin);
  config.subdivideDimension.value = {true, true, true};
  config.cutoff.value = cutoff;
  config.verletSkinRadiusPerTimestep.value = 0.01;
  config.verletRebuildFrequency.value = 10;
  config.addSiteType(0, 0.1, 0.2, 1.);
  config.addSiteType(1, 1000., 0.4, 1.);
  config.addSiteType(2, 0.0001, 0.4, 1.);
  config.boundaryOption.value = {options::BoundaryTypeOption::reflective, options::BoundaryTypeOption::reflective,
                                 options::BoundaryTypeOption::reflective};

  RegularGridDecomposition domainDecomposition(config);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);
  auto particlePropertiesLibrary = std::make_shared<ParticlePropertiesLibraryType>(cutoff);

  autoPasContainer->setBoxMin(boxMin);
  autoPasContainer->setBoxMax(boxMax);
  autoPasContainer->setCutoff(cutoff);
  autoPasContainer->setVerletSkinPerTimestep(config.verletSkinRadiusPerTimestep.value);
  autoPasContainer->setVerletRebuildFrequency(config.verletRebuildFrequency.value);
  autoPasContainer->init();

  particlePropertiesLibrary->addSiteType(0, 1, 0.2, 1.);
  particlePropertiesLibrary->addSiteType(1, 1000., 0.4, 1.);
  particlePropertiesLibrary->addSiteType(2, 0.0001, 0.4, 1.);
  particlePropertiesLibrary->calculateMixingCoefficients();

  // Add molecules

  // Molecule 0: All sites experience repulsion.
  particlePropertiesLibrary->addMolType(0, {0, 1}, {{-0.01, 0., 0.}, {0.01, 0., 0.}}, {1., 1., 1.});
  autopas::MultisiteMoleculeLJ mol0({0.05, 2.5, 2.5}, {0., 0., 0.}, {1., 0., 0., 0.}, {0., 0., 0.}, 0, 0);
  autoPasContainer->addParticle(mol0);

  // Molecule 1: One site experiences attraction, one site experiences repulsion, overall molecule experiences repulsion.
  particlePropertiesLibrary->addMolType(1, {0, 1}, {{-0.01, 0., 0.}, {0.01, 0., 0.}}, {1., 1., 1.});
  autopas::MultisiteMoleculeLJ mol1({0.2, 2.5, 2.5}, {0., 0., 0.}, {1., 0., 0., 0.}, {0., 0., 0.}, 1, 1);
  autoPasContainer->addParticle(mol1);

  // Molecule 2: Both sites experience attraction.
  particlePropertiesLibrary->addMolType(2, {2, 0}, {{-0.01, 0., 0.}, {0.01, 0., 0.}}, {1., 1., 1.});
  autopas::MultisiteMoleculeLJ mol2({0.22, 2.5, 2.5}, {0., 0., 0.}, {1., 0., 0., 0.}, {0., 0., 0.}, 2, 2);
  autoPasContainer->addParticle(mol2);

  // Molecule 3: One site experiences attraction, one site experiences repulsion, overall molecule experiences attraction.
  particlePropertiesLibrary->addMolType(3, {0, 01}, {{-0.01, 0., 0.}, {0.01, 0., 0.}}, {1., 1., 1.});
  autopas::MultisiteMoleculeLJ mol3({0.3, 2.5, 2.5}, {0., 0., 0.}, {1., 0., 0., 0.}, {0., 0., 0.}, 3, 3);
  autoPasContainer->addParticle(mol3);

  // Reflect

  domainDecomposition.reflectParticlesAtBoundaries(*autoPasContainer, *particlePropertiesLibrary);

  // Test for correctness

  for (auto iter = autoPasContainer->begin(); iter.isValid(); ++iter) {
    switch (iter->getID()) {
      case 0:
      case 1:
        // Expect repulsion. Reflection on lower boundary is a positive force.
        EXPECT_TRUE(iter->getF()[0] > 0) << "Molecule " << iter->getID() << " is not experiencing repulsion";
        break;
      case 2:
      case 3:
        EXPECT_DOUBLE_EQ(iter->getF()[0], 0) << "Molecule " << iter->getID() << " is experiencing a non-zero force";
        break;
    }
  }


#endif


}