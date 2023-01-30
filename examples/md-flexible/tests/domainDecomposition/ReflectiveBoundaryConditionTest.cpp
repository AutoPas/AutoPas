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
 * the particle receives the correct force.
 */
TEST_P(ReflectiveBoundaryConditionTest, simpleReflectionTest) {
  // initialise AutoPas container & domainDecomposition
  MDFlexConfig config(0, nullptr);
  config.epsilonMap.value.clear();
  config.sigmaMap.value.clear();
  config.massMap.value.clear();

  const std::array<double, 3> boxMin = {0., 0., 0.};
  const std::array<double, 3> boxMax = {5., 5., 5.};

  config.boxMin.value = boxMin;
  config.boxMax.value = boxMax;
  const std::array<double, 3> boxLength = autopas::utils::ArrayMath::sub(boxMax, boxMin);
  config.subdivideDimension.value = {true, true, true};
  const double cutoff = 0.3;
  config.cutoff.value = cutoff;
  config.verletSkinRadiusPerTimestep.value = 0.02;
  config.verletRebuildFrequency.value = 10;
  const double sigma = 1.;
  config.addParticleType(0, 1., sigma, 1.);
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

  particlePropertiesLibrary->addType(0, 1., sigma, 1.);
  particlePropertiesLibrary->calculateMixingCoefficients();

  // get particle properties
  const std::array<double, 3> particlePosition = std::get<0>(GetParam());
  const std::array<double, 3> particleVelocity = std::get<1>(GetParam());

  // derive expected position
  auto forceFromReflection = [&](const std::array<double, 3> position, const int dimensionOfBoundary,
                                 const bool isUpper) {
    const auto distanceToBoundary = isUpper ? boxMax[dimensionOfBoundary] - position[dimensionOfBoundary]
                                            : position[dimensionOfBoundary] - boxMin[dimensionOfBoundary];
    const auto distanceToMirrorParticle = distanceToBoundary * 2.;
    const auto distanceSquared = distanceToMirrorParticle * distanceToMirrorParticle;

    const auto inverseDistanceSquared = 1. / distanceSquared;
    const auto lj2 = sigma * sigma * inverseDistanceSquared;
    const auto lj6 = lj2 * lj2 * lj2;
    const auto lj12 = lj6 * lj6;
    const auto lj12m6 = lj12 - lj6;
    const auto ljForceFactor = 24 * (lj12 + lj12m6) * inverseDistanceSquared;
    const auto force = ljForceFactor * distanceToMirrorParticle * (isUpper ? -1. : 1.);

    return force;
  };

  const auto expectedPosition = particlePosition;
  const auto expectedVelocity = particleVelocity;
  std::array<double, 3> expectedForce = {0., 0., 0.};
  for (int dimension = 0; dimension < 3; ++dimension) {
    if (particlePosition[dimension] < boxMin[dimension] + sixthRootOfTwo * sigma) {
      expectedForce[dimension] = forceFromReflection(particlePosition, dimension, false);
    } else if (particlePosition[dimension] > boxMax[dimension] - sixthRootOfTwo * sigma) {
      expectedForce[dimension] = forceFromReflection(particlePosition, dimension, true);
    }
  }

  // create particle and add to container
  // in the MPI case this is expected to be wrong for all but one rank
  if (domainDecomposition.isInsideLocalDomain(particlePosition)) {
    ParticleType particle;
    particle.setID(0);
    particle.setR(particlePosition);
    particle.setV(particleVelocity);
    particle.setF({0., 0., 0.});
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
    for (size_t i = 0; i < 3; ++i) {
      EXPECT_NEAR(reflectedPosition[i], expectedPosition[i], 1e-13) << "Unexpected position[" << i << "]";
      EXPECT_NEAR(reflectedVelocity[i], expectedVelocity[i], 1e-13) << "Unexpected velocity[" << i << "]";
      EXPECT_NEAR(reflectedForce[i], expectedForce[i], 1e-13) << "Unexpected force[" << i << "]";
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

INSTANTIATE_TEST_SUITE_P(
    TestSimpleReflections, ReflectiveBoundaryConditionTest,
    testing::Values(/*position*/ /*velocity*/ /*is reflected*/
                    std::make_tuple(std::array<double, 3>{0.005, 2.50, 2.50}, std::array<double, 3>{1, 1, -1},
                                    std::array<bool, 3>{false, false, false}),
                    std::make_tuple(std::array<double, 3>{0.005, 2.50, 2.50}, std::array<double, 3>{-1, 1, -1},
                                    std::array<bool, 3>{true, false, false}),
                    std::make_tuple(std::array<double, 3>{4.995, 2.50, 2.50}, std::array<double, 3>{1, 1, -1},
                                    std::array<bool, 3>{true, false, false}),
                    std::make_tuple(std::array<double, 3>{4.995, 2.50, 2.50}, std::array<double, 3>{-1, 1, -1},
                                    std::array<bool, 3>{false, false, false}),

                    std::make_tuple(std::array<double, 3>{2.50, 0.005, 2.50}, std::array<double, 3>{1, 1, -1},
                                    std::array<bool, 3>{false, false, false}),
                    std::make_tuple(std::array<double, 3>{2.50, 0.005, 2.50}, std::array<double, 3>{1, -1, -1},
                                    std::array<bool, 3>{false, true, false}),
                    std::make_tuple(std::array<double, 3>{2.50, 4.995, 2.50}, std::array<double, 3>{1, 1, -1},
                                    std::array<bool, 3>{false, true, false}),
                    std::make_tuple(std::array<double, 3>{2.50, 4.995, 2.50}, std::array<double, 3>{1, -1, -1},
                                    std::array<bool, 3>{false, false, false}),

                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 0.005}, std::array<double, 3>{1, -1, 1},
                                    std::array<bool, 3>{false, false, false}),
                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 0.005}, std::array<double, 3>{1, -1, -1},
                                    std::array<bool, 3>{false, false, true}),
                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 4.995}, std::array<double, 3>{1, -1, 1},
                                    std::array<bool, 3>{false, false, true}),
                    std::make_tuple(std::array<double, 3>{2.50, 2.50, 4.995}, std::array<double, 3>{1, -1, -1},
                                    std::array<bool, 3>{false, false, false}))
    //    ,ReflectiveBoundaryConditionTest::PrintToStringParamName());
);

/**
 * Implements the reflective boundary zoning test.
 * @param particlePosition
 * @param particleType Must be either 0 or 1. 0 corresponds to a sigma of 0.5, 1 corresponds to a sigma of 1.0.
 */
void testReflectiveBoundaryZoning(const std::array<double, 3> particlePosition, int particleTypeID) {
  if (particleTypeID != 0 and particleTypeID != 1) {
    std::cerr << "testReflectiveBoundaryZoning only takes particle types of 0 or 1 only!";
  }
  MDFlexConfig config(0, nullptr);
  config.epsilonMap.value.clear();
  config.sigmaMap.value.clear();
  config.massMap.value.clear();

  const std::array<double, 3> boxMin = {0., 0., 0.};
  const std::array<double, 3> boxMax = {5., 5., 5.};
  const std::array<double, 2> sigmas = {0.1, 0.2};
  const double cutoff = 0.3;

  config.boxMin.value = boxMin;
  config.boxMax.value = boxMax;
  const std::array<double, 3> boxLength = autopas::utils::ArrayMath::sub(boxMax, boxMin);
  config.subdivideDimension.value = {true, true, true};
  config.cutoff.value = cutoff;
  config.verletSkinRadiusPerTimestep.value = 0.01;
  config.verletRebuildFrequency.value = 10;
  config.addParticleType(0, 1., sigmas[0], 1.);
  config.addParticleType(1, 1., sigmas[1], 1.);
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

  particlePropertiesLibrary->addType(0, 1., sigmas[0], 1.);
  particlePropertiesLibrary->addType(1, 1., sigmas[1], 1.);
  particlePropertiesLibrary->calculateMixingCoefficients();

  std::array<bool, 3> expectReflection = {false, false, false};
  for (int dim = 0; dim < 3; dim++) {
    if (particlePosition[dim] < boxMin[dim] + sixthRootOfTwo * sigmas[particleTypeID] / 2.) {
      expectReflection[dim] = true;
    } else if (particlePosition[dim] > boxMax[dim] - sixthRootOfTwo * sigmas[particleTypeID] / 2.) {
      expectReflection[dim] = true;
    }
  }

  ParticleType particle;
  particle.setR(particlePosition);
  particle.setF({0., 0., 0.});
  particle.setTypeId(particleTypeID);

  autoPasContainer->addParticle(particle);

  domainDecomposition.reflectParticlesAtBoundaries(*autoPasContainer, *particlePropertiesLibrary);

  auto returnedParticle = autoPasContainer->begin();

  const auto reflectedForce = returnedParticle->getF();

  for (int dim = 0; dim < 3; dim++) {
    if (expectReflection[dim]) {
      EXPECT_NE(reflectedForce[dim], 0.) << "Particle does not experience reflective force in the " << dim << " dimension when it should.\n"
          << "Position = " << autopas::utils::ArrayUtils::to_string(particlePosition) << "; Actual Force = " << autopas::utils::ArrayUtils::to_string(reflectedForce) << ";";
    } else {
      EXPECT_DOUBLE_EQ(reflectedForce[dim], 0.) << "Particle experiences reflective force in the " << dim << " dimension when it shouldn't.\n"
                                                << "Position = " << autopas::utils::ArrayUtils::to_string(particlePosition) << "; Actual Force = " << autopas::utils::ArrayUtils::to_string(reflectedForce) << ";";
    }
  }
}

/**
 * Tests that reflective boundaries are applied only to particles near enough to the boundary that a repulsive force is
 * applied.
 *
 * Four particles are created, with two different sigmas, as follows
 * * One, with a small sigma, is placed close to the boundary and is expected to be reflected.
 * * One, with a large sigma, is placed further from the boundary, but such that it is still expected to be reflected.
 * * One, with a small sigma, is placed at the same distance as the previous, but is not expected to be reflected.
 * * One, with a large sigma, is placed such that it is not expected to be reflected.
 *
 * This is repeated for all boundaries.
 */
TEST_F(ReflectiveBoundaryConditionTest, reflectiveZoningTest) {
  testReflectiveBoundaryZoning({0.05, 2.5,  2.5},  0);
  testReflectiveBoundaryZoning({2.5,  0.05, 2.5},  0);
  testReflectiveBoundaryZoning({2.5,  2.5,  0.05}, 0);
  testReflectiveBoundaryZoning({4.95, 2.5,  2.5},  0);
  testReflectiveBoundaryZoning({2.5,  4.95, 2.5},  0);
  testReflectiveBoundaryZoning({2.5,  2.5,  4.95}, 0);

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

  testReflectiveBoundaryZoning({0.15, 2.5,  2.5}, 1);
  testReflectiveBoundaryZoning({2.5,  0.15, 2.5}, 1);
  testReflectiveBoundaryZoning({2.5,  2.5,  0.15}, 1);
  testReflectiveBoundaryZoning({4.85, 2.5,  2.5}, 1);
  testReflectiveBoundaryZoning({2.5,  4.85, 2.5}, 1);
  testReflectiveBoundaryZoning({2.5,  2.5,  4.85}, 1);
}
