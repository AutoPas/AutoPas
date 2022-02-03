/**
 * @file MixedBoundaryConditionTest.cpp
 * @author S. J. Newcome
 * @date 21/01/2022
*/
#include "MixedBoundaryConditionTest.h"
#include "src/TypeDefinitions.h"
#include "src/configuration/MDFlexConfig.h"
#include "src/domainDecomposition/DomainTools.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"

namespace {
void initializeAutoPasContainer(RegularGridDecomposition::SharedAutoPasContainer &autoPasContainer,
                                MDFlexConfig &configuration) {
  autoPasContainer->setAllowedCellSizeFactors(*configuration.cellSizeFactors.value);
  autoPasContainer->setAllowedContainers(configuration.containerOptions.value);
  autoPasContainer->setAllowedDataLayouts(configuration.dataLayoutOptions.value);
  autoPasContainer->setAllowedNewton3Options(configuration.newton3Options.value);
  autoPasContainer->setAllowedTraversals(configuration.traversalOptions.value);
  autoPasContainer->setAllowedLoadEstimators(configuration.loadEstimatorOptions.value);
  autoPasContainer->setBoxMin(configuration.boxMin.value);
  autoPasContainer->setBoxMax(configuration.boxMax.value);
  autoPasContainer->setCutoff(configuration.cutoff.value);
  autoPasContainer->setRelativeOptimumRange(configuration.relativeOptimumRange.value);
  autoPasContainer->setMaxTuningPhasesWithoutTest(configuration.maxTuningPhasesWithoutTest.value);
  autoPasContainer->setRelativeBlacklistRange(configuration.relativeBlacklistRange.value);
  autoPasContainer->setEvidenceFirstPrediction(configuration.evidenceFirstPrediction.value);
  autoPasContainer->setExtrapolationMethodOption(configuration.extrapolationMethodOption.value);
  autoPasContainer->setNumSamples(configuration.tuningSamples.value);
  autoPasContainer->setMaxEvidence(configuration.tuningMaxEvidence.value);
  autoPasContainer->setSelectorStrategy(configuration.selectorStrategy.value);
  autoPasContainer->setTuningInterval(configuration.tuningInterval.value);
  autoPasContainer->setTuningStrategyOption(configuration.tuningStrategyOption.value);
  autoPasContainer->setMPIStrategy(configuration.mpiStrategyOption.value);
  autoPasContainer->setVerletClusterSize(configuration.verletClusterSize.value);
  autoPasContainer->setVerletRebuildFrequency(configuration.verletRebuildFrequency.value);
  autoPasContainer->setVerletSkin(configuration.verletSkinRadius.value);
  autoPasContainer->setAcquisitionFunction(configuration.acquisitionFunctionOption.value);
  autoPasContainer->init();
}
}  // namespace

/**
 * Simple test designed to show that reflection at a corner between 2 reflective faces works correctly.
 * Places particles within reflective skin and tests for correct reflective behaviour including for:
 * - cases where the particle's velocity is away from the boundary
 * - cases where the particle is in a corner, including
 *   + reflective/reflective
 *   + reflective/periodic
 *   + reflective/reflective/periodic
 * - cases for both upper and lower boundaries
 */
TEST_F(MixedBoundaryConditionTest, testSimpleReflection) {
  // initialise AutoPas container & domainDecomposition
  const std::array<double,3> boxMin = {0.,0.,0.};
  const std::array<double,3> boxMax = {5.,5.,5.};
  const std::array<double,3> boxLength = autopas::utils::ArrayMath::sub(boxMax, boxMin);
  const std::array<bool,3> subdivideDimension = {true,true,true};
  const double cutoffWidth = 2.;
  const double skinWidth = 0.2;
  const std::array<options::BoundaryTypeOption,3> boundaryConditions =
      {options::BoundaryTypeOption::periodic,options::BoundaryTypeOption::reflective,options::BoundaryTypeOption::reflective};

  RegularGridDecomposition domainDecomposition(boxMin, boxMax, subdivideDimension,
                                               cutoffWidth, skinWidth,boundaryConditions);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

  //initializeAutoPasContainer(autoPasContainer, configuration);
  autoPasContainer->setBoxMin(boxMin);
  autoPasContainer->setBoxMax(boxMax);
  autoPasContainer->setCutoff(cutoffWidth);
  autoPasContainer->setVerletSkin(skinWidth);
  autoPasContainer->init();

  const std::vector<std::array<double, 3>> particlePositions = {
      {1.0,0.05,1.0}, {2.0,0.05,1.0}, {1.0,0.05,0.05}, {2.0,0.05,0.05}, {4.95,0.05,0.05},
      {4.0,4.95,4.0}, {3.0,4.95,4.0}, {4.0,4.95,4.95}, {3.0,4.95,4.95}, {0.05,4.95,4.95}
  };
  const std::vector<std::array<double, 3>> particleVelocities = {
      {0.0,-1.0,1.0}, {0.0,1.0,1.0}, {0.0,-1.0,-0.5}, {0.0,-1.0,0.5}, {1.0,-0.5,-0.2},
      {0.0,1.0,-1.0}, {0.0,-1.0,-1.0}, {0.0,1.0,0.5}, {0.0,1.0,-0.5}, {-1.0,0.5,0.2}
  };

  size_t id = 0;
  for (const auto &position : particlePositions) {
    ParticleType particle;
    particle.setID(id);
    particle.setR(particlePositions[id]);
    particle.setV(particleVelocities[id]);

    autoPasContainer->addParticle(particle);

    ++id;
  }

  EXPECT_NO_THROW(domainDecomposition.reflectParticlesAtBoundaries(autoPasContainer));

  // derive the expected velocities by reflecting the velocities in the appropriate directions
  std::vector<std::array<double, 3>> expectedVelocities = particleVelocities;
  expectedVelocities[0][1] *= -1;
  // particle 1 should not reflect
  expectedVelocities[2][1] *= -1; expectedVelocities[2][2] *= -1;
  expectedVelocities[3][1] *= -1;
  expectedVelocities[4][1] *= -1; expectedVelocities[4][2] *= -1;
  expectedVelocities[5][1] *= -1;
  // particle 6 should not reflect
  expectedVelocities[7][1] *= -1; expectedVelocities[7][2] *= -1;
  expectedVelocities[8][1] *= -1;
  expectedVelocities[9][1] *= -1; expectedVelocities[9][2] *= -1;



  for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    const auto id = particle->getID();
    const auto &velocity = particle->getV();
    EXPECT_NEAR(velocity[0], expectedVelocities[id][0],1e-13);
    EXPECT_NEAR(velocity[1], expectedVelocities[id][1],1e-13);
    EXPECT_NEAR(velocity[2], expectedVelocities[id][2],1e-13);
  }

}

/**
 * Designed to test that exchangeMigratingParticles and exchangeHaloParticles in the mixed boundary case
 * Note: this is not designed to replace the more extensive tests in RegularGridDecompositionTest, but to test the
 * periodic BC in the mixed case.
 * Places particles in reflective skin and (primarily) outside of periodic boundary to test that particles
 * translated (as a result of the periodic boundary) and halo particles have the correct reflections.
 */
TEST_F(MixedBoundaryConditionTest, testPeriodic) {
  // initialise AutoPas container & domainDecomposition
  const std::array<double,3> boxMin = {0.,0.,0.};
  const std::array<double,3> boxMax = {5.,5.,5.};
  const std::array<double,3> boxLength = {boxMax[0] - boxMin[0],boxMax[1] - boxMin[1],boxMax[2] - boxMin[2]};
  const std::array<bool,3> subdivideDimension = {true,true,true};
  const double cutoffWidth = 2.;
  const double skinWidth = 0.2;
  const std::array<options::BoundaryTypeOption,3> boundaryConditions =
      {options::BoundaryTypeOption::periodic,options::BoundaryTypeOption::reflective,options::BoundaryTypeOption::reflective};

  RegularGridDecomposition domainDecomposition(boxMin, boxMax, subdivideDimension,
                                               cutoffWidth, skinWidth,boundaryConditions);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

  //initializeAutoPasContainer(autoPasContainer, configuration);
  autoPasContainer->setBoxMin(boxMin);
  autoPasContainer->setBoxMax(boxMax);
  autoPasContainer->setCutoff(cutoffWidth);
  autoPasContainer->setVerletSkin(skinWidth);
  autoPasContainer->init();

  // the last position is purposfully nonsense - it is just to confirm periodic BCs aren't being applied to reflective boundaries
  const std::vector<std::array<double, 3>> particlePositions = {
      {-0.05,0.05,1.0}, {-0.05,0.05,2.0}, {-0.05,0.05,3.0}, {4.95,4.95,4.95}, {5.05,4.95,4.95}, {5.05,5.05,4.95}
  };

  const std::vector<std::array<double, 3>> particleVelocities = {
      {-1.0,-1.0,0.0}, {-1.0,1.0,0.0}, {1.0,-1.0,0.0}, {1.0,1.0,1.0}, {1.0,1.0,1.0}, {1.0,1.0,1.0}
  };

  size_t id = 0;
  for (const auto position : particlePositions) {
    ParticleType particle;
    particle.setID(id);
    autoPasContainer->addParticle(particle);

    ++id;
  }

  for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    const auto id = particle->getID();
    particle->setR(particlePositions[id]);
    particle->setV(particleVelocities[id]);

  }

  EXPECT_NO_THROW(domainDecomposition.exchangeMigratingParticles(autoPasContainer));
  EXPECT_NO_THROW(domainDecomposition.reflectParticlesAtBoundaries(autoPasContainer));
  EXPECT_NO_THROW(domainDecomposition.exchangeHaloParticles(autoPasContainer));

//  const std::vector<std::array<double, 3>> expectedPositions = {
//      {4.95,0.05,1.0}, {4.95,0.05,2.0}, {4.95,0.05,3.0}, {4.95,4.95,4.95}, {0.05,4.95,4.95}, {0.05,5.05,4.95}
//  };
  // derive the expected positions
  auto expectedPositions = particlePositions;
  expectedPositions[0][0] += boxLength[0];
  expectedPositions[1][0] += boxLength[0];
  expectedPositions[2][0] += boxLength[0];
  // particle 3 is not translated across periodic boundary
  expectedPositions[4][0] -= boxLength[0];
  expectedPositions[5][0] -= boxLength[0];

  // derive the expected velocities
  auto expectedVelocities = particleVelocities;
  expectedVelocities[0][1] *= -1;
  // particle 1 is not reflected
  expectedVelocities[2][1] *= -1;
  expectedVelocities[3][1] *= -1; expectedVelocities[3][2] *= -1;
  expectedVelocities[4][1] *= -1; expectedVelocities[4][2] *= -1;
  expectedVelocities[5][2] *= -1;

  // derive expected halo particle positions
  // except for particle 3, all particles are beyond the periodic boundary, thus they are translated to the other side
  // of the box. Halo particles are created and translated back across the box to the original position. Thus, the
  // expected halo positions should map the original position.
  // Particle 3 is within the box, thus is not translated. A halo particle is still created by translating the position
  // to the other side of the box.
  auto expectedHaloPositions = particlePositions;
  expectedHaloPositions[3][0] -= boxLength[0];

  // derive expected halo particle velocities
  // these should match completely the velocities of original particles
  const auto expectedHaloVelocities = expectedVelocities;


  for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    const auto id = particle->getID();
    const auto position = particle->getR();
    const auto velocity = particle->getV();
    EXPECT_NEAR(position[0], expectedPositions[id][0],1e-13);
    EXPECT_NEAR(position[1], expectedPositions[id][1],1e-13);
    EXPECT_NEAR(position[2], expectedPositions[id][2],1e-13);
    EXPECT_NEAR(velocity[0], expectedVelocities[id][0],1e-13);
    EXPECT_NEAR(velocity[1], expectedVelocities[id][1],1e-13);
    EXPECT_NEAR(velocity[2], expectedVelocities[id][2],1e-13);
  }

  for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::halo); particle.isValid(); ++particle) {
    const auto id = particle->getID();
    const auto position = particle->getR();
    const auto velocity = particle->getV();
    EXPECT_NEAR(position[0], expectedHaloPositions[id][0],1e-13);
    EXPECT_NEAR(position[1], expectedHaloPositions[id][1],1e-13);
    EXPECT_NEAR(position[2], expectedHaloPositions[id][2],1e-13);
    EXPECT_NEAR(velocity[0], expectedHaloVelocities[id][0],1e-13);
    EXPECT_NEAR(velocity[1], expectedHaloVelocities[id][1],1e-13);
    EXPECT_NEAR(velocity[2], expectedHaloVelocities[id][2],1e-13);
  }
}

/**
 * Checks that the logic used to decide whether or not to implement a certain boundary condition is not only working in
 * the case where the boundaries are exclusively periodic or reflective. This should imply that there won't be any
 * problems with adding additional boundary condition types in the future.
 *
 * This test does this by confirming that no periodic or reflective behaviour occurs when no boundary conditions are
 * used. Two sets of particles are placed, one for the upper and lower boundaries, each of two particles - one just
 * inside the boundary, the other just outside - and nothing should change.
 */
TEST_F(MixedBoundaryConditionTest, testNoBoundary) {
  // initialise AutoPas container & domainDecomposition
  const std::array<double,3> boxMin = {0.,0.,0.};
  const std::array<double,3> boxMax = {5.,5.,5.};
  const std::array<bool,3> subdivideDimension = {true,true,true};
  const double cutoffWidth = 2.;
  const double skinWidth = 0.2;
  const std::array<options::BoundaryTypeOption,3> boundaryConditions =
      {options::BoundaryTypeOption::none,options::BoundaryTypeOption::none,options::BoundaryTypeOption::none};

  RegularGridDecomposition domainDecomposition(boxMin, boxMax, subdivideDimension,
                                               cutoffWidth, skinWidth,boundaryConditions);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

  //initializeAutoPasContainer(autoPasContainer, configuration);
  autoPasContainer->setBoxMin(boxMin);
  autoPasContainer->setBoxMax(boxMax);
  autoPasContainer->setCutoff(cutoffWidth);
  autoPasContainer->setVerletSkin(skinWidth);
  autoPasContainer->init();

  const std::vector<std::array<double, 3>> particlePositions = {
      {-0.05, 2.5, 2.5}, {0.05, 2.5, 2.5}, {4.95, 2.5, 2.5}, {5.05, 2.5, 2.5}
  };

  const std::vector<std::array<double, 3>> particleVelocities = {
      {-1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}
  };

  size_t id = 0;
  for (const auto position : particlePositions) {
    ParticleType particle;
    particle.setID(id);
    autoPasContainer->addParticle(particle);
    ++id;
  }

  for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    const auto id = particle->getID();
    particle->setR(particlePositions[id]);
    particle->setV(particleVelocities[id]);
  }

  EXPECT_NO_THROW(domainDecomposition.exchangeMigratingParticles(autoPasContainer));
  EXPECT_NO_THROW(domainDecomposition.reflectParticlesAtBoundaries(autoPasContainer));
  EXPECT_NO_THROW(domainDecomposition.exchangeHaloParticles(autoPasContainer));

  // check positions and velocities have not changed
  for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    const auto id = particle->getID();
    const auto position = particle->getR();
    const auto velocity = particle->getV();
    EXPECT_NEAR(position[0], particlePositions[id][0],1e-13);
    EXPECT_NEAR(position[1], particlePositions[id][1],1e-13);
    EXPECT_NEAR(position[2], particlePositions[id][2],1e-13);
    EXPECT_NEAR(velocity[0], particleVelocities[id][0],1e-13);
    EXPECT_NEAR(velocity[1], particleVelocities[id][1],1e-13);
    EXPECT_NEAR(velocity[2], particleVelocities[id][2],1e-13);
  }

  // check that there are no halo particles
  EXPECT_EQ(0, autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo));

}
