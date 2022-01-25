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
 * Simple test designed to show that reflection at a corner between 2 reflective faces works correctly
 */
TEST_F(MixedBoundaryConditionTest, testSimpleReflection) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "reflectionTest.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> localBoxMin = configuration.boxMin.value;
  std::array<double, 3> localBoxMax = configuration.boxMax.value;

  RegularGridDecomposition domainDecomposition(configuration.boxMin.value, configuration.boxMax.value,
                                               configuration.subdivideDimension.value, configuration.cutoff.value,
                                               configuration.verletSkinRadius.value, configuration.reflWidth.value,
                                               configuration.boundaryOption.value);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

  initializeAutoPasContainer(autoPasContainer, configuration);

  std::vector<std::vector<double>> particlePositions = {
      {1.0,0.05,1.0}, {2.0,0.05,1.0}, {1.0,0.05,0.05}, {2.0,0.05,0.05}, {4.95,0.05,0.05},
      {4.0,4.95,4.0}, {3.0,4.95,4.0}, {4.0,4.95,4.95}, {3.0,4.95,4.95}, {0.05,4.95,4.95}
  };
  std::vector<std::vector<double>> particleVelocities = {
      {0.0,-1.0,1.0}, {0.0,1.0,1.0}, {0.0,-1.0,-0.5}, {0.0,-1.0,0.5}, {1.0,-0.5,-0.2},
      {0.0,1.0,-1.0}, {0.0,-1.0,-1.0}, {0.0,1.0,0.5}, {0.0,1.0,-0.5}, {-1.0,0.5,0.2}
  };

  size_t id = 0;
  for (const auto position : particlePositions) {
    ParticleType particle;
    particle.setID(id);
    particle.setR({position[0], position[1], position[2]});
    particle.setV({particleVelocities[id][0], particleVelocities[id][1], particleVelocities[id][2]});

    autoPasContainer->addParticle(particle);

    ++id;
  }

  EXPECT_NO_THROW(domainDecomposition.reflectParticlesAtBoundaries(autoPasContainer));

  std::vector<std::array<double, 3>> expectedVelocities = {
      {0.0,1.0,1.0}, {0.0,1.0,1.0}, {0.0,1.0,0.5}, {0.0,1.0,0.5}, {1.0,0.5,0.2},
      {0.0,-1.0,-1.0}, {0.0,-1.0,-1.0}, {0.0,-1.0,-0.5}, {0.0,-1.0,-0.5}, {-1.0,-0.5,-0.2}
  };

  for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    const auto id = particle->getID();
    const auto velocity = particle->getV();
    EXPECT_NEAR(velocity[0], expectedVelocities[id][0],1e-13);
    EXPECT_NEAR(velocity[1], expectedVelocities[id][1],1e-13);
    EXPECT_NEAR(velocity[2], expectedVelocities[id][2],1e-13);
  }

}

/**
 * Designed to test that exchangeMigratingParticles and exchangeHaloParticles in the mixed boundary case
 * Note: this is not designed to replace the more extensive tests in RegularGridDecompositionTest, but to test the
 * periodic BC in the mixed case.
 */
TEST_F(MixedBoundaryConditionTest, testPeriodic) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                        std::string(YAMLDIRECTORY) + "reflectionTest.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);

  std::array<double, 3> localBoxMin = configuration.boxMin.value;
  std::array<double, 3> localBoxMax = configuration.boxMax.value;

  RegularGridDecomposition domainDecomposition(configuration.boxMin.value, configuration.boxMax.value,
                                               configuration.subdivideDimension.value, configuration.cutoff.value,
                                               configuration.verletSkinRadius.value, configuration.reflWidth.value,
                                               configuration.boundaryOption.value);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

  initializeAutoPasContainer(autoPasContainer, configuration);

  // the last position is purposfully nonsense - it is just to confirm periodic BCs aren't being applied to reflective boundaries
  std::vector<std::vector<double>> particlePositions = {
      {-0.05,0.05,1.0}, {-0.05,0.05,2.0}, {-0.05,0.05,3.0}, {4.95,4.95,4.95}, {5.05,4.95,4.95}, {5.05,5.05,4.95}
  };

  std::vector<std::vector<double>> particleVelocities = {
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
    particle->setR({particlePositions[id][0], particlePositions[id][1], particlePositions[id][2]});
    particle->setV({particleVelocities[id][0], particleVelocities[id][1], particleVelocities[id][2]});

  }

  EXPECT_NO_THROW(domainDecomposition.exchangeMigratingParticles(autoPasContainer));
  EXPECT_NO_THROW(domainDecomposition.reflectParticlesAtBoundaries(autoPasContainer));
  EXPECT_NO_THROW(domainDecomposition.exchangeHaloParticles(autoPasContainer));

  std::vector<std::vector<double>> expectedPositions = {
      {4.95,0.05,1.0}, {4.95,0.05,2.0}, {4.95,0.05,3.0}, {4.95,4.95,4.95}, {0.05,4.95,4.95}, {0.05,5.05,4.95}
  };
  std::vector<std::vector<double>> expectedVelocities = {
      {-1.0,1.0,0.0}, {-1.0,1.0,0.0}, {1.0,1.0,0.0}, {1.0,-1.0,-1.0}, {1.0,-1.0, -1.0}, {1.0,1.0,-1.0}
  };


  std::vector<std::vector<double>> expectedHaloPositions = {
      {-0.05,0.05,1.0}, {-0.05,0.05,2.0}, {-0.05,0.05,3.0}, {-0.05,4.95,4.95}, {5.05,4.95,4.95}, {5.05,5.05,4.95}
  };
  std::vector<std::vector<double>> expectedHaloVelocities = {
      {-1.0,1.0,0.0}, {-1.0,1.0,0.0}, {1.0,1.0,0.0}, {1.0,-1.0,-1.0}, {1.0,-1.0, -1.0}, {1.0,1.0,-1.0}
  };


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

