/**
 * @file TestRegularGridDecomposition.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "TestRegularGridDecomposition.h"

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapMPI.h"
#include "src/ParticleAttributes.h"
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

TEST_F(TestRegularGridDecomposition, testGetLocalDomain) {
  std::array<double, 3> globalBoxMin = {1.0, 1.0, 1.0};
  std::array<double, 3> globalBoxMax = {10.0, 10.0, 10.0};

  RegularGridDecomposition domainDecomposition(globalBoxMin, globalBoxMax, 0, 0);

  std::array<double, 3> globalBoxExtend = autopas::utils::ArrayMath::sub(globalBoxMax, globalBoxMin);

  int numberOfProcesses;
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numberOfProcesses);

  std::array<int, 3> decomposition;
  DomainTools::generateDecomposition(numberOfProcesses, decomposition);

  std::array<double, 3> expectedLocalBoxExtend = autopas::utils::ArrayMath::div(
      globalBoxExtend, {(double)decomposition[0], (double)decomposition[1], (double)decomposition[2]});

  std::array<double, 3> resultingLocalBoxExtend =
      autopas::utils::ArrayMath::sub(domainDecomposition.getLocalBoxMax(), domainDecomposition.getLocalBoxMin());

  EXPECT_NEAR(expectedLocalBoxExtend[0], resultingLocalBoxExtend[0], 1e-10);
  EXPECT_NEAR(expectedLocalBoxExtend[1], resultingLocalBoxExtend[1], 1e-10);
  EXPECT_NEAR(expectedLocalBoxExtend[2], resultingLocalBoxExtend[2], 1e-10);
}

TEST_F(TestRegularGridDecomposition, testExchangeHaloParticles) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename", std::string(YAMLDIRECTORY) + "haloParticleTest.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);
  std::cout << "Configuration: " << configuration.cutoff.value << ", " << configuration.verletSkinRadius.value << std::endl;
  std::cout << "BoxMin: " << autopas::utils::ArrayUtils::to_string(configuration.boxMin.value) << std::endl;
  std::cout << "BoxMax: " << autopas::utils::ArrayUtils::to_string(configuration.boxMax.value) << std::endl;

  std::array<double, 3> localBoxMin = configuration.boxMin.value;
  std::array<double, 3> localBoxMax = configuration.boxMax.value;

  RegularGridDecomposition domainDecomposition(configuration.boxMin.value, configuration.boxMax.value, configuration.cutoff.value, configuration.verletSkinRadius.value);

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

  initializeAutoPasContainer(autoPasContainer, configuration);

  std::vector<std::vector<double>> particlePositions = {
    {1.5, 1.5, 1.5}, {5.0, 1.5, 1.5}, {8.5, 1.5,  1.5},
    {1.5, 5.0, 1.5}, {5.0, 5.0, 1.5}, {8.5, 5.0,  1.5},
    {1.5, 8.5, 1.5}, {5.0, 8.5, 1.5}, {8.5, 8.5,  1.5},

    {1.5, 1.5, 5.0}, {5.0, 1.5, 5.0}, {8.5, 1.5,  5.0},
    {1.5, 5.0, 5.0}, {5.0, 5.0, 5.0}, {8.5, 5.0,  5.0},
    {1.5, 8.5, 5.0}, {5.0, 8.5, 5.0}, {8.5, 8.5,  5.0},

    {1.5, 1.5, 8.5}, {5.0, 1.5, 8.5}, {8.5, 1.5,  8.5},
    {1.5, 5.0, 8.5}, {5.0, 5.0, 8.5}, {8.5, 5.0,  8.5},
    {1.5, 8.5, 8.5}, {5.0, 8.5, 8.5}, {8.5, 8.5,  8.5}
  };
  
  size_t id = 0;
  for (const auto position : particlePositions) {
      ParticleType particle;
      particle.setID(id);
      particle.setR({position[0], position[1], position[2]});

      autoPasContainer->addParticle(particle);

      ++id;
  }

  EXPECT_NO_THROW(domainDecomposition.exchangeHaloParticles(autoPasContainer));

  for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::halo); particle.isValid(); ++particle) {
    std::cout << particle->getID() << ", " << autopas::utils::ArrayUtils::to_string(particle->getR()) << std::endl;
  }

  const size_t haloParticleCount = autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo);
  EXPECT_EQ(haloParticleCount,98);
}

TEST_F(TestRegularGridDecomposition, testExchangeMigratingParticles) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename", std::string(YAMLDIRECTORY) + "cubeGrid.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};
  MDFlexConfig configuration(3, argv);

  std::array<double, 3> globalBoxMin = {1.0, 1.0, 1.0};
  std::array<double, 3> globalBoxMax = {10.0, 10.0, 10.0};

  RegularGridDecomposition domainDecomposition(globalBoxMin, globalBoxMax, 0, 0);

  std::array<double, 3> localBoxMin = domainDecomposition.getLocalBoxMin();
  std::array<double, 3> localBoxMax = domainDecomposition.getLocalBoxMax();
  for (int i = 0; i < localBoxMin.size(); ++i) {
    configuration.boxMin.value[i] = localBoxMin[i];
    configuration.boxMax.value[i] = localBoxMax[i];
  }

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

  initializeAutoPasContainer(autoPasContainer, configuration);

  for (auto &particle : configuration.getParticles()) {
    if (domainDecomposition.isInsideLocalDomain(particle.getR())) {
      autoPasContainer->addParticle(particle);
    }
  }

  for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    std::array<double, 3> deltaPosition({0.01, 0.0, 0.0});
    particle->addR(deltaPosition);
  }

  EXPECT_NO_THROW(domainDecomposition.exchangeMigratingParticles(autoPasContainer));
}
