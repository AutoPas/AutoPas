/**
 * @file TestRegularGridDecomposition.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "TestRegularGridDecomposition.h"

#include "autopas/utils/WrapMPI.h"
#include "src/ParticleAttributes.h"
#include "src/ParticleSerializationTools.h"
#include "src/TypeDefinitions.h"
#include "src/configuration/MDFlexConfig.h"
#include "src/domainDecomposition/DomainTools.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"

namespace {
std::vector<double> sub(std::vector<double> a, std::vector<double> b) {
  std::vector<double> difference = a;
  if (a.size() == b.size()) {
    for (int i = 0; i < a.size(); ++i) {
      difference[i] = a[i] - b[i];
    }
  }
  return difference;
}

std::vector<double> div(std::vector<double> a, std::vector<int> b) {
  std::vector<double> result = a;
  if (a.size() == b.size()) {
    for (int i = 0; i < a.size(); ++i) {
      result[i] = a[i] / (double)b[i];
    }
  }
  return result;
}

void initializeAutoPasContainer(RegularGridDecomposition::SharedAutoPasContainer &autoPasContainer, MDFlexConfig &configuration) {
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
  // This using directive is necessary, because 'autopas::AUTOPAS_...' variables defined in WrapMPI.h do not exist
  // when compiling with MPI. When compiling without MPI the namespace prefix needs to be used.
  using namespace autopas;

  std::vector<double> globalBoxMin = {1.0, 1.0, 1.0};
  std::vector<double> globalBoxMax = {10.0, 10.0, 10.0};

  RegularGridDecomposition domainDecomposition(3, globalBoxMin, globalBoxMax, 0, 0);

  std::vector<double> globalBoxExtend = sub(globalBoxMax, globalBoxMin);

  int numberOfProcesses;
  AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numberOfProcesses);

  std::vector<int> decomposition;
  DomainTools::generateDecomposition(numberOfProcesses, 3, decomposition);

  std::vector<double> expectedLocalBoxExtend = div(globalBoxExtend, decomposition);

  std::vector<double> resultingLocalBoxExtend =
      sub(domainDecomposition.getLocalBoxMax(), domainDecomposition.getLocalBoxMin());

  EXPECT_NEAR(expectedLocalBoxExtend[0], resultingLocalBoxExtend[0], 1e-10);
  EXPECT_NEAR(expectedLocalBoxExtend[1], resultingLocalBoxExtend[1], 1e-10);
  EXPECT_NEAR(expectedLocalBoxExtend[2], resultingLocalBoxExtend[2], 1e-10);
}

TEST_F(TestRegularGridDecomposition, testExchangeHaloParticles) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename", std::string(YAMLDIRECTORY) + "cubeGrid.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);

  std::vector<double> globalBoxMin = {1.0, 1.0, 1.0};
  std::vector<double> globalBoxMax = {10.0, 10.0, 10.0};

  RegularGridDecomposition domainDecomposition(3, globalBoxMin, globalBoxMax, 0, 0);

  std::vector<double> localBoxMin = domainDecomposition.getLocalBoxMin();
  std::vector<double> localBoxMax = domainDecomposition.getLocalBoxMax();
  for (int i = 0; i < localBoxMin.size(); ++i) {
    configuration.boxMin.value[i] = localBoxMin[i];
    configuration.boxMax.value[i] = localBoxMax[i];
  }

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

  initializeAutoPasContainer(autoPasContainer, configuration);

  for (auto &particle : configuration.getParticles()) {
    if (domainDecomposition.isInsideLocalDomain(particle.position)) {
      autoPasContainer->addParticle(ParticleSerializationTools::convertParticleAttributesToParticle(particle));
    }
  }

  EXPECT_NO_THROW(domainDecomposition.exchangeHaloParticles(autoPasContainer));
}

TEST_F(TestRegularGridDecomposition, testExchangeMigratingParticles) {
  std::vector<std::string> arguments = {"md-flexible", "--yaml-filename", std::string(YAMLDIRECTORY) + "cubeGrid.yaml"};

  char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

  MDFlexConfig configuration(3, argv);

  std::vector<double> globalBoxMin = {1.0, 1.0, 1.0};
  std::vector<double> globalBoxMax = {10.0, 10.0, 10.0};

  RegularGridDecomposition domainDecomposition(3, globalBoxMin, globalBoxMax, 0, 0);

  std::vector<double> localBoxMin = domainDecomposition.getLocalBoxMin();
  std::vector<double> localBoxMax = domainDecomposition.getLocalBoxMax();
  for (int i = 0; i < localBoxMin.size(); ++i) {
    configuration.boxMin.value[i] = localBoxMin[i];
    configuration.boxMax.value[i] = localBoxMax[i];
  }

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

  initializeAutoPasContainer(autoPasContainer, configuration);

  for (auto &particle : configuration.getParticles()) {
    if (domainDecomposition.isInsideLocalDomain(particle.position)) {
      autoPasContainer->addParticle(ParticleSerializationTools::convertParticleAttributesToParticle(particle));
    }
  }

  for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    std::array<double, 3> deltaPosition({0.01, 0.0, 0.0});
    particle->addR(deltaPosition);
  }

  EXPECT_NO_THROW(domainDecomposition.exchangeMigratingParticles(autoPasContainer));
}
