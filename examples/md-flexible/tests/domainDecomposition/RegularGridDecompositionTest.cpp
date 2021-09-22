/**
 * @file RegularGridDecompositionTest.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "RegularGridDecompositionTest.h"

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapMPI.h"
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

TEST_F(RegularGridDecompositionTest, testGetLocalDomain) {
  std::array<double, 3> globalBoxMin = {1.0, 1.0, 1.0};
  std::array<double, 3> globalBoxMax = {10.0, 10.0, 10.0};
  std::array<bool, 3> subdivideDimension = {true, true, true};

  RegularGridDecomposition domainDecomposition(globalBoxMin, globalBoxMax, subdivideDimension, 0, 0);

  std::array<double, 3> globalBoxExtend = autopas::utils::ArrayMath::sub(globalBoxMax, globalBoxMin);

  int numberOfProcesses;
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numberOfProcesses);

  std::array<int, 3> decomposition;
  DomainTools::generateDecomposition(numberOfProcesses, subdivideDimension, decomposition);

  std::array<double, 3> expectedLocalBoxExtend = autopas::utils::ArrayMath::div(
      globalBoxExtend, {(double)decomposition[0], (double)decomposition[1], (double)decomposition[2]});

  std::array<double, 3> resultingLocalBoxExtend =
      autopas::utils::ArrayMath::sub(domainDecomposition.getLocalBoxMax(), domainDecomposition.getLocalBoxMin());

  EXPECT_NEAR(expectedLocalBoxExtend[0], resultingLocalBoxExtend[0], 1e-10);
  EXPECT_NEAR(expectedLocalBoxExtend[1], resultingLocalBoxExtend[1], 1e-10);
  EXPECT_NEAR(expectedLocalBoxExtend[2], resultingLocalBoxExtend[2], 1e-10);
}

/**
 * This test is designed to check if halo particles are properly being created.
 * It uses a very specific set of particles create a controlled test case.
 * For more information see the comments in the test.
 */
TEST_F(RegularGridDecompositionTest, testExchangeHaloParticles) {
  int numberOfRanks;
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &numberOfRanks);

  std::cout << "NumberOfRanks: " << numberOfRanks << std::endl;

  if (numberOfRanks != 1) {
    EXPECT_EQ(true, true);
  } else {
    std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                          std::string(YAMLDIRECTORY) + "particleExchangeTest.yaml"};

    char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

    MDFlexConfig configuration(3, argv);

    std::array<double, 3> localBoxMin = configuration.boxMin.value;
    std::array<double, 3> localBoxMax = configuration.boxMax.value;

    RegularGridDecomposition domainDecomposition(configuration.boxMin.value, configuration.boxMax.value,
                                                 configuration.subdivideDimension.value, configuration.cutoff.value,
                                                 configuration.verletSkinRadius.value);

    auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

    initializeAutoPasContainer(autoPasContainer, configuration);

    // Setup 27 particles of which 26 will be relevant during halo update. Imagine a rubik's cube where each cell
    // contains a single particle. This layout contains 8 particles with 3 adjacent cell which is outside the cube,
    // 12 particles with two adjacent cells which are outside the cube and 6 particles with a single adjacent cell
    // outside the cube.
    std::vector<std::vector<double>> particlePositions = {
        {1.5, 1.5, 1.5}, {5.0, 1.5, 1.5}, {8.5, 1.5, 1.5}, {1.5, 5.0, 1.5}, {5.0, 5.0, 1.5}, {8.5, 5.0, 1.5},
        {1.5, 8.5, 1.5}, {5.0, 8.5, 1.5}, {8.5, 8.5, 1.5}, {1.5, 1.5, 5.0}, {5.0, 1.5, 5.0}, {8.5, 1.5, 5.0},
        {1.5, 5.0, 5.0}, {5.0, 5.0, 5.0}, {8.5, 5.0, 5.0}, {1.5, 8.5, 5.0}, {5.0, 8.5, 5.0}, {8.5, 8.5, 5.0},
        {1.5, 1.5, 8.5}, {5.0, 1.5, 8.5}, {8.5, 1.5, 8.5}, {1.5, 5.0, 8.5}, {5.0, 5.0, 8.5}, {8.5, 5.0, 8.5},
        {1.5, 8.5, 8.5}, {5.0, 8.5, 8.5}, {8.5, 8.5, 8.5}};

    size_t id = 0;
    for (const auto position : particlePositions) {
      ParticleType particle;
      particle.setID(id);
      particle.setR({position[0], position[1], position[2]});

      autoPasContainer->addParticle(particle);

      ++id;
    }

    EXPECT_NO_THROW(domainDecomposition.exchangeHaloParticles(autoPasContainer));

    // The resulting haloParticleCount has to be 98 because, the 8 corner particles have to be replicated 7 times each,
    // resulting in 56 halo particles for the corner cells. The 12 particles contained
    // in the edge cells of the rubik's cube need to be replicated 3 times each raising the total number of halo
    // particles to 90. The remaining 6 particles with a single adjacent cell only produce a single halo particle each.
    // Therefor the total amount of particles is 98.
    std::vector<std::vector<double>> expectedHaloParticlePositions = {
        {-2.725, -2.725, -2.725}, {1.5, -2.725, -2.725},    {5, -2.725, -2.725},   {8.5, -2.725, -2.725},
        {12.725, -2.725, -2.725}, {-2.725, 1.5, -2.725},    {1.5, 1.5, -2.725},    {5, 1.5, -2.725},
        {8.5, 1.5, -2.725},       {12.725, 1.5, -2.725},    {-2.725, 5, -2.725},   {1.5, 5, -2.725},
        {5, 5, -2.725},           {8.5, 5, -2.725},         {12.725, 5, -2.725},   {-2.725, 8.5, -2.725},
        {1.5, 8.5, -2.725},       {5, 8.5, -2.725},         {8.5, 8.5, -2.725},    {12.725, 8.5, -2.725},
        {-2.725, 12.725, -2.725}, {1.5, 12.725, -2.725},    {5, 12.725, -2.725},   {8.5, 12.725, -2.725},
        {12.725, 12.725, -2.725}, {-2.725, -2.725, 1.5},    {1.5, -2.725, 1.5},    {5, -2.725, 1.5},
        {8.5, -2.725, 1.5},       {12.725, -2.725, 1.5},    {-2.725, 1.5, 1.5},    {12.725, 1.5, 1.5},
        {-2.725, 5, 1.5},         {12.725, 5, 1.5},         {-2.725, 8.5, 1.5},    {12.725, 8.5, 1.5},
        {-2.725, 12.725, 1.5},    {1.5, 12.725, 1.5},       {5, 12.725, 1.5},      {8.5, 12.725, 1.5},
        {12.725, 12.725, 1.5},    {-2.725, -2.725, 5},      {1.5, -2.725, 5},      {5, -2.725, 5},
        {8.5, -2.725, 5},         {12.725, -2.725, 5},      {-2.725, 1.5, 5},      {12.725, 1.5, 5},
        {-2.725, 5, 5},           {12.725, 5, 5},           {-2.725, 8.5, 5},      {12.725, 8.5, 5},
        {-2.725, 12.725, 5},      {1.5, 12.725, 5},         {5, 12.725, 5},        {8.5, 12.725, 5},
        {12.725, 12.725, 5},      {-2.725, -2.725, 8.5},    {1.5, -2.725, 8.5},    {5, -2.725, 8.5},
        {8.5, -2.725, 8.5},       {12.725, -2.725, 8.5},    {-2.725, 1.5, 8.5},    {12.725, 1.5, 8.5},
        {-2.725, 5, 8.5},         {12.725, 5, 8.5},         {-2.725, 8.5, 8.5},    {12.725, 8.5, 8.5},
        {-2.725, 12.725, 8.5},    {1.5, 12.725, 8.5},       {5, 12.725, 8.5},      {8.5, 12.725, 8.5},
        {12.725, 12.725, 8.5},    {-2.725, -2.725, 12.725}, {1.5, -2.725, 12.725}, {5, -2.725, 12.725},
        {8.5, -2.725, 12.725},    {12.725, -2.725, 12.725}, {-2.725, 1.5, 12.725}, {1.5, 1.5, 12.725},
        {5, 1.5, 12.725},         {8.5, 1.5, 12.725},       {12.725, 1.5, 12.725}, {-2.725, 5, 12.725},
        {1.5, 5, 12.725},         {5, 5, 12.725},           {8.5, 5, 12.725},      {12.725, 5, 12.725},
        {-2.725, 8.5, 12.725},    {1.5, 8.5, 12.725},       {5, 8.5, 12.725},      {8.5, 8.5, 12.725},
        {12.725, 8.5, 12.725},    {-2.725, 12.725, 12.725}, {1.5, 12.725, 12.725}, {5, 12.725, 12.725},
        {8.5, 12.725, 12.725},    {12.725, 12.725, 12.725}};

    const size_t haloParticleCount = autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo);

    EXPECT_EQ(haloParticleCount, expectedHaloParticlePositions.size());

    size_t index = 0;
    for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::halo); particle.isValid(); ++particle) {
      const auto particlePosition = particle->getR();
      EXPECT_NEAR(particlePosition[0], expectedHaloParticlePositions[index][0], 1e-13);
      EXPECT_NEAR(particlePosition[1], expectedHaloParticlePositions[index][1], 1e-13);
      EXPECT_NEAR(particlePosition[2], expectedHaloParticlePositions[index][2], 1e-13);

      ++index;
    }
  }
}

/**
 * This test is designed to check if particles are properly being migrated.
 * It uses a very specific set of particles create a controlled test case.
 * For more information see the comments in the test.
 */
TEST_F(RegularGridDecompositionTest, testExchangeMigratingParticles) {
  int numberOfRanks;
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &numberOfRanks);

  if (numberOfRanks != 1) {
    EXPECT_EQ(true, true);
  } else {
    std::vector<std::string> arguments = {"md-flexible", "--yaml-filename",
                                          std::string(YAMLDIRECTORY) + "particleExchangeTest.yaml"};
    char *argv[3] = {arguments[0].data(), arguments[1].data(), arguments[2].data()};

    MDFlexConfig configuration(3, argv);

    std::array<double, 3> localBoxMin = configuration.boxMin.value;
    std::array<double, 3> localBoxMax = configuration.boxMax.value;

    RegularGridDecomposition domainDecomposition(configuration.boxMin.value, configuration.boxMax.value,
                                                 configuration.subdivideDimension.value, configuration.cutoff.value,
                                                 configuration.verletSkinRadius.value);

    auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

    initializeAutoPasContainer(autoPasContainer, configuration);

    // Setup 27 particles. Imagine a rubik's cube where each cell contains a single particle.
    std::vector<std::vector<double>> particlePositions = {
        {1.5, 1.5, 1.5}, {5.0, 1.5, 1.5}, {8.5, 1.5, 1.5}, {1.5, 5.0, 1.5}, {5.0, 5.0, 1.5}, {8.5, 5.0, 1.5},
        {1.5, 8.5, 1.5}, {5.0, 8.5, 1.5}, {8.5, 8.5, 1.5}, {1.5, 1.5, 5.0}, {5.0, 1.5, 5.0}, {8.5, 1.5, 5.0},
        {1.5, 5.0, 5.0}, {5.0, 5.0, 5.0}, {8.5, 5.0, 5.0}, {1.5, 8.5, 5.0}, {5.0, 8.5, 5.0}, {8.5, 8.5, 5.0},
        {1.5, 1.5, 8.5}, {5.0, 1.5, 8.5}, {8.5, 1.5, 8.5}, {1.5, 5.0, 8.5}, {5.0, 5.0, 8.5}, {8.5, 5.0, 8.5},
        {1.5, 8.5, 8.5}, {5.0, 8.5, 8.5}, {8.5, 8.5, 8.5}};

    size_t id = 0;
    for (const auto position : particlePositions) {
      ParticleType particle;
      particle.setID(id);
      particle.setR({position[0], position[1], position[2]});

      autoPasContainer->addParticle(particle);

      ++id;
    }

    // Move particles outside the simulation box to make them migrate.
    // Particles in corner cells (of the rubik's cube) will be moved diagonally in all dimensions.
    // Particles in edge cells will be sifted diagonally in two dimiensions.
    // Particles in surface cells which are neither a corner nor a edge will be moved along a single dimension.
    // Particles which are not in a surface cell will not be moved at all.
    for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
      auto position = particle->getR();
      if (position[0] < 2.5) {
        position[0] -= 3.0;
      } else if (position[0] > 7.5) {
        position[0] += 3.0;
      }
      if (position[1] < 2.5) {
        position[1] -= 3.0;
      } else if (position[1] > 7.5) {
        position[1] += 3.0;
      }
      if (position[2] < 2.5) {
        position[2] -= 3.0;
      } else if (position[2] > 7.5) {
        position[2] += 3.0;
      }
      particle->setR(position);
    }

    auto emigrants = autoPasContainer->updateContainer();
    EXPECT_NO_THROW(domainDecomposition.exchangeMigratingParticles(autoPasContainer, emigrants));

    std::vector<std::array<double, 3>> expectedPositionsAfterMigration = {
        {9.725, 9.725, 9.725}, {5, 9.725, 9.725}, {0.275, 9.725, 9.725}, {9.725, 5, 9.725},
        {5, 5, 9.725},         {0.275, 5, 9.725}, {9.725, 0.275, 9.725}, {5, 0.275, 9.725},
        {0.275, 0.275, 9.725}, {9.725, 9.725, 5}, {5, 9.725, 5},         {0.275, 9.725, 5},
        {9.725, 5, 5},         {5, 5, 5},         {0.275, 5, 5},         {9.725, 0.275, 5},
        {5, 0.275, 5},         {0.275, 0.275, 5}, {9.725, 9.725, 0.275}, {5, 9.725, 0.275},
        {0.275, 9.725, 0.275}, {9.725, 5, 0.275}, {5, 5, 0.275},         {0.275, 5, 0.275},
        {9.725, 0.275, 0.275}, {5, 0.275, 0.275}, {0.275, 0.275, 0.275}};

    for (auto particle = autoPasContainer->begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
      const auto id = particle->getID();
      const auto position = particle->getR();
      EXPECT_NEAR(position[0], expectedPositionsAfterMigration[id][0], 1e-13);
      EXPECT_NEAR(position[1], expectedPositionsAfterMigration[id][1], 1e-13);
      EXPECT_NEAR(position[2], expectedPositionsAfterMigration[id][2], 1e-13);
    }
  }
}
