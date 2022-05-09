/**
 * @file RegularGridDecompositionTest.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "RegularGridDecompositionTest.h"

#include <gmock/gmock-matchers.h>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapMPI.h"
#include "src/TypeDefinitions.h"
#include "src/configuration/MDFlexConfig.h"
#include "src/domainDecomposition/DomainTools.h"
#include "src/domainDecomposition/RegularGridDecomposition.h"

extern template class autopas::AutoPas<ParticleType>;

TEST_F(RegularGridDecompositionTest, testGetLocalDomain) {
  MDFlexConfig configuration(0, nullptr);
  configuration.boxMin.value = {1.0, 1.0, 1.0};
  configuration.boxMax.value = {10.0, 10.0, 10.0};
  configuration.subdivideDimension.value = {true, true, true};
  configuration.verletSkinRadius.value = 0;
  configuration.cutoff.value = 2.0;
  configuration.boundaryOption.value = {options::BoundaryTypeOption::periodic, options::BoundaryTypeOption::periodic,
                                        options::BoundaryTypeOption::periodic};

  RegularGridDecomposition domainDecomposition(configuration);

  const std::array<double, 3> globalBoxExtend =
      autopas::utils::ArrayMath::sub(configuration.boxMax.value, configuration.boxMin.value);

  const int numberOfProcesses = []() {
    int result;
    autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &result);
    return result;
  }();

  const std::array<int, 3> decomposition =
      DomainTools::generateDecomposition(numberOfProcesses, configuration.subdivideDimension.value);

  const std::array<double, 3> expectedLocalBoxExtend = autopas::utils::ArrayMath::div(
      globalBoxExtend, autopas::utils::ArrayUtils::static_cast_array<double>(decomposition));
  // make sure expectations make sense
  ASSERT_THAT(expectedLocalBoxExtend, ::testing::Each(::testing::Gt(1e-10)));

  const std::array<double, 3> resultingLocalBoxExtend =
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
  int numberOfRanks{};
  int myRank{};
  autopas::AutoPas_MPI_Comm_rank(AUTOPAS_MPI_COMM_WORLD, &myRank);
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numberOfRanks);

  std::cout << "My Rank: " << myRank << std::endl;

  MDFlexConfig configuration(0, nullptr);
  configuration.boxMin.value = {0., 0., 0.};
  configuration.cutoff.value = 2.5;
  configuration.verletSkinRadius.value = 0.5;
  const double interactionLength = configuration.cutoff.value + configuration.verletSkinRadius.value;
  // make sure evey rank gets exactly 3x3x3 cells
  const double localBoxLength = 3. * interactionLength;
  const double globalBoxLength = localBoxLength * numberOfRanks;
  configuration.boxMax.value = {globalBoxLength, globalBoxLength, globalBoxLength};
  configuration.boundaryOption.value = {options::BoundaryTypeOption::periodic, options::BoundaryTypeOption::periodic,
                                        options::BoundaryTypeOption::periodic};

  RegularGridDecomposition domainDecomposition(configuration);
  const auto &localBoxMin = domainDecomposition.getLocalBoxMin();
  const auto &localBoxMax = domainDecomposition.getLocalBoxMax();

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>();
  autoPasContainer->setBoxMin(localBoxMin);
  autoPasContainer->setBoxMax(localBoxMax);
  autoPasContainer->setCutoff(configuration.cutoff.value);
  autoPasContainer->setVerletSkin(configuration.verletSkinRadius.value);
  autoPasContainer->init();

  // Setup 27 particles of which 26 will be relevant during halo update. Imagine a rubik's cube where each cell
  // contains a single particle. This layout contains 8 particles with 3 adjacent cell which is outside the cube,
  // 12 particles with two adjacent cells which are outside the cube and 6 particles with a single adjacent cell
  // outside the cube. The 'first cell' is the leftmost in all dimensions.
  const auto particlePositions = [&]() {
    using autopas::utils::ArrayMath::mulScalar;
    using autopas::utils::ArrayMath::add;
    std::vector<std::array<double, 3>> positions{};
    positions.reserve(27);
    const auto midOfLocalBox = mulScalar(localBoxMax, 0.5);
    // halo particles should be placed cutoff/2 inside from the box border
    const auto midToParitcle1D = (localBoxLength - configuration.cutoff.value) / 2.;
    size_t id = 0;
    for (double z = -1; z <= 1; ++z) {
      for (double y = -1; y <= 1; ++y) {
        for (double x = -1; x <= 1; ++x) {
          const auto relativePosition = mulScalar(std::array<double, 3>{x, y, z}, midToParitcle1D);
          const auto absolutePosition = add(midOfLocalBox, relativePosition);
          positions.push_back(absolutePosition);

          ParticleType particle;
          particle.setID(id++);
          particle.setR(absolutePosition);
          autoPasContainer->addParticle(particle);
        }
      }
    }
    return positions;
  }();
  ASSERT_THAT(particlePositions, ::testing::SizeIs(27)) << "Test setup faulty!";
  ASSERT_EQ(autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned), 27)
      << "Not all setup particles added to the container!";

  const auto leavingParticles = autoPasContainer->updateContainer();
  ASSERT_EQ(leavingParticles.size(), 0) << "All particles should have been created inside the container!";

  // halos are generated here, so this is what we actually test
  domainDecomposition.exchangeHaloParticles(autoPasContainer);

  EXPECT_EQ(autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned), 27)
      << "Owned particles missing after halo exchange!";

  // expect particles to be all around the box since every particle creates multiple halo particles.
  const auto expectedHaloParticlePositions = [&]() {
    using autopas::utils::ArrayMath::mulScalar;
    using autopas::utils::ArrayMath::add;
    using autopas::utils::ArrayMath::addScalar;

    std::vector<std::array<double, 3>> positions{};
    positions.reserve(98);
    const auto midOfLocalBox = mulScalar(localBoxMax, 0.5);
    const auto midToParitcle1DNear = (localBoxLength - configuration.cutoff.value) / 2.;
    const auto midToParitcle1DFar = (localBoxLength + configuration.cutoff.value) / 2.;

    const std::vector<double> distances{-midToParitcle1DFar, -midToParitcle1DNear, 0., midToParitcle1DNear,
                                        midToParitcle1DFar};
    const auto relLocalBoxMin = mulScalar<double, 3>({1., 1., 1.}, -localBoxLength / 2.);
    const auto relLocalBoxMax = addScalar(relLocalBoxMin, localBoxLength);
    for (double z : distances) {
      for (double y : distances) {
        for (double x : distances) {
          const std::array<double, 3> relativePosition{x, y, z};
          // only add a particle if the position is in the halo region (=not in the inner box).
          // This is equivalent to at least one of x/y/z has to be ==abs(midToParticle1DFar)
          if (autopas::utils::notInBox(relativePosition, relLocalBoxMin, relLocalBoxMax)) {
            const auto absolutePosition = add(midOfLocalBox, relativePosition);
            positions.push_back(absolutePosition);
          }
        }
      }
    }
    return positions;
  }();
  ASSERT_EQ(expectedHaloParticlePositions.size(), 98) << "Expectation setup faulty!";

  EXPECT_EQ(autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo),
            expectedHaloParticlePositions.size());

  for (auto particleIter = autoPasContainer->begin(autopas::IteratorBehavior::halo); particleIter.isValid();
       ++particleIter) {
    EXPECT_THAT(expectedHaloParticlePositions, ::testing::Contains(particleIter->getR()));
  }
}

/**
 * This test is designed to check if particles are properly being migrated.
 * It uses a very specific set of particles create a controlled test case.
 * For more information see the comments in the test.
 */
TEST_F(RegularGridDecompositionTest, testExchangeMigratingParticles) {
  GTEST_SKIP() << "THIS TEST IS CURRENTLY BROKEN AND WILL BE FIXED IN https://github.com/AutoPas/AutoPas/pull/628/";

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

    RegularGridDecomposition domainDecomposition(configuration);

    auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>(std::cout);

    //    initializeAutoPasContainer(autoPasContainer, configuration);

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
