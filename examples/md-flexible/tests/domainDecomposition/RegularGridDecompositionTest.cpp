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

/**
 * Generate a simulation setup depending on the number of MPI ranks available.
 * The domain will be a stacked tower along the x-axis.
 * @return tuple(autoPasContainer, domainDecomposition)
 */
auto initDomain(options::BoundaryTypeOption boundaryType) {
  const int numberOfProcesses = []() {
    int result;
    autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &result);
    return result;
  }();

  MDFlexConfig configuration(0, nullptr);
  configuration.boxMin.value = {0., 0., 0.};
  configuration.cutoff.value = 2.5;
  configuration.verletSkinRadius.value = 1.0;
  configuration.verletRebuildFrequency.value = 2;
  const double interactionLength = configuration.cutoff.value + configuration.verletSkinRadius.value;
  // make sure evey rank gets exactly 3x3x3 cells
  const double localBoxLength = 3. * interactionLength;
  const double globalBoxLength = localBoxLength * numberOfProcesses;
  configuration.subdivideDimension.value = {true, false, false};
  configuration.boxMax.value = {globalBoxLength, localBoxLength, localBoxLength};
  configuration.boundaryOption.value = {boundaryType, boundaryType, boundaryType};

  auto domainDecomposition = std::make_shared<RegularGridDecomposition>(configuration);
  const auto &localBoxMin = domainDecomposition->getLocalBoxMin();
  const auto &localBoxMax = domainDecomposition->getLocalBoxMax();

  auto autoPasContainer = std::make_shared<autopas::AutoPas<ParticleType>>();
  autoPasContainer->setBoxMin(localBoxMin);
  autoPasContainer->setBoxMax(localBoxMax);
  autoPasContainer->setCutoff(configuration.cutoff.value);
  autoPasContainer->setVerletSkin(configuration.verletSkinRadius.value);
  autoPasContainer->setVerletRebuildFrequency(configuration.verletRebuildFrequency.value);
  autoPasContainer->setAllowedContainers({autopas::ContainerOption::directSum});
  autoPasContainer->setAllowedTraversals({autopas::TraversalOption::ds_sequential});
  autoPasContainer->init();

  return std::make_tuple(autoPasContainer, domainDecomposition);
}

/**
 * Setup 27 particle positions INSIDE the given container.
 *
 * Particles are placed near each corner (8), the middle of each edge (12), and the center of each face (6) such that
 * they are cutoff/2 away from the nearby borders. In addition, a further particle is placed in the center of the
 * domain.
 * @param autoPasContainer
 * @return Vector of generated positions.
 */
auto generatePositionsInsideDomain(const autopas::AutoPas<ParticleType> &autoPasContainer) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto &localBoxMin = autoPasContainer.getBoxMin();
  const auto &localBoxMax = autoPasContainer.getBoxMax();
  const auto localBoxLength = localBoxMax - localBoxMin;

  std::vector<std::array<double, 3>> positions{};
  positions.reserve(27);
  const auto midOfLocalBox = localBoxMin + (localBoxLength * 0.5);
  // particles should be placed cutoff/2 inside from the box border
  const auto midToParticle = (localBoxLength - autoPasContainer.getCutoff()) * 0.5;
  for (double z = -1; z <= 1; ++z) {
    for (double y = -1; y <= 1; ++y) {
      for (double x = -1; x <= 1; ++x) {
        const auto relativePosition = std::array<double, 3>{x, y, z} * midToParticle;
        const auto absolutePosition = midOfLocalBox + relativePosition;
        positions.push_back(absolutePosition);
      }
    }
  }
  return positions;
}

/**
 * Setup 98 particle positions OUTSIDE the given container.
 *
 * Particles are placed around each corner (8*7), the middle of every edge (12*3), and face (6*1).
 * The positions are cutoff/2 away from each border.
 *
 * @param autoPasContainer
 * @return Vector of generated positions.
 */
auto generatePositionsOutsideDomain(const autopas::AutoPas<ParticleType> &autoPasContainer) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto &localBoxMin = autoPasContainer.getBoxMin();
  const auto &localBoxMax = autoPasContainer.getBoxMax();
  const auto localBoxLength = localBoxMax - localBoxMin;

  std::vector<std::array<double, 3>> positions{};
  positions.reserve(98);
  const auto midOfLocalBox = localBoxMin + (localBoxLength * 0.5);
  const auto midToParticleNear = (localBoxLength - autoPasContainer.getCutoff()) * 0.5;
  const auto midToParticleFar = (localBoxLength + autoPasContainer.getCutoff()) * 0.5;

  const std::array<std::vector<double>, 3> distances = {{
      {-midToParticleFar[0], -midToParticleNear[0], 0., midToParticleNear[0], midToParticleFar[0]},
      {-midToParticleFar[1], -midToParticleNear[1], 0., midToParticleNear[1], midToParticleFar[1]},
      {-midToParticleFar[2], -midToParticleNear[2], 0., midToParticleNear[2], midToParticleFar[2]},
  }};
  const auto relLocalBoxMin = localBoxLength * -.5;
  const auto relLocalBoxMax = relLocalBoxMin + localBoxLength;
  for (double z : distances[2]) {
    for (double y : distances[1]) {
      for (double x : distances[0]) {
        const std::array<double, 3> relativePosition{x, y, z};
        // only add a particle if the position is in the halo region (=not in the inner box).
        // This is equivalent to at least one of x/y/z has to be ==abs(midToParticle1DFar)
        if (autopas::utils::notInBox(relativePosition, relLocalBoxMin, relLocalBoxMax)) {
          const auto absolutePosition = midOfLocalBox + relativePosition;
          positions.push_back(absolutePosition);
        }
      }
    }
  }
  return positions;
}

/**
 * Check if RegularGridDecomposition correctly splits the global domain into equal sized subdomains.
 */
TEST_F(RegularGridDecompositionTest, testGetLocalDomain) {
  using namespace autopas::utils::ArrayMath::literals;
  auto [autoPasContainer, domainDecomposition] = initDomain(options::BoundaryTypeOption::periodic);

  const auto globalBoxExtend = domainDecomposition->getGlobalBoxMax() - domainDecomposition->getGlobalBoxMin();
  const auto decomposition = domainDecomposition->getDecomposition();

  const std::array<double, 3> expectedLocalBoxExtend =
      globalBoxExtend / autopas::utils::ArrayUtils::static_cast_copy_array<double>(decomposition);
  // make sure expectations make sense
  ASSERT_THAT(expectedLocalBoxExtend, ::testing::Each(::testing::Gt(1e-10)));

  const std::array<double, 3> resultingLocalBoxExtend =
      domainDecomposition->getLocalBoxMax() - domainDecomposition->getLocalBoxMin();

  EXPECT_NEAR(expectedLocalBoxExtend[0], resultingLocalBoxExtend[0], 1e-10);
  EXPECT_NEAR(expectedLocalBoxExtend[1], resultingLocalBoxExtend[1], 1e-10);
  EXPECT_NEAR(expectedLocalBoxExtend[2], resultingLocalBoxExtend[2], 1e-10);
}

/**
 * Check if RegularGridDecomposition::exchangeHaloParticles generates all expected halo particles.
 * Create a grid of particles at the border of each rank, then check if halo particles appear in the expected positions
 * and nowhere else.
 */
TEST_P(RegularGridDecompositionTest, testExchangeHaloParticles) {
  const auto boundaryType = GetParam();
  auto [autoPasContainer, domainDecomposition] = initDomain(boundaryType);
  const auto &localBoxMin = autoPasContainer->getBoxMin();
  const auto &localBoxMax = autoPasContainer->getBoxMax();

  const auto particlePositions = generatePositionsInsideDomain(*autoPasContainer);
  ASSERT_THAT(particlePositions, ::testing::SizeIs(27)) << "Test setup faulty!";
  {
    size_t id = 0;
    for (const auto &pos : particlePositions) {
#if MD_FLEXIBLE_MODE == MULTISITE
      ParticleType particle(pos, {0., 0., 0.}, {0.7071067811865475, 0.7071067811865475, 0., 0.}, {0., 0., 0.}, id++);
#else
      ParticleType particle(pos, {0., 0., 0.}, id++);
#endif
      autoPasContainer->addParticle(particle);
    }
  }
  ASSERT_EQ(autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned), 27)
      << "Not all setup particles added to the container!";

  const auto leavingParticles = autoPasContainer->updateContainer();
  ASSERT_EQ(leavingParticles.size(), 0) << "All particles should have been created inside the container!";

  // halos are generated here, so this is what we actually test
  domainDecomposition->exchangeHaloParticles(*autoPasContainer);

  EXPECT_EQ(autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned), 27)
      << "Owned particles missing after halo exchange!";

  if (boundaryType == options::BoundaryTypeOption::periodic) {
    // expect particles to be all around the box since every particle creates multiple halo particles.
    const auto expectedHaloParticlePositions = generatePositionsOutsideDomain(*autoPasContainer);
    ASSERT_EQ(expectedHaloParticlePositions.size(), 98) << "Expectation setup faulty!";

    EXPECT_EQ(autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo),
              expectedHaloParticlePositions.size());

    for (auto particleIter = autoPasContainer->begin(autopas::IteratorBehavior::halo); particleIter.isValid();
         ++particleIter) {
      EXPECT_THAT(expectedHaloParticlePositions, ::testing::Contains(particleIter->getR()));
    }
  } else {
    const auto globalMin = domainDecomposition->getGlobalBoxMin();
    const auto globalMax = domainDecomposition->getGlobalBoxMax();

    const int numberOfProcesses = []() {
      int result;
      autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &result);
      return result;
    }();

    // Expect halos only between domains, not at global boundaries, when they are not periodic.
    const auto haloNumber = autopas::utils::ArrayMath::isNearRel(localBoxMin, globalMin) or
                                    autopas::utils::ArrayMath::isNearRel(localBoxMax, globalMax)
                                ? (numberOfProcesses == 1 ? 0 : 9)
                                : 18;

    EXPECT_EQ(autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::halo), haloNumber);
  }
}

INSTANTIATE_TEST_SUITE_P(TestHaloParticles, RegularGridDecompositionTest,
                         testing::Values(options::BoundaryTypeOption::periodic, options::BoundaryTypeOption::reflective,
                                         options::BoundaryTypeOption::none));

/**
 * This test is designed to check if particles are properly being migrated.
 * It uses a very specific set of particles create a controlled test case.
 * For more information see the comments in the test.
 */
TEST_F(RegularGridDecompositionTest, testExchangeMigratingParticles) {
  auto [autoPasContainer, domainDecomposition] = initDomain(options::BoundaryTypeOption::periodic);

  const auto positionsOutsideSubdomain = generatePositionsOutsideDomain(*autoPasContainer);
  ASSERT_THAT(positionsOutsideSubdomain, ::testing::SizeIs(98));
  // generate 98 particles at positions inside the subdomain. Otherwise, we cannot add them.
  // Then move them to the positions slightly outside the domain.
  {
    size_t id = 0;
    for (const auto &_ : positionsOutsideSubdomain) {
#if MD_FLEXIBLE_MODE == MULTISITE
      ParticleType p(domainDecomposition->getLocalBoxMin(), {0., 0., 0.},
                     {0.7071067811865475, 0.7071067811865475, 0., 0.}, {0., 0., 0.}, id++);
#else
      ParticleType p(domainDecomposition->getLocalBoxMin(), {0., 0., 0.}, id++);
#endif
      autoPasContainer->addParticle(p);
    }
    autoPasContainer->forEach([&](auto &p) { p.setR(positionsOutsideSubdomain[p.getID()]); });
  }

  // test
  auto emigrants = autoPasContainer->updateContainer();
  ASSERT_THAT(emigrants, ::testing::SizeIs(positionsOutsideSubdomain.size()));
  domainDecomposition->exchangeMigratingParticles(*autoPasContainer, emigrants);

  // derive expectations
  const auto expectedPositions = generatePositionsInsideDomain(*autoPasContainer);
  // check expectations
  EXPECT_EQ(autoPasContainer->getNumberOfParticles(autopas::IteratorBehavior::owned), positionsOutsideSubdomain.size());
  autoPasContainer->forEach([&](const auto &p) { EXPECT_THAT(expectedPositions, ::testing::Contains(p.getR())); });
}
