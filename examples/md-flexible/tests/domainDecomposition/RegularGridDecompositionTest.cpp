/**
 * @file RegularGridDecompositionTest.cpp
 * @author J. Körner
 * @date 27/05/21
 */
#include "RegularGridDecompositionTest.h"

#include <gmock/gmock-matchers.h>

#include <algorithm>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapMPI.h"
#include "autopas/utils/WrapOpenMP.h"
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

/**
 * Verifies multi-threaded RegularGridDecomposition::exchangeHaloParticles with LinkedCells (multiple cells per thread),
 * checking both that parallel collection + corner-halo filtering matches a 1-thread reference and that shrinking the
 * active halo count on a subsequent step does not leak stale particles from the high-water-mark _haloParticles buffer.
 */
TEST_F(RegularGridDecompositionTest, testExchangeHaloParticlesParallelAndHighWaterMark) {
  int numRanks = 1;
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);
  if (numRanks != 1) {
    GTEST_SKIP() << "This test targets single-rank multi-threaded periodic halo exchange.";
  }

  const int savedMaxThreads = autopas::autopas_get_max_threads();
  struct ThreadCountRestorer {
    int n;
    ~ThreadCountRestorer() { autopas::autopas_set_num_threads(n); }
  } restorer{savedMaxThreads};

  auto makeSetup = []() {
    MDFlexConfig config(0, nullptr);
    config.boxMin.value = {0., 0., 0.};
    config.boxMax.value = {10., 10., 10.};
    config.cutoff.value = 2.0;
    config.verletSkinRadius.value = 0.5;
    config.verletRebuildFrequency.value = 2;
    config.subdivideDimension.value = {true, true, true};
    config.boundaryOption.value = {options::BoundaryTypeOption::periodic, options::BoundaryTypeOption::periodic,
                                   options::BoundaryTypeOption::periodic};

    auto decomp = std::make_shared<RegularGridDecomposition>(config);
    auto container = std::make_shared<autopas::AutoPas<ParticleType>>();
    container->setBoxMin(decomp->getLocalBoxMin());
    container->setBoxMax(decomp->getLocalBoxMax());
    container->setCutoff(config.cutoff.value);
    container->setVerletSkin(config.verletSkinRadius.value);
    container->setVerletRebuildFrequency(config.verletRebuildFrequency.value);
    container->setAllowedContainers({autopas::ContainerOption::linkedCells});
    container->setAllowedTraversals({autopas::TraversalOption::lc_c08});
    container->init();

    // Place 1 particle at the center of every cell in the 4x4x4 owned cell grid (64 particles total).
    size_t id = 0;
    for (int z = 0; z < 4; ++z) {
      for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
          const std::array<double, 3> pos{1.25 + 2.5 * x, 1.25 + 2.5 * y, 1.25 + 2.5 * z};
#if MD_FLEXIBLE_MODE == MULTISITE
          ParticleType p(pos, {0., 0., 0.}, {0.7071067811865475, 0.7071067811865475, 0., 0.}, {0., 0., 0.}, id++);
#else
          ParticleType p(pos, {0., 0., 0.}, id++);
#endif
          container->addParticle(p);
        }
      }
    }
    return std::make_pair(container, decomp);
  };

  auto collectSortedHalos = [](autopas::AutoPas<ParticleType> &container) {
    std::vector<std::pair<size_t, std::array<double, 3>>> halos;
    for (auto iter = container.begin(autopas::IteratorBehavior::halo); iter.isValid(); ++iter) {
      halos.emplace_back(iter->getID(), iter->getR());
    }
    std::sort(halos.begin(), halos.end());
    return halos;
  };

  // Reference run with 1 thread.
  autopas::autopas_set_num_threads(1);
  auto [refContainer, refDecomp] = makeSetup();
  (void)refContainer->updateContainer();
  refDecomp->exchangeHaloParticles(*refContainer);
  const auto refHalosStep1 = collectSortedHalos(*refContainer);
  ASSERT_GT(refHalosStep1.size(), 100ul);

  // Multi-threaded run with 4 threads (triggers non-empty _threadHaloBuffers on multiple threads).
  autopas::autopas_set_num_threads(4);
  auto [parContainer, parDecomp] = makeSetup();
  (void)parContainer->updateContainer();
  parDecomp->exchangeHaloParticles(*parContainer);
  const auto parHalosStep1 = collectSortedHalos(*parContainer);

  EXPECT_EQ(parHalosStep1.size(), refHalosStep1.size());
  EXPECT_EQ(parHalosStep1, refHalosStep1);

  // Step 2: Delete half of the owned particles so haloCount shrinks below _haloParticles.size() (high-water mark).
  for (auto iter = refContainer->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    if (iter->getID() % 2 == 0) {
      refContainer->deleteParticle(iter);
    }
  }
  for (auto iter = parContainer->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    if (iter->getID() % 2 == 0) {
      parContainer->deleteParticle(iter);
    }
  }

  autopas::autopas_set_num_threads(1);
  (void)refContainer->updateContainer();
  refDecomp->exchangeHaloParticles(*refContainer);
  const auto refHalosStep2 = collectSortedHalos(*refContainer);
  ASSERT_LT(refHalosStep2.size(), refHalosStep1.size());

  autopas::autopas_set_num_threads(4);
  (void)parContainer->updateContainer();
  parDecomp->exchangeHaloParticles(*parContainer);
  const auto parHalosStep2 = collectSortedHalos(*parContainer);

  EXPECT_EQ(parContainer->getNumberOfParticles(autopas::IteratorBehavior::halo), refHalosStep2.size());
  EXPECT_EQ(parHalosStep2, refHalosStep2);
}

/**
 * Verifies that multi-threaded RegularGridDecomposition::reflectParticlesAtBoundaries produces identical
 * repulsive forces on particles distributed across multiple LinkedCells cells near all 6 reflective walls
 * compared to single-threaded execution.
 */
TEST_F(RegularGridDecompositionTest, testReflectParticlesAtBoundariesParallel) {
  int numRanks = 1;
  autopas::AutoPas_MPI_Comm_size(AUTOPAS_MPI_COMM_WORLD, &numRanks);
  if (numRanks != 1) {
    GTEST_SKIP() << "This test targets single-rank multi-threaded reflective boundaries.";
  }

  const int savedMaxThreads = autopas::autopas_get_max_threads();
  struct ThreadCountRestorer {
    int n;
    ~ThreadCountRestorer() { autopas::autopas_set_num_threads(n); }
  } restorer{savedMaxThreads};

  auto makeSetup = []() {
    MDFlexConfig config(0, nullptr);
    config.epsilonMap.value.clear();
    config.sigmaMap.value.clear();
    config.massMap.value.clear();
    config.boxMin.value = {0., 0., 0.};
    config.boxMax.value = {5., 5., 5.};
    config.cutoff.value = 1.0;
    config.verletSkinRadius.value = 0.25;
    config.verletRebuildFrequency.value = 10;
    config.subdivideDimension.value = {true, true, true};
    config.boundaryOption.value = {options::BoundaryTypeOption::reflective, options::BoundaryTypeOption::reflective,
                                   options::BoundaryTypeOption::reflective};
    config.addSiteType(0, 1.0);
    config.addLJParametersToSite(0, 1.0, 1.0);

    auto decomp = std::make_shared<RegularGridDecomposition>(config);
    auto container = std::make_shared<autopas::AutoPas<ParticleType>>();
    container->setBoxMin(decomp->getLocalBoxMin());
    container->setBoxMax(decomp->getLocalBoxMax());
    container->setCutoff(config.cutoff.value);
    container->setVerletSkin(config.verletSkinRadius.value);
    container->setVerletRebuildFrequency(config.verletRebuildFrequency.value);
    container->setAllowedContainers({autopas::ContainerOption::linkedCells});
    container->setAllowedTraversals({autopas::TraversalOption::lc_c08});
    container->init();

    auto ppl = std::make_shared<ParticlePropertiesLibraryType>(config.cutoff.value);
    ppl->addSiteType(0, 1.0);
    ppl->addLJParametersToSite(0, 1.0, 1.0);
#if MD_FLEXIBLE_MODE == MULTISITE
    ppl->addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.});
#endif
    ppl->calculateMixingCoefficients();

    // Place particles within the reflective skin (distance 0.35 < 2^(1/6)*0.5) near all 6 faces across a 4x4 grid.
    size_t id = 0;
    const std::array<double, 4> gridCoords{0.6, 1.85, 3.15, 4.4};
    for (int dim = 0; dim < 3; ++dim) {
      for (double wallCoord : {0.35, 4.65}) {
        for (double u : gridCoords) {
          for (double v : gridCoords) {
            std::array<double, 3> pos{};
            pos[dim] = wallCoord;
            pos[(dim + 1) % 3] = u;
            pos[(dim + 2) % 3] = v;
#if MD_FLEXIBLE_MODE == MULTISITE
            ParticleType p(pos, {0., 0., 0.}, {1., 0., 0., 0.}, {0., 0., 0.}, id++);
#else
            ParticleType p(pos, {0., 0., 0.}, id++);
#endif
            container->addParticle(p);
          }
        }
      }
    }
    return std::make_tuple(container, decomp, ppl);
  };

  auto collectForces = [](autopas::AutoPas<ParticleType> &container) {
    std::vector<std::pair<size_t, std::array<double, 3>>> forces;
    for (auto iter = container.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      forces.emplace_back(iter->getID(), iter->getF());
    }
    std::sort(forces.begin(), forces.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
    return forces;
  };

  autopas::autopas_set_num_threads(1);
  auto [refContainer, refDecomp, refPpl] = makeSetup();
  refDecomp->reflectParticlesAtBoundaries(*refContainer, *refPpl);
  const auto refForces = collectForces(*refContainer);

  autopas::autopas_set_num_threads(4);
  auto [parContainer, parDecomp, parPpl] = makeSetup();
  parDecomp->reflectParticlesAtBoundaries(*parContainer, *parPpl);
  const auto parForces = collectForces(*parContainer);

  ASSERT_EQ(parForces.size(), refForces.size());
  ASSERT_EQ(parForces.size(), 96ul);
  for (size_t i = 0; i < refForces.size(); ++i) {
    EXPECT_EQ(parForces[i].first, refForces[i].first);
    const double forceNormSq = autopas::utils::ArrayMath::dot(refForces[i].second, refForces[i].second);
    EXPECT_GT(forceNormSq, 0.0);
    for (size_t d = 0; d < 3; ++d) {
      EXPECT_NEAR(parForces[i].second[d], refForces[i].second[d], 1e-12);
    }
  }
}
