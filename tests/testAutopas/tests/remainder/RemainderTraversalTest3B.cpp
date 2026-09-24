/**
 * @file RemainderTraversalTest3B.cpp
 * @author muehlhaeusser
 * @date 23.09.2023
 */

#include "RemainderTraversalTest3B.h"

#include "autopas/LogicHandler.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/Configuration.h"
#include "molecularDynamicsLibrary/AxilrodTellerMutoFunctor.h"
#include "testingHelpers/NumThreadGuard.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * Can test all individual steps of computeInteractions. Expects exactly three particles that form an equilateral
 * triangle of side length sqrt(2).
 *
 * @note Tests invoking this function should have AutoPas logger instantiated (e.g. by inheriting from AutoPasTestBase).
 * @note Buffers need to have at least one (empty) cell. They must not be empty.
 *
 * @param particlesContainerOwned
 * @param particlesContainerHalo
 * @param particlesBuffers
 * @param particlesHaloBuffers
 * @param n3 Newton3 on or off
 * @param dataLayout
 * @param containerConfig
 */
void testIterateTriwiseSteps(std::vector<Molecule> &particlesContainerOwned,
                             std::vector<Molecule> &particlesContainerHalo,
                             std::vector<autopas::FullParticleCell<Molecule>> &particlesBuffers,
                             std::vector<autopas::FullParticleCell<Molecule>> &particlesHaloBuffers,
                             autopas::Newton3Option n3, autopas::DataLayoutOption dataLayout,
                             const ContainerConfiguration &containerConfig) {
  // sanity check that there are exactly three particles in the test
  const auto numParticlesInBuffers = std::transform_reduce(particlesBuffers.begin(), particlesBuffers.end(), 0,
                                                           std::plus<>(), [](const auto &cell) { return cell.size(); });
  const auto numParticlesHaloBuffers =
      std::transform_reduce(particlesHaloBuffers.begin(), particlesHaloBuffers.end(), 0, std::plus<>(), [](auto &cell) {
        // guarantee that all halo particles are actually tagged as such
        for (auto &p : cell) {
          p.setOwnershipState(autopas::OwnershipState::halo);
        }
        return cell.size();
      });
  ASSERT_EQ(
      particlesContainerOwned.size() + particlesContainerHalo.size() + numParticlesInBuffers + numParticlesHaloBuffers,
      3)
      << "This test expects exactly three particles!";

  constexpr unsigned int verletRebuildFrequency = 10;
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
      .cutoff = 2.5,
      .verletSkin = 0.5,
  };
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  // Get a search space consisting of all valid configurations with the appropriate container, csf, data layout, and n3.
  const auto searchSpace = generateAllValidConfigurations(
      autopas::InteractionTypeOption::triwise, {containerConfig.container}, autopas::TraversalOption::getAllOptions(),
      autopas::LoadEstimatorOption::getAllOptions(), {dataLayout}, {n3}, {containerConfig.cellSizeFactor},
      {autopas::OpenMPKindOption::omp_static}, {1}, autopas::VectorizationPatternOption::getAllOptions(),
      /*throwIfNone*/ false);
  // Sanity check
  if (searchSpace.empty()) {
    FAIL() << "No valid configuration exists for " << containerConfig << " with Data Layout " << dataLayout
           << " and Newton3 " << n3 << ". Only combinations with a valid configuration should be tested!";
  }
  auto tunerManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tunerManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::triwise);

  autopas::LogicHandler<Molecule> logicHandler(tunerManager, logicHandlerInfo, verletRebuildFrequency, "",
                                               autoTunerInfo.aosSortingThreshold, autoTunerInfo.soaSortingThreshold);

  // create a functor that calculates globals!
  mdLib::AxilrodTellerMutoFunctor<Molecule, /*mixing*/ false, autopas::FunctorN3Modes::Both, /*globals*/ true> functor(
      logicHandlerInfo.cutoff);
  constexpr double nu = 1.;
  functor.setParticleProperties(nu);

  // Not every valid configuration is applicable to this domain (e.g. lc_c04 requires csf >= 1). The tuner rejects the
  // inapplicable ones, but if none are applicable, something is wrong.
  auto configIter = searchSpace.begin();
  while (configIter != searchSpace.end()) {
    auto [traversalPtr, _] = logicHandler.isConfigurationApplicable(*configIter, functor);
    if (traversalPtr) {
      break;
    }
    ++configIter;
  }
  if (configIter == searchSpace.end()) {
    FAIL() << "No valid configuration for " << containerConfig << " with Data Layout " << dataLayout << " and Newton3 "
           << n3 << " is applicable to the domain. This test potentially needs fixing.";
  }

  // Add particles. Calling add(Halo)Particle on a fresh logicHandler should place the particles directly in the
  // container.
  auto &container = logicHandler.getContainer();
  for (const auto &p : particlesContainerOwned) {
    logicHandler.addParticle(p);
  }
  for (const auto &p : particlesContainerHalo) {
    logicHandler.addHaloParticle(p);
  }
  logicHandler.setParticleBuffers(particlesBuffers, particlesHaloBuffers);

  ASSERT_EQ(container.size(), 3 - numParticlesInBuffers - numParticlesHaloBuffers)
      << "Not all particles were added to the container! ParticlesBuffers(" << numParticlesInBuffers << ") HaloBuffer("
      << numParticlesHaloBuffers << ")";

  // do the actual test
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::triwise);
  // All configurations in the search space share the container, so selecting one should never replace the container.
  ASSERT_EQ(&logicHandler.getContainer(), &container)
      << "The container was replaced while selecting a configuration, invalidating the container reference!";

  // Expected displacements of particles. Signs can differ
  using namespace autopas::utils::ArrayMath::literals;
  constexpr std::array<double, 3> drij{-1., 1., 0.};
  constexpr std::array<double, 3> drki{1., 0., -1.};

  constexpr double distSquared = 2.;  // Distance between every particle pair should be sqrt(2)
  constexpr double cosNum = -1.;      // e.g. dot(drij, drki)
  constexpr double cosAll = cosNum * cosNum * cosNum;
  constexpr double distSix = distSquared * distSquared * distSquared;
  const double invdr5 = nu / (distSix * distSix * std::sqrt(distSix));

  auto expectedFi = drij * (cosNum * cosNum - distSquared * distSquared + 5.0 * cosAll / distSquared) +
                    drki * (-cosNum * cosNum + distSquared * distSquared - 5.0 * cosAll / distSquared);
  expectedFi *= 3.0 * invdr5;

  // Same absolute force for all three particles if they were placed in an equilateral triangle
  const auto expectedAbsForce = autopas::utils::ArrayMath::L2Norm(expectedFi);

  std::array<double, 3> totalObservedForce = {0., 0., 0.};
  using autopas::utils::ArrayUtils::operator<<;
  for (auto iter = container.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto particleForceL2 = autopas::utils::ArrayMath::L2Norm(iter->getF());
    totalObservedForce += iter->getF();
    EXPECT_NEAR(particleForceL2, expectedAbsForce, 1e-12)
        << "Force for particle " << iter->getID() << " in the container is wrong!";
  }
  const auto &[logicHandlerParticleBuffers, logicHandlerHaloBuffers] = logicHandler.getParticleBuffers();
  for (auto &particlesBuffer : logicHandlerParticleBuffers) {
    for (const auto &p : particlesBuffer) {
      const auto particleForceL2 = autopas::utils::ArrayMath::L2Norm(p.getF());
      totalObservedForce += p.getF();
      EXPECT_NEAR(particleForceL2, expectedAbsForce, 1e-12)
          << "Force for particle " << p.getID() << " in the particle buffer is wrong!";
    }
  }
  if (numParticlesHaloBuffers == 0 and particlesContainerHalo.empty()) {
    for (size_t dim = 0; dim < totalObservedForce.size(); ++dim) {
      EXPECT_NEAR(totalObservedForce[dim], 0, 1e-12)
          << "p1.f[" << dim << "] + p2.f[" << dim << "] + p3.f[" << dim << "] does not add up to zero!";
    }
  }
  // if halo particles are involved only expect one or two thirds of the potential energy
  const double energyFactor =
      (3.0 - static_cast<double>(numParticlesHaloBuffers + particlesContainerHalo.size())) / 3.0;
  const double expectedPotentialEnergy = energyFactor * (distSix - 3.0 * cosAll) * invdr5;

  EXPECT_NEAR(expectedPotentialEnergy, functor.getPotentialEnergy(), 1e-12);
}

TEST_P(RemainderTraversalDirectTest3B, testRemainderTraversalDirectly_container_container_3B) {
  const auto &[containerConfig, dataLayout, n3] = GetParam();
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{6., 5., 5.}, {0., 0., 0.}, 0, 0},
      Molecule{{5., 6., 5.}, {0., 0., 0.}, 1, 0},
      Molecule{{5., 5., 6.}, {0., 0., 0.}, 2, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  testIterateTriwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);
}

TEST_P(RemainderTraversalDirectTest3B, testRemainderTraversalDirectly_container_containerHalo_3B) {
  const auto &[containerConfig, dataLayout, n3] = GetParam();
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{1.5, 0.5, 0.5}, {0., 0., 0.}, 0, 0},
      Molecule{{0.5, 1.5, 0.5}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{
      Molecule{{0.5, 0.5, -0.5}, {0., 0., 0.}, 2, 0},
  };
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  testIterateTriwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);

  std::vector<Molecule> particlesContainerOwned2{
      Molecule{{1.5, 0.5, 0.5}, {0., 0., 0.}, 0, 0},
  };
  std::vector<Molecule> particlesContainerHalo2{
      Molecule{{0.5, -0.5, 0.5}, {0., 0., 0.}, 1, 0},
      Molecule{{0.5, 0.5, -0.5}, {0., 0., 0.}, 2, 0},
  };
  testIterateTriwiseSteps(particlesContainerOwned2, particlesContainerHalo2, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);
}

TEST_P(RemainderTraversalDirectTest3B, testRemainderTraversalDirectly_particleBuffer_container_3B) {
  const auto &[containerConfig, dataLayout, n3] = GetParam();
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{6., 5., 5.}, {0., 0., 0.}, 0, 0},
      Molecule{{5., 6., 5.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{5., 5., 6.}, {0., 0., 0.}, 2, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  testIterateTriwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);

  std::vector<Molecule> particlesContainerOwned2{
      Molecule{{6., 5., 5.}, {0., 0., 0.}, 0, 0},
  };
  particlesBuffers[0].addParticle(Molecule{{5., 6., 5.}, {0., 0., 0.}, 1, 0});
  testIterateTriwiseSteps(particlesContainerOwned2, particlesContainerHalo, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);
}

TEST_P(RemainderTraversalDirectTest3B, testRemainderTraversalDirectly_particleBuffer_containerHalo_3B) {
  const auto &[containerConfig, dataLayout, n3] = GetParam();
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{
      Molecule{{-0.5, 0.5, 0.5}, {0., 0., 0.}, 0, 0},
      Molecule{{0.5, -0.5, 0.5}, {0., 0., 0.}, 1, 0},
  };
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{0.5, 0.5, 1.5}, {0., 0., 0.}, 2, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  testIterateTriwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);

  std::vector<Molecule> particlesContainerHalo2{
      Molecule{{-0.5, 0.5, 0.5}, {0., 0., 0.}, 0, 0},
  };
  particlesBuffers[0].addParticle(Molecule{{0.5, 1.5, 0.5}, {0., 0., 0.}, 1, 0});
  testIterateTriwiseSteps(particlesContainerOwned, particlesContainerHalo2, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);
}

TEST_P(RemainderTraversalDirectTest3B, testRemainderTraversalDirectly_particleBufferA_particleBufferA_3B) {
  const auto &[containerConfig, dataLayout, n3] = GetParam();
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{6., 5., 5.}, {0., 0., 0.}, 0, 0});
  particlesBuffers[0].addParticle(Molecule{{5., 6., 5.}, {0., 0., 0.}, 1, 0});
  particlesBuffers[0].addParticle(Molecule{{5., 5., 6.}, {0., 0., 0.}, 2, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  testIterateTriwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);
}

// This test is not possible without OpenMP because there only exists one buffer
#ifdef AUTOPAS_USE_OPENMP
TEST_P(RemainderTraversalDirectTest3B, testRemainderTraversalDirectly_particleBufferA_particleBufferB_3B) {
  const auto &[containerConfig, dataLayout, n3] = GetParam();
  NumThreadGuard threadGuard(3);
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{3};
  particlesBuffers[0].addParticle(Molecule{{6., 5., 5.}, {0., 0., 0.}, 0, 0});
  particlesBuffers[1].addParticle(Molecule{{5., 6., 5.}, {0., 0., 0.}, 1, 0});
  particlesBuffers[2].addParticle(Molecule{{5., 5., 6.}, {0., 0., 0.}, 2, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(3);
  testIterateTriwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);
}
#endif

TEST_P(RemainderTraversalDirectTest3B, testRemainderTraversalDirectly_haloBuffer_container_3B) {
  const auto &[containerConfig, dataLayout, n3] = GetParam();
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{1.5, 0.5, 0.5}, {0., 0., 0.}, 0, 0},
      Molecule{{0.5, 1.5, 0.5}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  particlesHaloBuffers[0].addParticle(Molecule{{0.5, 0.5, -0.5}, {0., 0., 0.}, 2, 0});
  testIterateTriwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);

  std::vector<Molecule> particlesContainerOwned2{
      Molecule{{1.5, 0.5, 0.5}, {0., 0., 0.}, 0, 0},
  };
  particlesHaloBuffers[0].addParticle(Molecule{{0.5, -0.5, 0.5}, {0., 0., 0.}, 1, 0});
  testIterateTriwiseSteps(particlesContainerOwned2, particlesContainerHalo, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);
}

TEST_P(RemainderTraversalDirectTest3B, testRemainderTraversalDirectly_haloBuffer_particleBuffer_3B) {
  const auto &[containerConfig, dataLayout, n3] = GetParam();
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{1.5, 0.5, 0.5}, {0., 0., 0.}, 0, 0});
  particlesBuffers[0].addParticle(Molecule{{0.5, 1.5, 0.5}, {0., 0., 0.}, 1, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  particlesHaloBuffers[0].addParticle(Molecule{{0.5, 0.5, -0.5}, {0., 0., 0.}, 2, 0});
  testIterateTriwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);

  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers2{numBuffers};
  particlesBuffers2[0].addParticle(Molecule{{1.5, 0.5, 0.5}, {0., 0., 0.}, 0, 0});
  particlesHaloBuffers[0].addParticle(Molecule{{0.5, -0.5, 0.5}, {0., 0., 0.}, 1, 0});
  testIterateTriwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers2, particlesHaloBuffers, n3,
                          dataLayout, containerConfig);
}

namespace {
/**
 * Generates all combinations of valid container configurations, data layouts, and Newton3 options for which a valid
 * triwise configuration exists. Combinations without one (e.g. containers without triwise traversals) are skipped.
 *
 * @note Only AoS is tested, as the AxilrodTellerMutoFunctor (currently the only triwise functor) has no SoA
 * implementation. Once a triwise functor supports SoA, this should be extended to all data layouts.
 * @return
 */
std::vector<RemainderTraversalDirectTest3B::ParamType> getDirectTestParams() {
  std::vector<RemainderTraversalDirectTest3B::ParamType> params;
  for (const auto &containerConfig : generateAllValidContainerConfigurations()) {
    for (const auto dataLayout : std::set<autopas::DataLayoutOption>{autopas::DataLayoutOption::aos}) {
      for (const auto n3 : autopas::Newton3Option::getAllOptions()) {
        if (getArbitraryConfiguration(autopas::InteractionTypeOption::triwise, containerConfig.container, std::nullopt,
                                      std::nullopt, dataLayout, n3, containerConfig.cellSizeFactor,
                                      /*ompKind*/ std::nullopt, /*ompChunkSize*/ std::nullopt,
                                      /*vectorPattern*/ std::nullopt, /*throwIfNone*/ false)
                .has_value()) {
          params.emplace_back(containerConfig, dataLayout, n3);
        }
      }
    }
  }
  return params;
}
}  // namespace

INSTANTIATE_TEST_SUITE_P(Generated, RemainderTraversalDirectTest3B, ::testing::ValuesIn(getDirectTestParams()),
                         RemainderTraversalDirectTest3B::paramToString());

void testRemainderTraversal3B(const std::vector<Molecule> &particles, const std::vector<Molecule> &haloParticles,
                              std::vector<autopas::FullParticleCell<Molecule>> &particlesBuffer,
                              std::vector<autopas::FullParticleCell<Molecule>> &haloParticlesBuffer,
                              autopas::Newton3Option n3, autopas::DataLayoutOption dataLayout,
                              const ContainerConfiguration &containerConfig) {
  /// Setup AutoTuner
  constexpr unsigned int verletRebuildFrequency = 10;
  constexpr autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
      .cutoff = 2.5,
      .verletSkin = 0.5,
  };
  constexpr autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  // Get a search space consisting of all valid configurations with the appropriate container, csf, data layout, and n3.
  // From this, we will take an arbitrary configuration, but as we need an applicable configuration, we get all
  // candidates and let the logicHandler find one.
  const auto searchSpace = generateAllValidConfigurations(
      autopas::InteractionTypeOption::triwise, {containerConfig.container}, autopas::TraversalOption::getAllOptions(),
      autopas::LoadEstimatorOption::getAllOptions(), {dataLayout}, {n3}, {containerConfig.cellSizeFactor});
  auto tunerManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  tunerManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::triwise);

  autopas::LogicHandler<Molecule> logicHandler(tunerManager, logicHandlerInfo, verletRebuildFrequency, "",
                                               autoTunerInfo.aosSortingThreshold, autoTunerInfo.soaSortingThreshold);

  mdLib::AxilrodTellerMutoFunctor<Molecule> functor(logicHandlerInfo.cutoff);
  functor.setParticleProperties(1.);

  // Not every valid configuration is applicable to this domain (e.g. lc_c04 requires csf >= 1). The tuner rejects the
  // inapplicable ones, but if none are applicable, something is wrong. Similarly to the LogicHandler's configuration
  // selection, loop until an applicable configuration is found.
  auto configIter = searchSpace.begin();
  while (configIter != searchSpace.end()) {
    auto [traversalPtr, _] = logicHandler.isConfigurationApplicable(*configIter, functor);
    if (traversalPtr) {
      break;
    }
    ++configIter;
  }
  if (configIter == searchSpace.end()) {
    FAIL() << "No valid configuration for " << containerConfig << " with Data Layout " << dataLayout << " and Newton3 "
           << n3 << " is applicable to the domain. This test potentially needs fixing.";
  }

  // fill the container with the given particles
  for (const auto &p : particles) {
    logicHandler.addParticle(p);
  }
  ASSERT_EQ(logicHandler.getContainer().size(), particles.size())
      << "Container contains incorrect number of particles!";
  for (const auto &p : haloParticles) {
    logicHandler.addHaloParticle(p);
  }
  ASSERT_EQ(logicHandler.getContainer().size(), particles.size() + haloParticles.size())
      << "Container contains incorrect number of halo particles!";

  logicHandler.setParticleBuffers(particlesBuffer, haloParticlesBuffer);

  // do the actual test
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::triwise);

  // AutoPas assumes forces only are considered by the owning container. If we wished to add functionality to handle
  // e.g. half-shell approaches, this test should be changed.
  for (auto iter = logicHandler.getContainer().begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    EXPECT_THAT(iter->getF(), testing::Not(::testing::ElementsAreArray({0., 0., 0.})))
        << "Particle in container had no interaction!\n"
        << *iter;
  }
  const auto &[logicHandlerParticlesBuffer, logicHandlerHaloBuffer] = logicHandler.getParticleBuffers();
  for (const auto &buffer : logicHandlerParticlesBuffer) {
    for (const auto &p : buffer) {
      EXPECT_THAT(p.getF(), testing::Not(::testing::ElementsAreArray({0., 0., 0.})))
          << "Particle in particlesBuffer had no interaction!\n"
          << p;
    }
  }
}

/**
 * Add particles to the specified storage locations and check if they interact.
 */
TEST_P(RemainderTraversalTest3B, testRemainderTraversal3B) {
  /// SETUP
  const auto &[choiceA, choiceB, choiceC, containerConfig, dataLayout, n3] = GetParam();
  // helper buffers to set up the test
  std::vector<Molecule> containerParticles{};
  std::vector<Molecule> containerHaloParticles{};
  std::vector<autopas::FullParticleCell<Molecule>> bufferParticles{
      static_cast<size_t>(autopas::autopas_get_max_threads())};
  std::vector<autopas::FullParticleCell<Molecule>> bufferHaloParticles{
      static_cast<size_t>(autopas::autopas_get_max_threads())};

  auto addParticleToTest = [&](const auto &p, ParticleStorage storage, int subBuffer) {
    switch (storage) {
      case ParticleStorage::container:
        containerParticles.push_back(p);
        break;
      case ParticleStorage::containerHalo:
        containerHaloParticles.push_back(p);
        break;
      case ParticleStorage::buffer:
        bufferParticles[subBuffer].addParticle(p);
        break;
      case ParticleStorage::bufferHalo:
        bufferHaloParticles[subBuffer].addParticle(p);
        break;
    }
  };

  // create three particles. p1 always in the container, p2 and p3 outside if they are a halo.
  Molecule particle1{{8.5, 9.5, 9.5}, {0., 0., 0.}, 0, 0};
  addParticleToTest(particle1, choiceA, 0);
  Molecule particle2{{9.5, 8.5, 9.5}, {0., 0., 0.}, 1, 0};
  Molecule particle3{{9.5, 9.5, 8.5}, {0., 0., 0.}, 2, 0};
  if (choiceB == ParticleStorage::bufferHalo or choiceB == ParticleStorage::containerHalo) {
    particle2.setR({9.5, 10.5, 9.5});
  }
  addParticleToTest(particle2, choiceB, 0);
  if (choiceC == ParticleStorage::bufferHalo or choiceC == ParticleStorage::containerHalo) {
    particle3.setR({9.5, 9.5, 10.5});
  }
  addParticleToTest(particle3, choiceC, 0);

  // if we test buffers, make sure to also test that not only one sub buffer is used
  if (choiceC == ParticleStorage::buffer or choiceC == ParticleStorage::bufferHalo) {
    Molecule particle4{{9.5, 9.5, 9.5}, {0., 0., 0.}, 3, 0};
    if (choiceB == ParticleStorage::bufferHalo) {
      particle4.setR({9.5, 10.5, 10.5});
    }
    addParticleToTest(particle4, choiceC, autopas::autopas_get_max_threads() - 1);
  }

  /// TEST
  testRemainderTraversal3B(containerParticles, containerHaloParticles, bufferParticles, bufferHaloParticles, n3,
                           dataLayout, containerConfig);
}

namespace {
/**
 * Generates the test parameters as the cross product of the particle storage combinations under test with all
 * combinations of container configurations, data layouts, and Newton3 options for which a valid configuration exists.
 * @return
 */
std::vector<RemainderTraversalTest3B::ParamType> getTestParams() {
  const std::vector<std::tuple<ParticleStorage, ParticleStorage, ParticleStorage>> storageCombinations{
      {ParticleStorage::container, ParticleStorage::container, ParticleStorage::container},
      {ParticleStorage::container, ParticleStorage::container, ParticleStorage::containerHalo},
      {ParticleStorage::container, ParticleStorage::containerHalo, ParticleStorage::containerHalo},
      {ParticleStorage::container, ParticleStorage::container, ParticleStorage::buffer},
      {ParticleStorage::container, ParticleStorage::buffer, ParticleStorage::buffer},
      {ParticleStorage::container, ParticleStorage::container, ParticleStorage::bufferHalo},
      {ParticleStorage::container, ParticleStorage::bufferHalo, ParticleStorage::bufferHalo},
      {ParticleStorage::buffer, ParticleStorage::buffer, ParticleStorage::containerHalo},
      {ParticleStorage::buffer, ParticleStorage::containerHalo, ParticleStorage::containerHalo},
      {ParticleStorage::buffer, ParticleStorage::buffer, ParticleStorage::buffer},
      {ParticleStorage::buffer, ParticleStorage::buffer, ParticleStorage::bufferHalo},
      {ParticleStorage::buffer, ParticleStorage::bufferHalo, ParticleStorage::bufferHalo},
  };
  std::vector<RemainderTraversalTest3B::ParamType> params;
  for (const auto &[containerConfig, dataLayout, n3] : getDirectTestParams()) {
    for (const auto &[choiceA, choiceB, choiceC] : storageCombinations) {
      params.emplace_back(choiceA, choiceB, choiceC, containerConfig, dataLayout, n3);
    }
  }
  return params;
}
}  // namespace

INSTANTIATE_TEST_SUITE_P(Generated, RemainderTraversalTest3B, ::testing::ValuesIn(getTestParams()),
                         RemainderTraversalTest3B::threeParamToString());
