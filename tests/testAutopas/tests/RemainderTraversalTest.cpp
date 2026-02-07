/**
 * @file RemainderTraversalTest.cpp
 * @author F. Gratl
 * @date 28.11.2022
 */

#include "RemainderTraversalTest.h"

#include "autopas/LogicHandler.h"
#include "autopas/options/TuningMetricOption.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/Configuration.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/NumThreadGuard.h"
#include "testingHelpers/commonTypedefs.h"

/**
 * Can test all individual steps of iterate pairwise. Expects exactly two particles separated by 1.
 *
 * @note Tests invoking this function should have AutoPas logger instantiated (e.g. by inheriting from AutoPasTestBase).
 *
 * @note Buffers need to have at least one (empty) cell. They must not be empty.
 *
 * @param particlesContainerOwned
 * @param particlesBuffer
 * @param particlesHaloBuffer
 * @param n3 Newton3 on or off
 */
void testIteratePairwiseSteps(std::vector<Molecule> &particlesContainerOwned,
                              std::vector<Molecule> &particlesContainerHalo,
                              std::vector<autopas::FullParticleCell<Molecule>> &particlesBuffers,
                              std::vector<autopas::FullParticleCell<Molecule>> &particlesHaloBuffers,
                              autopas::Newton3Option n3, autopas::DataLayoutOption dataLayout) {
  // sanity check that there are exactly two particles in the test
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
      2)
      << "This test expects exactly two particles!";

  constexpr double cellSizeFactor = 1.;
  constexpr unsigned int verletRebuildFrequency = 10;
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
      .cutoff = 2.5,
      .verletSkin = 0.5,
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const std::set<autopas::Configuration> searchSpace(
      {{autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
        autopas::LoadEstimatorOption::none, dataLayout, n3, autopas::InteractionTypeOption::pairwise}});
  auto tunerManager = std::make_shared<autopas::TunerManager>(autoTunerInfo);
  tunerManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);
  autopas::LogicHandler<Molecule> logicHandler(tunerManager, logicHandlerInfo, verletRebuildFrequency, "");

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

  ASSERT_EQ(container.size(), 2 - numParticlesInBuffers - numParticlesHaloBuffers)
      << "Not all particles were added to the container! ParticlesBuffers(" << numParticlesInBuffers << ") HaloBuffer("
      << numParticlesHaloBuffers << ")";

  // create a functor that calculates globals!
  LJFunctorType</*shift*/ false, /*mixing*/ false, autopas::FunctorN3Modes::Both, /*globals*/ true> functor(
      logicHandlerInfo.cutoff);
  // Choose sigma != distance so we get Upot != 0
  constexpr double sigma = 2.;
  constexpr double epsilon = 1.;
  functor.setParticleProperties(24 * epsilon, sigma * sigma);
  // do the actual test
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);
  constexpr double expectedDist = 1.;
  const double expectedAbsForce =
      std::abs((24 * epsilon) / (expectedDist * expectedDist) *
               (std::pow(sigma / expectedDist, 6) - 2 * std::pow(sigma / expectedDist, 12)));
  std::array<double, 3> totalObservedForce = {0., 0., 0.};
  using autopas::utils::ArrayUtils::operator<<;
  using autopas::utils::ArrayMath::add;
  for (auto iter = container.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto particleForceL2 = autopas::utils::ArrayMath::L2Norm(iter->getF());
    totalObservedForce = add(totalObservedForce, iter->getF());
    EXPECT_NEAR(particleForceL2, expectedAbsForce, 1e-12)
        << "Force for particle " << iter->getID() << " in the container is wrong!";
  }
  const auto &[logicHandlerParticleBuffers, logicHandlerHaloBuffers] = logicHandler.getParticleBuffers();
  for (auto &particlesBuffer : logicHandlerParticleBuffers) {
    for (const auto &p : particlesBuffer) {
      const auto particleForceL2 = autopas::utils::ArrayMath::L2Norm(p.getF());
      totalObservedForce = add(totalObservedForce, p.getF());
      EXPECT_NEAR(particleForceL2, expectedAbsForce, 1e-12)
          << "Force for particle " << p.getID() << " in the particle buffer is wrong!";
    }
  }
  if (numParticlesHaloBuffers == 0 and particlesContainerHalo.empty()) {
    for (size_t dim = 0; dim < totalObservedForce.size(); ++dim) {
      EXPECT_NEAR(totalObservedForce[dim], 0, 1e-12)
          << "p1.f[" << dim << "] + p2.f[" << dim << "] does not add up to zero!";
    }
  }
  // if halo particles are involved only expect half the Upot
  const double expectedPotentialEnergy =
      4 * epsilon * (std::pow(sigma / expectedDist, 12.) - std::pow(sigma / expectedDist, 6.)) *
      ((numParticlesHaloBuffers != 0 or not particlesContainerHalo.empty()) ? 0.5 : 1);
  EXPECT_NEAR(expectedPotentialEnergy, functor.getPotentialEnergy(), 1e-12);

  if (::testing::Test::HasFailure()) {
    std::cout << "Failures occurred for DataLayout: " << dataLayout << std::endl;
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_container_container_NoN3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{6., 1., 1.}, {0., 0., 0.}, 0, 0},
      Molecule{{7., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::disabled, dataLayout);
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_container_containerHalo_NoN3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{0.5, 1., 1.}, {0., 0., 0.}, 0, 0},
  };
  std::vector<Molecule> particlesContainerHalo{
      Molecule{{-0.5, 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::disabled, dataLayout);
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_particleBuffer_container_NoN3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::disabled, dataLayout);
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_particleBuffer_containerHalo_NoN3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{
      Molecule{{-0.5, 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{0.5, 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::disabled, dataLayout);
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_particleBufferA_particleBufferA_NoN3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  particlesBuffers[0].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::disabled, dataLayout);
  }
}

// This test is not possible without OpenMP because there only exists one buffer
#ifdef AUTOPAS_USE_OPENMP
TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_particleBufferA_particleBufferB_NoN3) {
  NumThreadGuard threadGuard(2);
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{2};
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  particlesBuffers[1].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::disabled, dataLayout);
  }
}
#endif

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_haloBuffer_container_NoN3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  particlesHaloBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::disabled, dataLayout);
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_haloBuffer_particleBuffer_NoN3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  particlesHaloBuffers[0].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::disabled, dataLayout);
  }
}

/// Newton 3 enabled
TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_container_container_N3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{6., 1., 1.}, {0., 0., 0.}, 0, 0},
      Molecule{{7., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::enabled, dataLayout);
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_container_containerHalo_N3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{0.5, 1., 1.}, {0., 0., 0.}, 0, 0},
  };
  std::vector<Molecule> particlesContainerHalo{
      Molecule{{-0.5, 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::enabled, dataLayout);
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_particleBuffer_container_N3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::enabled, dataLayout);
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_particleBuffer_containerHalo_N3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{
      Molecule{{-0.5, 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{0.5, 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::enabled, dataLayout);
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_particleBufferA_particleBufferA_N3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  particlesBuffers[0].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::enabled, dataLayout);
  }
}

#ifdef AUTOPAS_USE_OPENMP
TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_particleBufferA_particleBufferB_N3) {
  NumThreadGuard threadGuard(2);
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  particlesBuffers[1].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::enabled, dataLayout);
  }
}
#endif

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_haloBuffer_container_N3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  particlesHaloBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::enabled, dataLayout);
  }
}

TEST_F(RemainderTraversalTest, testRemainderTraversalDirectly_haloBuffer_particleBuffer_N3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers{numBuffers};
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers{numBuffers};
  particlesHaloBuffers[0].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                             autopas::Newton3Option::enabled, dataLayout);
  }
}

void testRemainderTraversal(const std::vector<Molecule> &particles, const std::vector<Molecule> &haloParticles,
                            std::vector<autopas::FullParticleCell<Molecule>> &particlesBuffer,
                            std::vector<autopas::FullParticleCell<Molecule>> &haloParticlesBuffer,
                            autopas::DataLayoutOption dataLayout) {
  /// Setup AutoTuner
  constexpr double cellSizeFactor = 1.;
  constexpr unsigned int verletRebuildFrequency = 10;
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{9., 9., 9.},
      .cutoff = 2.5,
      .verletSkin = 0.5,
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 1000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};

  const std::set<autopas::Configuration> searchSpace(
      {{autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
        autopas::LoadEstimatorOption::none, dataLayout, autopas::Newton3Option::enabled,
        autopas::InteractionTypeOption::pairwise}});
  auto tunerManager = std::make_shared<autopas::TunerManager>(autoTunerInfo);
  tunerManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);
  autopas::LogicHandler<Molecule> logicHandler(tunerManager, logicHandlerInfo, verletRebuildFrequency, "");

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

  LJFunctorType<> functor(logicHandlerInfo.cutoff);
  functor.setParticleProperties(24, 1);
  // do the actual test
  logicHandler.computeInteractionsPipeline<decltype(functor)>(&functor, autopas::InteractionTypeOption::pairwise);

  for (const auto &p : logicHandler.getContainer()) {
    EXPECT_THAT(p.getF(), testing::Not(::testing::ElementsAreArray({0., 0., 0.})))
        << "Particle in container had no interaction!\n"
        << p << "\nwith Data Layout " << dataLayout;
  }
  const auto &[logicHandlerParticlesBuffer, logicHandlerHaloBuffer] = logicHandler.getParticleBuffers();
  for (const auto &buffer : logicHandlerParticlesBuffer) {
    for (const auto &p : buffer) {
      EXPECT_THAT(p.getF(), testing::Not(::testing::ElementsAreArray({0., 0., 0.})))
          << "Particle in particlesBuffer had no interaction!\n"
          << p << "\nwith Data Layout " << dataLayout;
    }
  }
}

/**
 * Add a particle to one storage location and one to another (or the same) and check if they interact.
 */
TEST_P(RemainderTraversalTest, testRemainderTraversal) {
  /// SETUP
  const auto &[choiceA, choiceB] = GetParam();
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

  // create two particles. p1 in the container p2 outside if it is a halo.
  Molecule particle1{{1., 1., 1.}, {0., 0., 0.}, 0, 0};
  addParticleToTest(particle1, choiceA, 0);
  Molecule particle2{{2., 1., 1.}, {0., 0., 0.}, 1, 0};
  if (choiceB == ParticleStorage::bufferHalo or choiceB == ParticleStorage::containerHalo) {
    particle2.setR({-1., 1., 1.});
  }
  addParticleToTest(particle2, choiceB, 0);
  // if we test buffers, make sure to also test that not only one sub buffer is used
  if (choiceB == ParticleStorage::buffer or choiceB == ParticleStorage::bufferHalo) {
    Molecule particle3{{1., 2., 1.}, {0., 0., 0.}, 1, 0};
    if (choiceB == ParticleStorage::bufferHalo) {
      particle3.setR({1., -1., 1.});
    }
    addParticleToTest(particle3, choiceB, autopas::autopas_get_max_threads() - 1);
  }

  /// TEST
  for (const auto dataLayout : autopas::DataLayoutOption::getAllOptions()) {
    testRemainderTraversal(containerParticles, containerHaloParticles, bufferParticles, bufferHaloParticles,
                           dataLayout);
  }
}

INSTANTIATE_TEST_SUITE_P(Generated, RemainderTraversalTest,
                         ::testing::ValuesIn(std::vector<std::tuple<ParticleStorage, ParticleStorage>>{
                             {ParticleStorage::container, ParticleStorage::container},
                             {ParticleStorage::container, ParticleStorage::containerHalo},
                             {ParticleStorage::container, ParticleStorage::buffer},
                             {ParticleStorage::container, ParticleStorage::bufferHalo},
                             {ParticleStorage::buffer, ParticleStorage::containerHalo},
                             {ParticleStorage::buffer, ParticleStorage::buffer},
                             {ParticleStorage::buffer, ParticleStorage::bufferHalo},
                         }),
                         RemainderTraversalTest::twoParamToString());
