/**
 * @file AutoTunerRemainderTraversalTest.cpp
 * @author F. Gratl
 * @date 28.11.2022
 */

#include "AutoTunerRemainderTraversalTest.h"

#include "autopas/selectors/AutoTuner.h"
#include "autopas/selectors/Configuration.h"
#include "autopas/selectors/tuningStrategy/FullSearch.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
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
                              autopas::Newton3Option n3) {
  // sanity check that there are exactly two particles in the test
  const auto numParticlesInBuffers =
      std::transform_reduce(particlesBuffers.begin(), particlesBuffers.end(), 0, std::plus<>(),
                            [](const auto &cell) { return cell.numParticles(); });
  const auto numParticlesHaloBuffers =
      std::transform_reduce(particlesHaloBuffers.begin(), particlesHaloBuffers.end(), 0, std::plus<>(), [](auto &cell) {
        // guarantee that all halo particles are actually tagged as such
        for (auto &p : cell) {
          p.setOwnershipState(autopas::OwnershipState::halo);
        }
        return cell.numParticles();
      });
  ASSERT_EQ(
      particlesContainerOwned.size() + particlesContainerHalo.size() + numParticlesInBuffers + numParticlesHaloBuffers,
      2)
      << "This test expects exactly two particles!";

  constexpr double cutoff = 2.5;
  constexpr double cellSizeFactor = 1.;
  const std::set<autopas::Configuration> confSet(
      {{autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
        autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, n3}});
  auto tuningStrategy = std::make_unique<autopas::FullSearch>(confSet);
  const std::array<double, 3> boxMin = {0., 0., 0.};
  const std::array<double, 3> boxMax = {10., 10., 10.};
  autopas::AutoTuner<Molecule> autoTuner{
      boxMin, boxMax, cutoff, 0.05, 4, std::move(tuningStrategy), 0.3, 0.0, autopas::SelectorStrategyOption::fastestAbs,
      1000,   3,      10};

  auto container = autoTuner.getContainer();
  for (const auto &p : particlesContainerOwned) {
    container->addParticle(p);
  }
  for (const auto &p : particlesContainerHalo) {
    container->addHaloParticle(p);
  }

  ASSERT_EQ(container->getNumberOfParticles(), 2 - numParticlesInBuffers - numParticlesHaloBuffers)
      << "Not all particles were added to the container! ParticlesBuffers(" << numParticlesInBuffers << ") HaloBuffer("
      << numParticlesHaloBuffers << ")";

  // create a functor that calculates globals!
  mdLib::LJFunctor<Molecule, /*shift*/ false, /*mixing*/ false, autopas::FunctorN3Modes::Both, /*globals*/ true>
      functor(cutoff);
  // Choose sigma != distance so we get Upot != 0
  constexpr double sigma = 2.;
  constexpr double epsilon = 1.;
  functor.setParticleProperties(24 * epsilon, sigma * sigma);
  // do the actual test
  autoTuner.iteratePairwise(&functor, false, particlesBuffers, particlesHaloBuffers);
  constexpr double expectedDist = 1.;
  const double expectedAbsForce =
      std::abs((24 * epsilon) / (expectedDist * expectedDist) *
               (std::pow(sigma / expectedDist, 6) - 2 * std::pow(sigma / expectedDist, 12)));
  std::array<double, 3> totalObservedForce = {0., 0., 0.};
  using autopas::utils::ArrayUtils::operator<<;
  using autopas::utils::ArrayMath::add;
  for (auto iter = container->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto particleForceL2 = autopas::utils::ArrayMath::L2Norm(iter->getF());
    totalObservedForce = add(totalObservedForce, iter->getF());
    EXPECT_NEAR(particleForceL2, expectedAbsForce, 1e-12)
        << "Force for particle " << iter->getID() << " in the container is wrong!";
  }
  for (size_t i = 0; i < particlesBuffers.size(); ++i) {
    for (const auto &p : particlesBuffers[i]) {
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
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_container_container_NoN3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{6., 1., 1.}, {0., 0., 0.}, 0, 0},
      Molecule{{7., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::disabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_container_containerHalo_NoN3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{0.5, 1., 1.}, {0., 0., 0.}, 0, 0},
  };
  std::vector<Molecule> particlesContainerHalo{
      Molecule{{-0.5, 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::disabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_particleBuffer_container_NoN3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::disabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_particleBuffer_containerHalo_NoN3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{
      Molecule{{-0.5, 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{0.5, 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::disabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_particleBufferA_particleBufferA_NoN3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  particlesBuffers[0].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::disabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_particleBufferA_particleBufferB_NoN3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  particlesBuffers[1].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::disabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_haloBuffer_container_NoN3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  particlesHaloBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::disabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_haloBuffer_particleBuffer_NoN3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  particlesHaloBuffers[0].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::disabled);
}

/// Newton 3 enabled
TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_container_container_N3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{6., 1., 1.}, {0., 0., 0.}, 0, 0},
      Molecule{{7., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::enabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_container_containerHalo_N3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{0.5, 1., 1.}, {0., 0., 0.}, 0, 0},
  };
  std::vector<Molecule> particlesContainerHalo{
      Molecule{{-0.5, 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::enabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_particleBuffer_container_N3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::enabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_particleBuffer_containerHalo_N3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{
      Molecule{{-0.5, 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{0.5, 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::enabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_particleBufferA_particleBufferA_N3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  particlesBuffers[0].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::enabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_particleBufferA_particleBufferB_N3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  particlesBuffers[1].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::enabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_haloBuffer_container_N3) {
  std::vector<Molecule> particlesContainerOwned{
      Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0},
  };
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  particlesHaloBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::enabled);
}

TEST_F(AutoTunerRemainderTraversalTest, testRemainderTraversalDirectly_haloBuffer_particleBuffer_N3) {
  std::vector<Molecule> particlesContainerOwned{};
  std::vector<Molecule> particlesContainerHalo{};
  std::vector<autopas::FullParticleCell<Molecule>> particlesBuffers(2);
  particlesBuffers[0].addParticle(Molecule{{1., 1., 1.}, {0., 0., 0.}, 0, 0});
  std::vector<autopas::FullParticleCell<Molecule>> particlesHaloBuffers(2);
  particlesHaloBuffers[0].addParticle(Molecule{{2., 1., 1.}, {0., 0., 0.}, 1, 0});
  testIteratePairwiseSteps(particlesContainerOwned, particlesContainerHalo, particlesBuffers, particlesHaloBuffers,
                           autopas::Newton3Option::enabled);
}

void testRemainderTraversal(const std::vector<Molecule> &particles, const std::vector<Molecule> &haloParticles,
                            std::vector<autopas::FullParticleCell<Molecule>> &particlesBuffer,
                            std::vector<autopas::FullParticleCell<Molecule>> &haloParticlesBuffer) {
  /// Setup AutoTuner
  constexpr double cutoff = 2.5;
  constexpr double cellSizeFactor = 1.;
  const std::set<autopas::Configuration> confSet(
      {{autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
        autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled}});
  auto tuningStrategy = std::make_unique<autopas::FullSearch>(confSet);
  const std::array<double, 3> boxMin = {0., 0., 0.};
  const std::array<double, 3> boxMax = {9., 9., 9.};
  autopas::AutoTuner<Molecule> autoTuner{
      boxMin, boxMax, cutoff, 0.05, 4, std::move(tuningStrategy), 0.3, 0.0, autopas::SelectorStrategyOption::fastestAbs,
      1000,   3,      10};

  // fill the container with the given particles
  for (const auto &p : particles) {
    autoTuner.getContainer()->addParticle(p);
  }
  ASSERT_EQ(autoTuner.getContainer()->getNumberOfParticles(), particles.size())
      << "Container contains incorrect number of particles!";
  for (const auto &p : haloParticles) {
    autoTuner.getContainer()->addHaloParticle(p);
  }
  ASSERT_EQ(autoTuner.getContainer()->getNumberOfParticles(), particles.size() + haloParticles.size())
      << "Container contains incorrect number of halo particles!";

  mdLib::LJFunctor<Molecule> functor(cutoff);
  functor.setParticleProperties(24, 1);
  // do the actual test
  autoTuner.iteratePairwise(&functor, false, particlesBuffer, haloParticlesBuffer);

  for (const auto &p : *(autoTuner.getContainer())) {
    EXPECT_THAT(p.getF(), testing::Not(::testing::ElementsAreArray({0., 0., 0.})))
        << "Particle in container had no interaction!\n"
        << p;
  }
  for (const auto &buffer : particlesBuffer) {
    for (const auto &p : buffer) {
      EXPECT_THAT(p.getF(), testing::Not(::testing::ElementsAreArray({0., 0., 0.})))
          << "Particle in particlesBuffer had no interaction!\n"
          << p;
    }
  }
}

/**
 * Add a particle to one storage location and one to another (or the same) and check if they interact.
 */
TEST_P(AutoTunerRemainderTraversalTest, testRemainderTraversal) {
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
    addParticleToTest(particle2, choiceB, autopas::autopas_get_max_threads() - 1);
  }

  /// TEST
  testRemainderTraversal(containerParticles, containerHaloParticles, bufferParticles, bufferHaloParticles);
}

INSTANTIATE_TEST_SUITE_P(Generated, AutoTunerRemainderTraversalTest,
                         ::testing::ValuesIn(std::vector<std::tuple<ParticleStorage, ParticleStorage>>{
                             {ParticleStorage::container, ParticleStorage::container},
                             {ParticleStorage::container, ParticleStorage::containerHalo},
                             {ParticleStorage::container, ParticleStorage::buffer},
                             {ParticleStorage::container, ParticleStorage::bufferHalo},
                             {ParticleStorage::buffer, ParticleStorage::containerHalo},
                             {ParticleStorage::buffer, ParticleStorage::buffer},
                             {ParticleStorage::buffer, ParticleStorage::bufferHalo},
                         }),
                         AutoTunerRemainderTraversalTest::twoParamToString());