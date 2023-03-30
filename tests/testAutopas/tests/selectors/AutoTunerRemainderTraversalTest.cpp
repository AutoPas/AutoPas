/**
 * @file AutoTunerRemainderTraversalTest.cpp
 * @author F. Gratl
 * @date 28.11.2022
 */

#include "AutoTunerRemainderTraversalTest.h"

#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/selectors/AutoTuner.h"
#include "autopas/selectors/Configuration.h"
#include "autopas/selectors/tuningStrategy/FullSearch.h"
#include "testingHelpers/commonTypedefs.h"

void testRemainderTraversal(const std::vector<Molecule> &particles, const std::vector<Molecule> &haloParticles,
                            std::vector<std::vector<Molecule>> &particlesBuffer,
                            std::vector<std::vector<Molecule>> &haloParticlesBuffer) {
  constexpr double cutoff = 2.5;
  constexpr double cellSizeFactor = 1.;
  const std::set<autopas::Configuration> confSet(
      {{autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
        autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled, 5}});
  auto tuningStrategy = std::make_unique<autopas::FullSearch>(confSet);
  const std::array<double, 3> boxMin =  {0., 0., 0.};
  const std::array<double, 3> boxMax = {9., 9., 9.};
  autopas::AutoTuner<Molecule> autoTuner{
      boxMin, boxMax, cutoff, 0.05, 4, std::move(tuningStrategy), 0.3, 0.0, autopas::SelectorStrategyOption::fastestAbs,
      1000,   3};

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

  autopas::LJFunctor<Molecule> functor(cutoff);
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
  const auto &[choiceA, choiceB] = GetParam();

  std::vector<Molecule> containerParticles{};
  std::vector<Molecule> containerHaloParticles{};
  std::vector<std::vector<Molecule>> bufferParticles{static_cast<size_t>(autopas::autopas_get_max_threads())};
  std::vector<std::vector<Molecule>> bufferHaloParticles{static_cast<size_t>(autopas::autopas_get_max_threads())};

  auto addParticleToTest = [&](const auto &p, ParticleStorage storage, int subBuffer) {
    switch (storage) {
      case ParticleStorage::container:
        containerParticles.push_back(p);
        break;
      case ParticleStorage::containerHalo:
        containerHaloParticles.push_back(p);
        break;
      case ParticleStorage::buffer:
        bufferParticles[subBuffer].push_back(p);
        break;
      case ParticleStorage::bufferHalo:
        bufferHaloParticles[subBuffer].push_back(p);
        break;
    }
  };

  Molecule particle1{{1., 1., 1.}, {0., 0., 0.}, 0, 0};
  addParticleToTest(particle1, choiceA, 0);
  Molecule particle2{{2., 1., 1.}, {0., 0., 0.}, 1, 0};
  if (choiceB == ParticleStorage::bufferHalo or choiceB == ParticleStorage::containerHalo) {
    particle2.setR({-1., 1., 1.});
  }
  addParticleToTest(particle2, choiceB, 0);
  // if we check buffers make sure to also check that not only one sub buffer is used
  if (choiceB == ParticleStorage::buffer or choiceB == ParticleStorage::bufferHalo) {
    Molecule particle3{{1., 2., 1.}, {0., 0., 0.}, 1, 0};
    if (choiceB == ParticleStorage::bufferHalo) {
      particle3.setR({1., -1., 1.});
    }
    addParticleToTest(particle2, choiceB, autopas::autopas_get_max_threads() - 1);
  }

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