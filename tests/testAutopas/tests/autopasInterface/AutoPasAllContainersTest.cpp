/**
 * @file AutoPasAllContainersTest.cpp
 * @author F. Gratl
 * @date 29.03.23
 */

#include "AutoPasAllContainersTest.h"

#include "autopas/AutoPasImpl.h"
#include "autopas/utils/generators/UniformGenerator.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<Molecule>;

INSTANTIATE_TEST_SUITE_P(Generated, AutoPasAllContainersTest,
                         ::testing::ValuesIn(generateAllValidContainerConfigurations()));

/**
 * Fill a container with randomly generated particles, invoke addParticles, then check if the expected amount was added.
 */
TEST_P(AutoPasAllContainersTest, addParticlesTest) {
  /// Config
  constexpr auto numMols = 100;
  constexpr std::array<double, 3> boxMin{0., 0., 0.};
  constexpr std::array<double, 3> boxMax{10., 10., 10.};

  /// Setup
  // initialize AutoPas
  autopas::AutoPas<Molecule> autoPas;
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  const auto config = GetParam();
  autoPas.setAllowedContainers({config.container});
  autoPas.setCellSizeFactor(config.cellSizeFactor);
  // Arbitrarily allow all valid traversals, so AutoPas can be initialized.
  autoPas.setAllowedTraversals(autopas::TraversalOption::getAllOptions());
  autoPas.setNumSamples(autoPas.getTuningInterval());  // to avoid warnings
  autoPas.init();
  // generate some random molecules
  std::vector<Molecule> mols;
  mols.reserve(numMols);
  std::mt19937 generator(42);
  for (size_t i = 0; i < numMols; ++i) {
    mols.emplace_back(autopas::generators::UniformGenerator::randomPosition(generator, boxMin, boxMax),
                      std::array<double, 3>{0., 0., 0.}, i, 0);
  }

  /// Test
  // add every second particle
  autoPas.addParticlesIf(mols, [](const auto &p) { return p.getID() % 2 == 0; });

  /// Expectations
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::owned), numMols / 2)
      << "AutoPas' counters are wrong!";
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::halo), 0) << "AutoPas' counters are wrong!";
  // count what is actually in the container
  size_t numParticlesFound = 0;
  size_t numHalosFound = 0;
  for (const auto &p : autoPas) {
    if (p.isOwned()) {
      ++numParticlesFound;
    }
    if (p.isHalo()) {
      ++numHalosFound;
    }
  }
  EXPECT_EQ(numParticlesFound, numMols / 2) << "Either iterators are broken or the number of particles is wrong!";
  EXPECT_EQ(numHalosFound, 0) << "Either iterators are broken or the number of particles is wrong!";
}

/**
 * Checks that deleting a particle that is stored in the container only marks it as deleted until the container is
 * rebuilt. This test checks that:
 * a) No particle is moved, i.e. every particle is still stored at the same address.
 * b) Exactly the deleted particles are dummies.
 * c) updateContainer() removes the dummies.
 */
TEST_P(AutoPasAllContainersTest, deleteParticleOnlyMarksTest) {
  /// Config
  constexpr unsigned long numParticles = 4;
  constexpr std::array<double, 3> boxMin{0., 0., 0.};
  constexpr std::array<double, 3> boxMax{10., 10., 10.};
  constexpr double cutoff = 1.;
  constexpr double skin = 0.2;

  /// Setup
  // initialize AutoPas
  autopas::AutoPas<Molecule> autoPas;
  autoPas.setBoxMin(boxMin);
  autoPas.setBoxMax(boxMax);
  autoPas.setCutoff(cutoff);
  autoPas.setVerletSkin(skin);
  // Generate a valid arbitrary full configuration for the container configuration under test.
  const auto config = GetParam().generateFullConfig(autopas::InteractionTypeOption::pairwise).value();
  autoPas.setAllowedContainers({config.container});
  autoPas.setAllowedTraversals({config.traversal});
  autoPas.setAllowedLoadEstimators({config.loadEstimator});
  autoPas.setAllowedDataLayouts({config.dataLayout});
  autoPas.setAllowedNewton3Options({config.newton3});
  autoPas.setAllowedCellSizeFactors(autopas::NumberSetFinite<double>(std::set<double>({config.cellSizeFactor})));
  autoPas.setAllowedVecPatterns({config.vecPattern});
  autoPas.setNumSamples(autoPas.getTuningInterval());  // to avoid warnings
  autoPas.init();
  // We assume that the all particles are in the same cell (e.g. that the tested cell size factor is not too small)
  // If so, prefer to reduce the particle placement so they all fit instead of removing that container configuration
  // from the test.
  ASSERT_GT(config.cellSizeFactor * (cutoff + skin), 0.4);
  for (unsigned long id = 0; id < numParticles; ++id) {
    autoPas.addParticle(Molecule({0.1 + 0.1 * static_cast<double>(id), 0.1, 0.1}, {0., 0., 0.}, id, 0));
  }
  std::map<Molecule *, unsigned long> idsByAddress;
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    idsByAddress[&*iter] = iter->getID();
  }
  ASSERT_EQ(idsByAddress.size(), numParticles);

  /// Test
  // Delete particle 0 via a reference to it.
  std::set<unsigned long> deletedIDs{0};
  for (const auto &[address, id] : idsByAddress) {
    if (id == 0) {
      // The return of deleteParticle is whether an Iterator/ptr pointing to that particle would remain valid.
      // This should always be false for particles in the container, as we always mark particles kept there as dummies
      // to delete them.
      EXPECT_FALSE(autoPas.deleteParticle(*address));
    }
  }
  // Delete particle 1 via an iterator. VerletClusterLists is left out, because it keeps new particles in a buffer until
  // its towers are built and deleting from this buffer via an iterator moves particles.
  if (config.container != autopas::ContainerOption::verletClusterLists) {
    deletedIDs.insert(1);
    for (auto iter = autoPas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      if (iter->getID() == 1) {
        autoPas.deleteParticle(iter);
      }
    }
  }

  /// Expectations
  size_t numParticlesFound = 0;
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
    ++numParticlesFound;
    const auto addressAndID = idsByAddress.find(&*iter);
    ASSERT_NE(addressAndID, idsByAddress.end()) << "Particle " << iter->getID() << " was moved.";
    EXPECT_EQ(iter->getID(), addressAndID->second);
    EXPECT_EQ(iter->isDummy(), deletedIDs.contains(iter->getID()));
  }
  EXPECT_EQ(numParticlesFound, numParticles);
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::owned), numParticles - deletedIDs.size());

  const auto leavingParticles = autoPas.updateContainer();
  EXPECT_TRUE(leavingParticles.empty());
  numParticlesFound = 0;
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
    ++numParticlesFound;
  }
  EXPECT_EQ(numParticlesFound, numParticles - deletedIDs.size());
}
