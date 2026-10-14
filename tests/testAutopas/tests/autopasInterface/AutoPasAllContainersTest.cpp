/**
 * @file AutoPasAllContainersTest.cpp
 * @author F. Gratl
 * @date 29.03.23
 */

#include "AutoPasAllContainersTest.h"

#include "autopas/AutoPasImpl.h"
#include "autopas/utils/generators/UniformGenerator.h"
#include "testingHelpers/EmptyPairwiseFunctor.h"
#include "testingHelpers/commonTypedefs.h"

extern template class autopas::AutoPas<Molecule>;
extern template bool autopas::AutoPas<Molecule>::computeInteractions(EmptyPairwiseFunctor<Molecule> *);

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
 * Checks that deleting a particle only marks it as deleted. This test checks that:
 * a) No particle is moved, i.e. every particle is still stored at the same address.
 * b) Exactly the deleted particles are dummies.
 * c) updateContainer() removes the dummies.
 */
TEST_P(AutoPasAllContainersTest, deleteParticleOnlyMarksTest) {
  /// Config
  // Particles 0-3 are stored in the container, 4-6 in the buffer for owned particles, and 7-9 in the buffer for halo
  // particles. In each group, the first particle is deleted via a reference to it and the second via an iterator.
  constexpr unsigned long numContainerParticles = 4;
  constexpr unsigned long numOwnedBufferParticles = 3;
  constexpr unsigned long numHaloBufferParticles = 3;
  constexpr unsigned long numParticles = numContainerParticles + numOwnedBufferParticles + numHaloBufferParticles;
  const std::set<unsigned long> deletedIDsViaReference{0, 4, 7};
  const std::set<unsigned long> deletedIDsViaIterator{1, 5, 8};
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
  // Rebuild in every iteration, so the updateContainer() at the end of the test removes all dummies.
  autoPas.setVerletRebuildFrequency(1);
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
  // We assume that the all particles in the container are in the same cell (e.g. that the tested cell size factor is
  // not too small) If so, prefer to reduce the particle placement so they all fit instead of removing that container
  // configuration from the test.
  ASSERT_GT(config.cellSizeFactor * (cutoff + skin), 0.4);
  for (unsigned long id = 0; id < numContainerParticles; ++id) {
    autoPas.addParticle(Molecule({0.1 + 0.1 * static_cast<double>(id), 0.1, 0.1}, {0., 0., 0.}, id, 0));
  }
  // Build the neighbor lists.
  ASSERT_TRUE(autoPas.updateContainer().empty());
  EmptyPairwiseFunctor<Molecule> functor;
  autoPas.computeInteractions(&functor);
  // Now we add buffer particles
  for (unsigned long id = numContainerParticles; id < numContainerParticles + numOwnedBufferParticles; ++id) {
    autoPas.addParticle(Molecule({0.1 + 0.1 * static_cast<double>(id), 0.1, 0.1}, {0., 0., 0.}, id, 0));
  }
  for (unsigned long id = numContainerParticles + numOwnedBufferParticles; id < numParticles; ++id) {
    autoPas.addHaloParticle(Molecule({0.1 + 0.1 * static_cast<double>(id), 0.1, -0.1}, {0., 0., 0.}, id, 0));
  }
  size_t numParticlesActuallyInContainer = 0;
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOrHalo | autopas::IteratorBehavior::containerOnly);
       iter.isValid(); ++iter) {
    ++numParticlesActuallyInContainer;
  }
  ASSERT_EQ(numParticlesActuallyInContainer, numContainerParticles)
      << "Particles that should be in the buffers are not.";
  std::map<Molecule *, unsigned long> idsByAddress;
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
    idsByAddress[&*iter] = iter->getID();
  }

  /// Test
  for (const auto &[address, id] : idsByAddress) {
    if (deletedIDsViaReference.contains(id)) {
      autoPas.deleteParticle(*address);
    }
  }
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOrHalo); iter.isValid(); ++iter) {
    if (deletedIDsViaIterator.contains(iter->getID())) {
      autoPas.deleteParticle(iter);
    }
  }
  std::set<unsigned long> deletedIDs{deletedIDsViaReference};
  deletedIDs.insert(deletedIDsViaIterator.begin(), deletedIDsViaIterator.end());

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
  // Four of the seven owned and two of the three halo particles were deleted.
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::owned), 3);
  EXPECT_EQ(autoPas.getNumberOfParticles(autopas::IteratorBehavior::halo), 1);

  const auto leavingParticles = autoPas.updateContainer();
  EXPECT_TRUE(leavingParticles.empty());
  numParticlesFound = 0;
  for (auto iter = autoPas.begin(autopas::IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
    ++numParticlesFound;
  }
  // updateContainer() now has removed all halo particles.
  EXPECT_EQ(numParticlesFound, 3);
}
