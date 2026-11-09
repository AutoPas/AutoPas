/**
 * @file AllContainersTests.cpp
 * @author humig
 * @date 08.07.2019
 */

#include "AllContainersTests.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"
#include "autopas/utils/inBox.h"
#include "testingHelpers/PeriodicBoundariesHelpers.h"

INSTANTIATE_TEST_SUITE_P(Generated, AllContainersTests,
                         ::testing::Combine(::testing::ValuesIn(generateAllValidContainerConfigurations())),
                         AllContainersTests::getParamToStringFunction());

INSTANTIATE_TEST_SUITE_P(Generated, AllContainersTestsBothUpdates,
                         ::testing::Combine(::testing::ValuesIn(generateAllValidContainerConfigurations()),
                                            ::testing::Bool()),
                         AllContainersTestsBothUpdates::getParamToStringFunction());

/**
 * Checks if ParticleContainerInterface::getNumParticle() returns the correct number of particles.
 */
TEST_P(AllContainersTests, testGetNumberOfParticles) {
  auto container = getInitializedContainer<ParticleFP64>(std::get<0>(GetParam()));
  EXPECT_EQ(container->size(), 0);

  const std::array<double, 3> r = {2, 2, 2};
  const ParticleFP64 p(r, {0., 0., 0.}, 0);
  container->addParticle(p);
  EXPECT_EQ(container->size(), 1);

  const std::array<double, 3> r2 = {1.5, 2, 2};
  const ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  container->addParticle(p2);
  EXPECT_EQ(container->size(), 2);
}

/**
 * Checks if ParticleContainerInterface::deleteAllParticles() deletes all particles.
 */
TEST_P(AllContainersTests, testDeleteAllParticles) {
  auto container = this->getInitializedContainer<ParticleFP64>(std::get<0>(GetParam()));
  EXPECT_EQ(container->size(), 0);

  const std::array<double, 3> r = {2, 2, 2};
  const ParticleFP64 p(r, {0., 0., 0.}, 0);
  container->addParticle(p);

  const std::array<double, 3> r2 = {1.5, 2, 2};
  const ParticleFP64 p2(r2, {0., 0., 0.}, 1);
  container->addParticle(p2);
  EXPECT_EQ(container->size(), 2);

  container->deleteAllParticles();
  EXPECT_EQ(container->size(), 0);
}

/**
 * Checks if addParticle() only accepts particles in the domain and throws otherwise, and if addHaloParticle() never
 * throws.
 */
TEST_P(AllContainersTests, testParticleAdding) {
  auto container = getInitializedContainer<ParticleFP64>(std::get<0>(GetParam()));
  int id = 1;
  for (double x : {boxMin[0] - 1.5, boxMin[0] - .5, boxMin[0], boxMin[0] + 5., boxMax[0] - 0.001, boxMax[0],
                   boxMax[0] + .5, boxMax[0] + 1.5}) {
    for (double y : {boxMin[1] - 1.5, boxMin[1] - .5, boxMin[1], boxMin[1] + 5., boxMax[1] - 0.001, boxMax[1],
                     boxMax[1] + .5, boxMax[1] + 1.5}) {
      for (double z : {boxMin[2] - 1.5, boxMin[2] - .5, boxMin[2], boxMin[2] + 5., boxMax[2] - 0.001, boxMax[2],
                       boxMax[2] + .5, boxMax[2] + 1.5}) {
        autopas::ParticleBaseFP64 p({x, y, z}, {0., 0., 0.},
                                    id++);  // not const as updating ownership (as containers no longer handle this)
        if (x == -1.5 or y == -1.5 or z == -1.5 or x == 11.5 or y == 11.5 or z == 11.5) {
          EXPECT_ANY_THROW(container->addParticle(p));  // outside, therefore not ok!
          p.setOwnershipState(autopas::OwnershipState::halo);
          EXPECT_NO_THROW(container->addHaloParticle(p));  // much outside, still ok because it is ignored!
        } else if (x == 10. or y == 10. or z == 10. or x == -.5 or y == -.5 or z == -.5 or x == 10.5 or y == 10.5 or
                   z == 10.5) {
          EXPECT_ANY_THROW(container->addParticle(p));  // outside, therefore not ok!
          p.setOwnershipState(autopas::OwnershipState::halo);
          EXPECT_NO_THROW(container->addHaloParticle(p));  // outside, therefore ok!
        } else {
          EXPECT_NO_THROW(container->addParticle(p));  // inside, therefore ok!
          p.setOwnershipState(autopas::OwnershipState::halo);
          EXPECT_ANY_THROW(
              // inside, and not ok, as halo particles cannot be added inside of the domain!
              container->addHaloParticle(p));
        }
      }
    }
  }
}

/**
 * Add particles in the halo in every direction. Then delete them via deleteHaloParticles() and check that they are
 * gone.
 */
TEST_P(AllContainersTests, testDeleteHaloParticles) {
  using namespace autopas::utils::ArrayMath::literals;

  auto container = getInitializedContainer<ParticleFP64>(std::get<0>(GetParam()));

  const std::array<double, 3> zeros{0, 0, 0};

  // counter and id for particles to be added
  size_t numParticles = 0;

  // calculate some distances needed later
  auto domainSize = container->getBoxMax() - container->getBoxMin();
  auto domainSizeHalf = domainSize * 0.5;
  auto domainCenter = container->getBoxMin() + domainSizeHalf;
  auto interactionLengthHalf = container->getInteractionLength() * 0.5;
  auto distCenterToMidHalo =
      domainSizeHalf + std::array<double, 3>{interactionLengthHalf, interactionLengthHalf, interactionLengthHalf};

  // setup: add some halo particles on every side and corner
  for (const int x : {-1, 0, 1}) {
    for (const int y : {-1, 0, 1}) {
      for (const int z : {-1, 0, 1}) {
        if (x == 0 and y == 0 and z == 0) {
          continue;
        }
        auto pos = domainCenter + std::array<double, 3>{distCenterToMidHalo[0] * x, distCenterToMidHalo[1] * y,
                                                        distCenterToMidHalo[2] * z};
        const ParticleFP64 p{pos, zeros, numParticles++, autopas::OwnershipState::halo};
        container->addHaloParticle(p);
      }
    }
  }
  // sanity checks
  ASSERT_GT(numParticles, 0);
  ASSERT_EQ(container->size(), numParticles);

  // actual test:
  container->deleteHaloParticles();
  ASSERT_EQ(container->size(), 0);
}

/**
 * Checks if updateContainer() deletes particles in halo.
 */
TEST_P(AllContainersTestsBothUpdates, testUpdateContainerHalo) {
  auto container = getInitializedContainer<ParticleFP64>(std::get<0>(GetParam()));
  const autopas::ParticleBaseFP64 p({boxMin[0] - 0.5, boxMin[1] - 0.5, boxMin[2] - 0.5}, {0, 0, 0}, 42,
                                    autopas::OwnershipState::halo);
  container->addHaloParticle(p);

  EXPECT_EQ(container->size(), 1);
  EXPECT_EQ(container->begin()->getID(), 42);

  auto invalidParticles = container->updateContainer(std::get<1>(GetParam()));

  // no particle should be returned
  EXPECT_EQ(invalidParticles.size(), 0);

  // no particle should remain
  auto iter = container->begin();
  EXPECT_FALSE(iter.isValid());
}

class TestParticle : public ParticleFP64 {
 public:
  static std::atomic<unsigned long> numParticles;

  TestParticle(std::array<double, 3> r, std::array<double, 3> v, unsigned long id) : ParticleFP64(r, v, id) {
    ++numParticles;
  }
  TestParticle(const TestParticle &testParticle) : ParticleFP64(testParticle) { ++numParticles; }
  ~TestParticle() override { --numParticles; }
  using SoAArraysType = autopas::utils::SoAType<TestParticle *, unsigned long /*id*/, double /*x*/, double /*y*/,
                                                double /*z*/, double /*fx*/, double /*fy*/, double /*fz*/,
                                                autopas::OwnershipState /*ownershipState*/>::Type;
  template <AttributeNames attribute, std::enable_if_t<attribute == AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    return this;
  }
  template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    return ParticleFP64::get<attribute>();
  }
};
// Initialize the TestParticle counter
std::atomic<unsigned long> TestParticle::numParticles{0};

/**
 * Checks if updateContainer deletes dummy particles.
 * @param previouslyOwned Specifies whether the particle was previously owned.
 */
void AllContainersTestsBothUpdates::testUpdateContainerDeletesDummy(bool previouslyOwned) {
  // We need the container to use TestParticle!
  auto container = getInitializedContainer<TestParticle>(std::get<0>(GetParam()));
  bool keepNeighborListsValid = std::get<1>(GetParam());

  // Add particle
  {
    const double pos = previouslyOwned ? 0.5 : -0.5;
    TestParticle p({pos, pos, pos}, {0, 0, 0}, 42);
    if (previouslyOwned) {
      p.setOwnershipState(autopas::OwnershipState::owned);
      container->addParticle(p);
    } else {
      p.setOwnershipState(autopas::OwnershipState::halo);
      container->addHaloParticle(p);
    }
  }
  // Mark particle as deleted{
  auto iter = container->begin();
  ASSERT_TRUE(iter.isValid());
  autopas::internal::markParticleAsDeleted(*iter);

  // Check that we do not iterate over it.
  iter = container->begin();
  ASSERT_FALSE(iter.isValid());

  // This should remove the dummy particle(s), while not returning it as invalid particle.
  auto invalidParticles = container->updateContainer(keepNeighborListsValid);

  // No particle should be returned
  EXPECT_EQ(invalidParticles.size(), 0);
  // The particle should no longer exist if the neighbor lists are not valid.
  EXPECT_EQ(TestParticle::numParticles, keepNeighborListsValid ? 1 : 0);

  // no particle should remain, therefore, the iterator should be invalid!
  EXPECT_FALSE(container->begin().isValid());

  container->deleteAllParticles();
  ASSERT_EQ(TestParticle::numParticles, 0);
}

/**
 * Checks if updateContainer deletes dummy particles.
 */
TEST_P(AllContainersTestsBothUpdates, testUpdateContainerDeletesPreviouslyOwnedDummy) {
  testUpdateContainerDeletesDummy(true);
}

TEST_P(AllContainersTestsBothUpdates, testUpdateContainerDeletesPreviouslyHaloDummy) {
  testUpdateContainerDeletesDummy(false);
}

/**
 * This test checks the correct behavior of the updateContainer(true) call.
 * Hereby we check, that:
 * a) Halo and leaving particles are converted to dummy particles and NOT completely removed.
 * b) Owned particles are not moved within the container, i.e., the pointer to the particle is still the same.
 */
TEST_P(AllContainersTests, testUpdateContainerKeepsNeighborListsValidIfSpecified) {
  auto container = getInitializedContainer<ParticleFP64>(std::get<0>(GetParam()));

  {
    const ParticleFP64 p({-.1, -.1, -.1}, {0., 0., 0.}, 0, autopas::OwnershipState::halo);
    container->addHaloParticle(p);
  }

  {
    const ParticleFP64 p({.02, .1, .1}, {0., 0., 0.}, 1);
    container->addParticle(p);
  }

  {
    const ParticleFP64 p({1.23, .1, .1}, {0., 0., 0.}, 2);
    container->addParticle(p);
  }
  struct Values {
    unsigned long id;
    int occurrences;
  };
  std::map<ParticleFP64 *, Values> previous_particles;
  for (auto &&p : *container) {
    // Iterates over owned and halo particles!
    previous_particles[&p] = {p.getID(), 0};
    if (p.getID() == 1) {
      // moves particle 1 outside the container (leaving particle!)
      p.addR({-0.04, 0., 0.});
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
      p.resetRAtRebuild();
#endif
    } else if (p.getID() == 2) {
      // Should normally move it to next linked cell, here the cell should NOT change!
      p.addR({0.04, 0., 0.});
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
      p.resetRAtRebuild();
#endif
    }
  }

  auto leavingParticles = container->updateContainer(true);

  // Particle 1 should be returned!
  ASSERT_EQ(leavingParticles.size(), 1ul);
  EXPECT_EQ(leavingParticles[0].getID(), 1);

  for (auto iter = container->begin(autopas::IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
    ASSERT_EQ(previous_particles.count(&*iter), 1ul);

    ++previous_particles[&*iter].occurrences;
    if (iter->getID() == 0 or iter->getID() == 1) {
      // particle 0 and 1 should now be dummy!
      EXPECT_TRUE(iter->isDummy());
    } else {
      // particle 2 should still be owned!
      EXPECT_TRUE(iter->isOwned());
    }
  }

  // Every previous particle should be visited exactly once!
  for (auto [key, value] : previous_particles) {
    EXPECT_EQ(value.occurrences, 1) << "Particle with id " << value.id << " at address " << key
                                    << " not visited correctly!";
  }
}

/**
 * Checks that updateHaloParticle() updates the correct particle if several periodic images of the same particle are
 * stored for the case that the domain is so thin that these images would be stored in neighboring cells.
 * All images are added and the container is built, the owned particles are moved, some of them across the periodic
 * boundaries, and the images of the new positions are updated as in an iteration without a rebuild.
 * This test then checks that:
 * a) An image is updated always if the same image was stored less than half the skin away, or a dummy particle that
 * arose from an owned particle leaving the box and being deleted may be used as an update candidate by some containers.
 * b) Exactly the updated images are stored afterwards, i.e. no image is lost.
 */
TEST_P(AllContainersTests, testUpdateHaloParticlePeriodicImagesInThinDomain) {
  using namespace autopas::utils::ArrayMath::literals;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayUtils::to_string;

  // With this box and cutoff, the cell based containers have only one cell per dimension inside of the box, unless the
  // cell size factor is smaller than one. In the CSF < 1 case, this test will probably trivially pass, but we keep it
  // in for completeness.
  boxMax = {3., 3., 3.};
  cutoff = 2.5;
  auto container = getInitializedContainer<ParticleFP64>(std::get<0>(GetParam()));

  // Start position and displacement of the owned particles. The index is the particle ID.
  const std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> startPositionsAndDisplacements{
      // Stays inside of the box.
      {{1.5, 1.5, 1.5}, {0.05, 0.05, 0.05}},
      // Leaves the box through a face.
      {{1.5, 2.97, 1.5}, {0., 0.05, 0.}},
      // Leaves the box through an edge.
      {{2.97, 2.97, 1.5}, {0.05, 0.05, 0.}},
      // Leaves the box through a corner.
      {{2.97, 2.97, 2.97}, {0.05, 0.05, 0.05}},
      // Stays inside of the box, but moves further than half the skin. None of its images may be updated.
      {{0.5, 0.5, 0.5}, {0.15, 0., 0.}},
  };

  const auto generateHaloImages = [&](const std::array<double, 3> &position, size_t id) {
    const ParticleFP64 haloParticle(position, {0., 0., 0.}, id, autopas::OwnershipState::halo);
    return generatePeriodicImagesOfSingleParticle(haloParticle, boxMin, boxMax, cutoff + skin);
  };

  // Add the owned particles and their periodic images.
  for (size_t id = 0; id < startPositionsAndDisplacements.size(); ++id) {
    const auto &startPosition = startPositionsAndDisplacements[id].first;
    container->addParticle(ParticleFP64(startPosition, {0., 0., 0.}, id));
    for (const auto &image : generateHaloImages(startPosition, id)) {
      container->addHaloParticle(image);
    }
  }
  rebuildNeighborLists(*container, std::get<0>(GetParam()));

  // Move the owned particles.
  for (auto iter = container->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(startPositionsAndDisplacements[iter->getID()].second);
  }

  const auto leavingParticles = container->updateContainer(true);
  EXPECT_EQ(leavingParticles.size(), 3ul);

  // Update the periodic images of all particles.
  std::vector<std::array<double, 3>> updatedPositions;
  for (size_t id = 0; id < startPositionsAndDisplacements.size(); ++id) {
    const auto &[startPosition, displacement] = startPositionsAndDisplacements[id];
    const auto previousImages = generateHaloImages(startPosition, id);
    // Particles that left the box reenter it through the opposite periodic boundary.
    const auto newPosition = wrapIntoBox(startPosition + displacement, boxMin, boxMax);
    // A particle that left the box is still stored as a dummy at the position it moved to, and can in some containers
    // e.g. VCL be a valid update candidate.
    const auto leavingDummyPosition = startPosition + displacement;
    const bool leftBox = not autopas::utils::inBox(leavingDummyPosition, boxMin, boxMax);
    for (const auto &image : generateHaloImages(newPosition, id)) {
      const bool sameImageWasStored =
          std::any_of(previousImages.begin(), previousImages.end(), [&](const auto &previousImage) {
            const auto distance = previousImage.getR() - image.getR();
            return dot(distance, distance) < skin * skin / 4.;
          });
      const auto distanceToLeavingDummy = leavingDummyPosition - image.getR();
      const bool leavingDummyIsNearby =
          leftBox and dot(distanceToLeavingDummy, distanceToLeavingDummy) < skin * skin / 4.;
      const bool updated = container->updateHaloParticle(image);
      if (sameImageWasStored) {
        EXPECT_TRUE(updated) << "Image of particle " << id << " at " << to_string(image.getR());
      } else if (leavingDummyIsNearby) {
        // If no same image was stored but leavingDummyIsNearby, we allow either updated to be true or false, as some
        // containers may legitimately update an owned particle that left the box, became a dummy, and is therefore a
        // valid candidate for the halo update of itself.
      } else {
        EXPECT_FALSE(updated) << "Image of particle " << id << " at " << to_string(image.getR());
      }
      if (updated) {
        updatedPositions.push_back(image.getR());
      }
    }
  }

  std::vector<std::array<double, 3>> storedPositions;
  for (auto iter = container->begin(autopas::IteratorBehavior::halo); iter.isValid(); ++iter) {
    storedPositions.push_back(iter->getR());
  }
  EXPECT_THAT(storedPositions, ::testing::UnorderedElementsAreArray(updatedPositions));
}

/**
 * Particle with an attribute that the base particle does not have, like the custom attributes a functor might need.
 */
class ParticleWithFakeAttribute : public ParticleFP64 {
 public:
  ParticleWithFakeAttribute(const std::array<double, 3> &r, const std::array<double, 3> &v, unsigned long id,
                            double fakeAttribute, autopas::OwnershipState ownershipState)
      : ParticleFP64(r, v, id, ownershipState), _fakeAttribute(fakeAttribute) {}

  enum AttributeNames : int { ptr, id, posX, posY, posZ, forceX, forceY, forceZ, fakeAttribute, ownershipState };

  using SoAArraysType =
      autopas::utils::SoAType<ParticleWithFakeAttribute *, unsigned long /*id*/, double /*x*/, double /*y*/,
                              double /*z*/, double /*fx*/, double /*fy*/, double /*fz*/, double /*fakeAttribute*/,
                              autopas::OwnershipState /*ownershipState*/>::Type;

  template <AttributeNames attribute, std::enable_if_t<attribute == AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() {
    return this;
  }

  template <AttributeNames attribute, std::enable_if_t<attribute != AttributeNames::ptr, bool> = true>
  constexpr typename std::tuple_element<attribute, SoAArraysType>::type::value_type get() const {
    if constexpr (attribute == AttributeNames::id) {
      return getID();
    } else if constexpr (attribute == AttributeNames::posX) {
      return getR()[0];
    } else if constexpr (attribute == AttributeNames::posY) {
      return getR()[1];
    } else if constexpr (attribute == AttributeNames::posZ) {
      return getR()[2];
    } else if constexpr (attribute == AttributeNames::forceX) {
      return getF()[0];
    } else if constexpr (attribute == AttributeNames::forceY) {
      return getF()[1];
    } else if constexpr (attribute == AttributeNames::forceZ) {
      return getF()[2];
    } else if constexpr (attribute == AttributeNames::fakeAttribute) {
      return _fakeAttribute;
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      return this->_ownershipState;
    } else {
      autopas::utils::ExceptionHandler::exception("ParticleWithFakeAttribute::get() unknown attribute {}", attribute);
    }
  }

  template <AttributeNames attribute>
  constexpr void set(typename std::tuple_element<attribute, SoAArraysType>::type::value_type value) {
    if constexpr (attribute == AttributeNames::id) {
      setID(value);
    } else if constexpr (attribute == AttributeNames::posX) {
      this->_r[0] = value;
    } else if constexpr (attribute == AttributeNames::posY) {
      this->_r[1] = value;
    } else if constexpr (attribute == AttributeNames::posZ) {
      this->_r[2] = value;
    } else if constexpr (attribute == AttributeNames::forceX) {
      this->_f[0] = value;
    } else if constexpr (attribute == AttributeNames::forceY) {
      this->_f[1] = value;
    } else if constexpr (attribute == AttributeNames::forceZ) {
      this->_f[2] = value;
    } else if constexpr (attribute == AttributeNames::fakeAttribute) {
      _fakeAttribute = value;
    } else if constexpr (attribute == AttributeNames::ownershipState) {
      this->_ownershipState = value;
    } else {
      autopas::utils::ExceptionHandler::exception("ParticleWithFakeAttribute::set() unknown attribute {}", attribute);
    }
  }

  [[nodiscard]] double getFakeAttribute() const { return _fakeAttribute; }

 private:
  double _fakeAttribute;
};

/**
 * Checks that updateHaloParticle() overwrites the stored halo particle with all attributes of the given particle,
 * including those that only a derived particle class has, and that the stored particle is a halo particle again.
 *
 * Strictly speaking, AutoPas only needs to make sure that all potential attributes listed by getNeededAttribute e.g.
 * something like velocity is not needed. For now, this isn't possible so this test is overly strict.
 */
TEST_P(AllContainersTests, testUpdateHaloParticleUpdatesAllAttributes) {
  using namespace autopas::utils::ArrayMath::literals;

  const auto containerConfig = std::get<0>(GetParam());
  auto container = getInitializedContainer<ParticleWithFakeAttribute>(containerConfig);

  constexpr unsigned long id = 0;
  const std::array<double, 3> position{boxMin[0] - 0.1, 0.5, 0.5};
  container->addHaloParticle(ParticleWithFakeAttribute(position, {0., 0., 0.}, id, 1., autopas::OwnershipState::halo));
  rebuildNeighborLists(*container, containerConfig);

  // As in an iteration without a rebuild, this marks the halo particle as a dummy.
  ASSERT_TRUE(container->updateContainer(true).empty());
  ASSERT_FALSE(container->begin(autopas::IteratorBehavior::halo).isValid());

  // New values for all attributes except id and ownership state. The particle moved less than skin/2, as it would
  // between two rebuilds.
  ParticleWithFakeAttribute updatedParticle(position + std::array<double, 3>{0., 0.4 * skin, 0.}, {1., 2., 3.}, id, 2.,
                                            autopas::OwnershipState::halo);
  updatedParticle.setF({4., 5., 6.});
  EXPECT_TRUE(container->updateHaloParticle(updatedParticle));

  size_t numHalosFound = 0;
  for (auto iter = container->begin(autopas::IteratorBehavior::halo); iter.isValid(); ++iter) {
    ++numHalosFound;
    EXPECT_EQ(iter->getID(), id);
    EXPECT_EQ(iter->getR(), updatedParticle.getR());
    EXPECT_EQ(iter->getV(), updatedParticle.getV());
    EXPECT_EQ(iter->getF(), updatedParticle.getF());
    EXPECT_EQ(iter->getFakeAttribute(), updatedParticle.getFakeAttribute());
  }
  EXPECT_EQ(numHalosFound, 1);
}

/**
 * Checks that deleting a particle only marks it as deleted until the container is rebuilt. This test checks that:
 * a) No particle is moved, i.e. every particle is still stored at the same address.
 * b) Exactly the deleted particles are dummies.
 * c) updateContainer(false) removes the dummies.
 */
TEST_P(AllContainersTests, testDeleteParticleOnlyMarks) {
  const auto containerConfig = std::get<0>(GetParam());
  auto container = getInitializedContainer<ParticleFP64>(containerConfig);

  // We assume that the all particles are in the same cell (e.g. that the tested cell size factor is not too small)
  // If so, prefer to reduce the particle placement so they all fit instead of removing that container configuration
  // from the test.
  ASSERT_THAT(container->getTraversalSelectorInfo().cellLength, ::testing::Each(::testing::Gt(0.4)));
  const unsigned long numParticles = 4;
  for (unsigned long id = 0; id < numParticles; ++id) {
    container->addParticle(ParticleFP64({0.1 + 0.1 * static_cast<double>(id), 0.1, 0.1}, {0., 0., 0.}, id));
  }
  std::map<ParticleFP64 *, unsigned long> idsByAddress;
  for (auto iter = container->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    idsByAddress[&*iter] = iter->getID();
  }
  ASSERT_EQ(idsByAddress.size(), numParticles);

  // Delete particle 0 via a reference to it.
  std::set<unsigned long> deletedIDs{0};
  for (const auto &[address, id] : idsByAddress) {
    if (id == 0) {
      container->deleteParticle(*address);
    }
  }
  // Delete particle 1 via an iterator.
  deletedIDs.insert(1);
  for (auto iter = container->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    if (iter->getID() == 1) {
      autopas::internal::deleteParticle(iter);
    }
  }

  size_t numParticlesFound = 0;
  for (auto iter = container->begin(autopas::IteratorBehavior::ownedOrHaloOrDummy); iter.isValid(); ++iter) {
    ++numParticlesFound;
    const auto addressAndID = idsByAddress.find(&*iter);
    ASSERT_NE(addressAndID, idsByAddress.end()) << "Particle " << iter->getID() << " was moved.";
    EXPECT_EQ(iter->getID(), addressAndID->second);
    EXPECT_EQ(iter->isDummy(), deletedIDs.contains(iter->getID()));
  }
  EXPECT_EQ(numParticlesFound, numParticles);
  EXPECT_EQ(container->getNumberOfParticles(autopas::IteratorBehavior::owned), numParticles - deletedIDs.size());

  const auto leavingParticles = container->updateContainer(false);
  EXPECT_TRUE(leavingParticles.empty());
  EXPECT_EQ(container->size(), numParticles - deletedIDs.size());
}

/**
 * Adds particles while an owned particle has already left the box. This mimics what can happen at a rebuild, when
 * AutoPas adds the particles from its buffer. The oddly specific scenario comes from this being a regression test for
 * bugs which used to occur for this scenario. Checks that adding neither throws nor loses a particle and that
 * updateContainer(false) returns the particle that left.
 */
TEST_P(AllContainersTests, testAddParticleWhileOwnedParticleLeftBox) {
  const auto containerConfig = std::get<0>(GetParam());
  auto container = getInitializedContainer<ParticleFP64>(containerConfig);

  constexpr unsigned long numInitialParticles = 10;
  for (unsigned long id = 0; id < numInitialParticles; ++id) {
    container->addParticle(ParticleFP64({0.05 + 0.1 * static_cast<double>(id), 0.5, 0.5}, {0., 0., 0.}, id));
  }
  // Particle 0 moves just outside the box, as if it moved there since the last container update. The address of
  // particle 1 is stored to check below that the adds reallocated the particle vector of LinkedCellsReferences.
  std::uintptr_t addressOfParticle1{0};
  for (auto iter = container->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    if (iter->getID() == 0) {
      iter->setR({boxMin[0] - 0.1 * skin, 0.5, 0.5});
    }
    if (iter->getID() == 1) {
      addressOfParticle1 = reinterpret_cast<std::uintptr_t>(&*iter);
    }
  }
  // Add more particles than are stored, so a container that stores all particles in one std::vector must reallocate it.
  constexpr unsigned long numAddedParticles = numInitialParticles + 1;
  for (unsigned long id = numInitialParticles; id < numInitialParticles + numAddedParticles; ++id) {
    EXPECT_NO_THROW(
        container->addParticle(ParticleFP64({1. + 0.1 * static_cast<double>(id), 1., 1.}, {0., 0., 0.}, id)));
  }
  // A reallocation copies every particle to new memory, so particle 1 must now be stored at a different address.
  if (containerConfig.container == autopas::ContainerOption::linkedCellsReferences) {
    for (auto iter = container->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      if (iter->getID() == 1) {
        EXPECT_NE(reinterpret_cast<std::uintptr_t>(&*iter), addressOfParticle1)
            << "The particle vector was not reallocated, so this test does not cover the scenario it is meant for.";
      }
    }
  }

  const auto leavingParticles = container->updateContainer(false);
  ASSERT_EQ(leavingParticles.size(), 1);
  EXPECT_EQ(leavingParticles.front().getID(), 0);
  EXPECT_EQ(container->getNumberOfParticles(autopas::IteratorBehavior::owned),
            numInitialParticles + numAddedParticles - 1);
}