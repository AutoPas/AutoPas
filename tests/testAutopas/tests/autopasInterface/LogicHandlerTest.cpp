/**
 * @file LogicHandlerTest.cpp
 * @author Manish
 * @date 13.05.24
 */

#include "LogicHandlerTest.h"

#include "autopas/LogicHandler.h"
#include "autopas/options/IteratorBehavior.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/ArbitraryConfigurations.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Return;

size_t LogicHandlerTest::initLogicHandler() {
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
      .cutoff = 2.5,
      .verletSkin = 0.5,
      .deltaT = 1,
  };
  const autopas::AutoTunerInfo autoTunerInfo{
      .tuningInterval = 10000,
      .maxSamples = 3,
  };
  autopas::AutoTuner::TuningStrategiesListType tuningStrategies{};
  constexpr double cellSizeFactor = 1.;
  constexpr unsigned int verletRebuildFrequency = 10;
  const std::set<autopas::Configuration> searchSpace(
      {arbitraryConfigurations::_arbitrary_config_2B_2, arbitraryConfigurations::_arbitrary_config_2B_6});
  _tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  _tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);
  _logicHandler = std::make_unique<autopas::LogicHandler<Molecule>>(
      _tuningManager, logicHandlerInfo, verletRebuildFrequency, "", autoTunerInfo.aosSortingThreshold,
      autoTunerInfo.soaSortingThreshold);
  // Depending on what type of neighbor list rebuilding is used.
#ifdef AUTOPAS_DYNAMIC_CONTAINERS_FPB
  return 1;
#elif defined(AUTOPAS_DYNAMIC_CONTAINERS_STANDARD)
  return 2;
#else
  return 0;
#endif
}

/**
 * Tests that a particle after moving more than half skin distance is added to the fast particles buffer (FPB). We also
 * check that a halo particle (corresponding to an owned particle in neighboring domain) is added to the FPB when it has
 * moved more than halfSkin.
 */
TEST_F(LogicHandlerTest, testOneParticleAddedtoFPB) {
  auto usingDynamicRebuilding = initLogicHandler();
  if (not(usingDynamicRebuilding == 1)) {
    GTEST_SKIP_("Using static or standard dynamic rebuilding.");
  }
  Molecule p1({0.5, 1., 1.}, {0., 0., 0.}, 0, 0);
  Molecule hP1({-0.5, 1., 1.}, {0., 0., 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  auto &container = _logicHandler->getContainer();
  std::array<double, 3> moveVec{0, container.getVerletSkin() * 0.1, 0};
  std::array<double, 3> moveVecHalf{0, container.getVerletSkin() * 0.5, 0};
  std::array<double, 3> moveVecCloseToSkin{0, container.getVerletSkin() * 0.48, 0};

  EXPECT_EQ(container.getNumberOfParticles(), 1) << "Only one owned particle has been added. \n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesOwned(), 1) << "Only one owned particle has been added. \n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0)
      << "Before the first update container call, all the particles would be added to the container.\n";

  // Iteration 0
  // The particle should stay in the container
  auto leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 0 \n";

  // UpdateContainer call removed the halo particles.
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 1)
      << "Halo particle is removed during updateContainer() call. \n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesHalo(), 0)
      << "Halo particle is removed during updateContainer() call. \n";

  // Halo particles is inserted.
  _logicHandler->addHaloParticle(hP1);
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 2)
      << "Halo particle is inserted in the container.\n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesHalo(), 1) << "Halo particle is inserted. \n";

  // Buffer should remain empty because the particle has not moved
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0)
      << "No particles in the buffer in iteration 0 because there was no movement \n";

  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Iteration 0 requires a rebuild. \n";

  // At the beginning, dynamic rebuild is not required, so false is expected
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << "No particles in the buffer => decision was not triggered, initial value is false. \n";

  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);
  // rebuild
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 1
  EXPECT_EQ(container.getNumberOfParticles(), 1) << "The particle stays in the container in iteration 1 \n";

  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 1 \n";

  EXPECT_EQ(_logicHandler->getNumberOfParticlesHalo(), 0)
      << "Halo particle is removed during updateContainer() call. \n";
  _logicHandler->addHaloParticle(hP1);
  EXPECT_EQ(_logicHandler->getNumberOfParticlesHalo(), 1) << "Halo particle is reinserted. \n";

  // Buffer should remain empty because the particle has not moved far enough
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0)
      << "No particles in the buffer in iteration 1 because the particle has not moved far enough \n";

  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "No rebuild is required in iteration 1 \n";
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << "No particles in the buffer => decision was not triggered, initial value is false. \n";

  // no rebuild
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 2
  // Particle moves more than half the skin
  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVecHalf);
  }
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added \n";

  // The particle has moved more than half a skin and enters the buffer
  leavingParticles = _logicHandler->updateContainer();

  EXPECT_EQ(_logicHandler->getNumberOfParticlesHalo(), 0)
      << "Halo particle is removed during updateContainer() call. \n";
  // halo particle moved 0.1 skin distance
  hP1.addR(moveVec);
  _logicHandler->addHaloParticle(hP1);
  EXPECT_EQ(_logicHandler->getNumberOfParticlesHalo(), 1) << "Halo particle is reinserted. \n";

  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays within the simulation boundary in iteration 2 \n";

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1) << "Fast particle is added into the buffer.\n";
  EXPECT_EQ(container.getNumberOfParticles(), 0) << "Container is empty as the particle has moved to the buffer. \n";
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 1)
      << "Halo particle is reinserted in the container. \n";

  auto bufferTriggersRebuild = _logicHandler->getDoDynamicRebuild();
  if (bufferTriggersRebuild) {
    ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Neighbor lists will be rebuilt.\n";
  } else {
    ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << " Neighbor lists remain valid. \n";
  }

  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  if (not bufferTriggersRebuild) {
    EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1) << "Fast particle remains in the buffer.\n";
  } else {
    EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0)
        << "Fast particle is moved back to container upon rebuilt.\n";
  }

  // Iteration 3
  leavingParticles = _logicHandler->updateContainer();

  EXPECT_EQ(_logicHandler->getNumberOfParticlesHalo(), 0)
      << "Halo particle is removed during updateContainer() call. \n";
  // halo particle moved 0.5 skin distance
  hP1.addR(moveVecHalf);
  _logicHandler->addHaloParticle(hP1);
  EXPECT_EQ(_logicHandler->getNumberOfParticlesHalo(), 1) << "Halo particle is reinserted. \n";
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 2)
      << "Halo particle is reinserted in the container as we are about to rebuild. \n";

  EXPECT_EQ(leavingParticles.size(), 0)
      << "The particle stays in the container as it has not moved since previous iteration. \n";
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << "Dynamic rebuild is false (either by default, or after rebuild due to fast particle in previous step) \n";
  // Tuning sample count is 3, so tuner should trigger a rebuild here.
  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Invalid because tuning triggers rebuild \n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 4
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0)
      << "There should be no particles in the buffer, as the list was rebuilt in previous iteration. \n";
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 2)
      << "One owned particle and old halo particle exist. \n";

  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 1)
      << "Halo particle is cleared, only the owned particle remain. \n";
  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "Rebuild should not be expected. \n";

  // halo particle has moved to closed to the halfSkin but not yet farther than that.
  hP1.addR(moveVecCloseToSkin);
  _logicHandler->addHaloParticle(hP1);

  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 2)
      << "Halo particle is reinserted into the container. \n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0) << "Halo particle is in the container. \n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 5
  // we should detect that next time halo will be further away than halfskin based on movement until previous iteration
  // and ensure that it will be added to the buffer.
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 1)
      << "Halo particle is cleared, only the owned particle remain. \n";
  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "Rebuild should not be expected. \n";

  hP1.addR(moveVec);
  _logicHandler->addHaloParticle(hP1);

  // updateContainer call should have deleted the previous halo particle from container with a different ID. New halo
  // should be added to the buffer.
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 1)
      << "Halo particle cannot be inserted into the container using updateHaloParticle() call. \n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1) << "Halo particle is inserted in the buffer. \n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);
}

/**
 * Tests that a particle moving across the periodic boundary is added to the buffer. As the rebuild decision when using
 * Fast Particles Buffer (FPB) is made during the updateContainer() call but migrating particles enter the domain after
 * that, the estimation method for predicting particles in the buffer due to migrating particles (based on data from
 * previous iteration) is tested.
 */
TEST_F(LogicHandlerTest, testParticleInContainerMoveAcrossPeriodicBoundaryAddedtoBuffer) {
  auto usingDynamicRebuilding = initLogicHandler();
  if (not(usingDynamicRebuilding == 1)) {
    GTEST_SKIP_("Using static or standard dynamic rebuilding.");
  }
  auto boxMaxY = _logicHandler->getContainer().getBoxMax()[1];
  auto boxMinY = _logicHandler->getContainer().getBoxMin()[1];
  auto &container = _logicHandler->getContainer();
  auto skin = container.getVerletSkin();
  std::array<double, 3> moveVec{0, skin * 0.3, 0};
  // periodic boundary shift
  std::array<double, 3> shiftVecPeriodicY{0, boxMinY - boxMaxY, 0};

  Molecule p1({0.5, boxMaxY - skin * 0.15, 1.}, {0., 0., 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added \n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0) << "No particles in the buffer yet.\n";

  // Iteration 0
  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Iteration 0 requires a rebuild. \n";
  auto leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "No particle movement at the start. \n";

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);
  functor.setParticleProperties(24.0, 1);

  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 1
  // 0.3 skin + (boxMaxY - 0.15 skin) = boxMaxY + 0.15 skim -> particle is outside the boundary
  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 1) << " Exactly one particle has moved outside the periodic boundary. \n";
  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "No rebuild in iteration 1 \n";
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 0) << "No particle left in the container \n";

  // shifting particle position to replicate periodic boundary effect
  for (auto particle : leavingParticles) {
    particle.addR(shiftVecPeriodicY);
    _logicHandler->addParticle(particle);
  }

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1) << "Particle is added into the buffer.\n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 2

  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << " No particles leave the simulation domain in this iteration. \n";

  // Buffer contains one migrating particle and neighbor lists remain valid.
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild()) << "Particle stays in buffer without triggering rebuild. \n";
  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "Neighbor lists remain valid. \n";

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1)
      << "Particle stays in buffer without triggering rebuild. \n";
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 0)
      << "Rebuild is not triggered, so container remains empty. \n";

  // The estimated number of particles in the buffer is 2 because, at the point of decision-making, the buffer contains
  // 1 particle, and we expect the buffer filling to increase by one additional particle due to migration, as observed
  // in the previous iteration.
  EXPECT_EQ(_logicHandler->getNumParticlesBufferEstimate(), 2)
      << "The expected estimate of particles in the buffer is 2 \n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 3
  // Rebuilding due to tuning
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << " No particles leave the simulation domain in this iteration. \n";

  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Neighbor lists remain valid. \n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0) << "Particle moved back to container during rebuild. \n";
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "Particle moved back to container during rebuild. \n";
}

/**
 * Tests standard dynamic rebuilding functionalities for one particle case.
 * Dynamic rebuild should be triggered only when particle moves more than skin/2.
 * As the information about movement of haloParticles is not available during the updateContainer() call where we check
 * the displacements, we simply force it to the buffer if it was close to halfSkin displacement in the previous
 * iteration.
 */
TEST_F(LogicHandlerTest, testOneParticleForDynamicRebuild) {
  auto usingDynamicRebuilding = initLogicHandler();
  if (not(usingDynamicRebuilding == 2)) {
    GTEST_SKIP_("Currently using the static or dynamic rebuilding with fast particles buffer (FPB) approach.");
  }
  Molecule p1({0.5, 1., 1.}, {0., 0., 0.}, 0, 0);
  Molecule hP1({-0.5, 1., 1.}, {0., 0., 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  auto &container = _logicHandler->getContainer();

  std::array<double, 3> moveVec{0, container.getVerletSkin() * 0.3, 0};
  std::array<double, 3> moveVecHalf{0, container.getVerletSkin() * 0.5, 0};
  std::array<double, 3> moveVecCloseToSkin{0, container.getVerletSkin() * 0.48, 0};

  EXPECT_EQ(container.getNumberOfParticles(), 1) << "Only one owned particle has been added. \n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesOwned(), 1) << "Only one owned particle has been added. \n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0)
      << "Before the first update container call, all the particles would be added to the container.\n";

  // Iteration 0
  // The particle should stay in the container
  auto leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 0 \n";

  // UpdateContainer call removed the halo particles.
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 1)
      << "Halo particle is removed during updateContainer() call. \n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesHalo(), 0)
      << "Halo particle is removed during updateContainer() call. \n";

  // Halo particles is inserted.
  _logicHandler->addHaloParticle(hP1);
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 2)
      << "Halo particle is inserted in the container.\n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesHalo(), 1) << "Halo particle is inserted. \n";

  // Buffer should remain empty because the particle has not moved
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0)
      << "No particles in the buffer in iteration 0 because there was no movement \n";

  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Iteration 0 requires a rebuild. \n";

  // At the beginning, dynamic rebuild is not required, so false is expected
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild()) << "No particle has moved, so no rebuilding is required. \n";

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 1

  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 0 \n";

  // Halo particles is reinserted.
  _logicHandler->addHaloParticle(hP1);

  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "Only one particle should be in the container. \n";
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 2)
      << "Halo particle is reinserted in the container.\n";

  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << "Particle has moved but not more than half the skin, so dynamic rebuild is not required. \n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 2

  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVecHalf);
  }
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 0 \n";

  // Halo particles is reinserted.
  _logicHandler->addHaloParticle(hP1);

  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "Only one particle should be in the container. \n";
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 2)
      << "Halo particle is reinserted in the container.\n";

  ASSERT_TRUE(_logicHandler->getDoDynamicRebuild())
      << "Particle has moved more than half the skin, so dynamic rebuild is required. \n";

  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 3
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 0 \n";

  // Halo particles is reinserted.
  hP1.addR(moveVecHalf);
  _logicHandler->addHaloParticle(hP1);

  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "Only one particle should be in the container. \n";
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 2)
      << "Halo particle is reinserted in the container.\n";

  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << "Particle has moved but not more than half the skin, so dynamic rebuild is not required. \n";

  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Tuning triggers the rebuilding. \n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 4
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 0 \n";

  hP1.addR(moveVecCloseToSkin);
  _logicHandler->addHaloParticle(hP1);
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0) << "Halo particle is in the container. \n";
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "Only one particle should be in the container. \n";
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 2)
      << "Halo particle is reinserted in the container.\n";

  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << "Particle has moved but not more than half the skin, so dynamic rebuild is not required. \n";
  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "No rebuilding is expected. \n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Iteration 5
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 0 \n";

  hP1.addR(moveVec);
  _logicHandler->addHaloParticle(hP1);
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1)
      << "Halo particle should be added to the buffer during updateContainer() call, as it was close to halfSkin "
         "earlier. \n";
  EXPECT_EQ(container.getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo), 1)
      << "Halo particle is reinserted in the buffer.\n";
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << "Halo particle has moved more than half the skin, but it is now in a buffer. \n";
  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "No rebuilding is expected. \n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);
}

/**
 * Tests that the Velocity Method correctly estimates the rebuild frequency using the formula skin/vmax/deltaT/2 for one
 * particle system.
 */
TEST_F(LogicHandlerTest, testVelocityMethodForOneParticle) {
  auto usingDynamicRebuilding = initLogicHandler();
  if (usingDynamicRebuilding == 0) {
    GTEST_SKIP_("Currently using the static rebuilding approach with fixed rebuild frequency.");
  }
  Molecule p1({0.5, 1., 1.}, {0., 0.2, 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  // Iteration 0
  _logicHandler->updateContainer();
  // Iteration 1
  _logicHandler->updateContainer();
  // Iteration 2
  _logicHandler->updateContainer();

  double skin = 0.5;
  double deltaT = 1.;
  double velocityMethodRFEstimate = _logicHandler->getVelocityMethodRFEstimate(skin, deltaT);
  // only one particle in the container --> vmax = 0.2
  // set vmax according to the maximum velocity of the particles that were added to the container before
  double vmax = 0.2;
  EXPECT_NEAR(velocityMethodRFEstimate, skin / vmax / deltaT / 2, 1.25);
}

/**
 * Tests that the Velocity Method correctly estimates the rebuild frequency using the formula skin/vmax/deltaT/2 for
 * several particle system.
 */
TEST_F(LogicHandlerTest, testVelocityMethodForMultipleParticles) {
  auto usingDynamicRebuilding = initLogicHandler();
  if (usingDynamicRebuilding == 0) {
    GTEST_SKIP_("Currently using the static rebuilding approach with fixed rebuild frequency.");
  }
  Molecule p1({0.5, 1., 1.}, {0., 0.2, 0.}, 0, 0);
  Molecule p2({1.5, 1., 1.}, {0., 0.4, 0.3}, 0, 0);
  Molecule p3({1.5, 2., 2.}, {0., 0., 0.3}, 0, 0);
  _logicHandler->addParticle(p1);
  _logicHandler->addParticle(p2);
  _logicHandler->addParticle(p3);
  // Iteration 0
  _logicHandler->updateContainer();
  // Iteration 1
  _logicHandler->updateContainer();
  // Iteration 2
  _logicHandler->updateContainer();

  // test the method with different skin length
  auto skin = 0.5;
  auto deltaT = 1.;
  auto velocityMethodRFEstimate = _logicHandler->getVelocityMethodRFEstimate(skin, deltaT);

  // 3 particles in the container --> vmax = (0.3*0.3 + 0.4*0.4)^0.5 = 0.5
  // set vmax according to the maximum velocity of the particles that were added to the container before
  auto vmax = 0.5;
  EXPECT_NEAR(velocityMethodRFEstimate, skin / vmax / deltaT / 2, 0.5);
}

/**
 * Tests standard dynamic rebuilding functionalities for one particle moving across the periodic boundary and added to
 * the container when entering from other side. Particle is added to the container when rebuild is expected and hence
 * should not affect simulation pipeline anyway.
 */
TEST_F(LogicHandlerTest, testParticleInContainerMoveAcrossPeriodicBoundaryForDynamicRebuild) {
  auto usingDynamicRebuilding = initLogicHandler();
  if (not(usingDynamicRebuilding == 2)) {
    GTEST_SKIP_(
        "Currently using the static rebuilding or dynamic rebuilding with fast particles buffer (FPB) approach.");
  }
  auto boxMaxY = _logicHandler->getContainer().getBoxMax()[1];
  auto boxMinY = _logicHandler->getContainer().getBoxMin()[1];
  auto &container = _logicHandler->getContainer();
  auto skin = container.getVerletSkin();
  std::array<double, 3> moveVec{0, skin * 0.3, 0};
  // periodic boundary shift
  std::array<double, 3> shiftVecPeriodicY{0, boxMinY - boxMaxY, 0};

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);

  Molecule p1({0.5, boxMaxY - skin * 0.15, 1.}, {0., 0., 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added. \n";

  // In the beginning, dynamic build is not required, so we expect false
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << " Particle has not moved yet, so no dynamic rebuild required. \n";

  // At the start, neighborListsAreValid is false
  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "In the first iteration, neighbor lists are invalid.";

  // Moving particle outside the boundary
  // 0.3 skin + (boxMaxY - 0.15 skin) = boxMaxY + 0.15 skim -> particle is outside the boundary
  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }

  // First iteration - expecting a rebuild
  auto leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 1) << "No particles left the container. \n";
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 0) << "No particle left in the container \n";

  // shifting particle position to replicate periodic boundary effect
  for (auto particle : leavingParticles) {
    // boxMaxY + 0.15 skim - boxMaxY = 0.15 skim -> particle is inside the domain
    particle.addR(shiftVecPeriodicY);
    _logicHandler->addParticle(particle);
  }
  // As neighbor lists are invalid, particle is added to the container directly, as it will be rebuilt soon
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "One particle added on the other side of periodic boundary \n";

  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_TRUE(_logicHandler->getDoDynamicRebuild())
      << " Particle has moved across the periodic boundary but not more than half the skin, so dynamic rebuild is not "
         "required. But because we have not shifted the rAtRebuild, we expect otherwise. \n";

  // Rebuild neighbor lists to make them valid for the next iteration
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // Second iteration
  leavingParticles = _logicHandler->updateContainer();
  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "Neighbor lists were rebuilt in first iteration.";
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << " Displacement since last rebuild was reset during the rebuild in previous iteration and particle did not "
         "move in this iteration, so no dynamic rebuild is still not required. \n";
}

/**
 * Tests standard dynamic rebuilding functionalities for one particle moving across the periodic boundary and added to
 * the buffer when entering from other side. Particle is added to the buffer as the neighbor lists for the container is
 * valid. Dynamic rebuild should not be triggered even when particle moves more than skin/2 as it is now in buffer and
 * hence should not affect simulation pipeline anyway.
 */
TEST_F(LogicHandlerTest, testParticleInBufferMoveAcrossPeriodicBoundaryForDynamicRebuild) {
  auto usingDynamicRebuilding = initLogicHandler();
  if (not(usingDynamicRebuilding == 2)) {
    GTEST_SKIP_(
        "Currently using the static rebuilding or dynamic rebuilding with fast particles buffer (FPB) approach.");
  }
  auto boxMaxY = _logicHandler->getContainer().getBoxMax()[1];
  auto boxMinY = _logicHandler->getContainer().getBoxMin()[1];
  auto &container = _logicHandler->getContainer();
  auto skin = container.getVerletSkin();
  std::array<double, 3> moveVec{0, skin * 0.3, 0};
  // periodic boundary shift
  std::array<double, 3> shiftVecPeriodicY{0, boxMinY - boxMaxY, 0};

  Molecule p1({0.5, boxMaxY - skin * 0.15, 1.}, {0., 0., 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added \n";

  // In the beginning, dynamic build is not required, so we expect false
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << " Particle has not moved yet, so no dynamic rebuild required. \n";

  // In the first step, neighborListsAreValid is false
  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "In the first iteration, neighbor lists are invalid.";

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);
  functor.setParticleProperties(24.0, 1);

  // iterate once so that neighbor lists are valid
  auto leavingParticles = _logicHandler->updateContainer();
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // REMOVED: ASSERT_TRUE(_logicHandler->neighborListsAreValid())
  // Calling this mid-step evaluates requiresRebuilding(0), which is true, and falsely mutates the valid list!

  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << " The neighbor list is rebuilt as previously neighbor lists were invalid.\n";

  // 0.3 skin + (boxMaxY - 0.15 skin) = boxMaxY + 0.15 skim -> particle is outside the boundary and has travelled less
  // than skin/2
  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 1) << "Exactly one particle has left the container \n";
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 0) << "No particle left in the container \n";

  // shifting particle position to replicate periodic boundary effect
  for (auto particle : leavingParticles) {
    particle.addR(shiftVecPeriodicY);
    _logicHandler->addParticle(particle);
  }
  // As neighbor lists are valid, particle is added to the buffer, so container will have 0 particles
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 0)
      << "Particle added to the buffer, so container must still be empty. \n";
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1) << "Exactly one particle is now added to the buffer. \n";

  EXPECT_EQ(_logicHandler->getNumberOfParticlesOwned(), 1)
      << "One particle added on the other side of periodic boundary \n";

  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);
  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << " Particle has moved across the periodic boundary and more than half the skin, but as it is in the buffer, it "
         "doesn't affect the container. \n";
}