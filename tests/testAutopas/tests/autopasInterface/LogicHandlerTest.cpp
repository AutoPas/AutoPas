/**
 * @file LogicHandlerTest.cpp
 * @author Manish
 * @date 13.05.24
 */

#include "LogicHandlerTest.h"

#include "autopas/LogicHandler.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/ArbitraryConfigurations.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;
using ::testing::AtLeast;
using ::testing::Return;

void LogicHandlerTest::initLogicHandler() {
  const autopas::LogicHandlerInfo logicHandlerInfo{
      .boxMin{0., 0., 0.},
      .boxMax{10., 10., 10.},
      .cutoff = 2.5,
      .verletSkin = 0.5,
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
}

#ifdef AUTOPAS_DYNAMIC_CONTAINERS_FPB
/**
 * Tests that a particle after moving more than half skin distance is added to the fast particle buffer without
 * triggering a rebuild.
 */
TEST_F(LogicHandlerTest, testOneParticleForDynamicRebuild) {
  initLogicHandler();
  Molecule p1({0.5, 1., 1.}, {0., 0., 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  auto &container = _logicHandler->getContainer();
  std::array<double, 3> moveVec{0, container.getVerletSkin() * 0.1, 0};
  std::array<double, 3> moveVec_half{0, container.getVerletSkin() * 0.5, 0};

  // 0 ITERATION
  EXPECT_EQ(container.getNumberOfParticles(), 1) << "Only one particle has been added \n";

  // The particle should stay in the container
  auto leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 0 \n";

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

  // 1 ITERATION
  EXPECT_EQ(container.getNumberOfParticles(), 1) << "The particle stays in the container in iteration 1 \n";

  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 1 \n";

  // Buffer should remain empty because the particle has not moved far enough
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0)
      << "No particles in the buffer in iteration 1 because the particle has not moved far enough \n";

  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "No rebuild is required in iteration 1 \n";
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << "No particles in the buffer => decision was not triggered, initial value is false. \n";

  // no rebuild
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // 2 ITERATION
  // Particle moves more than half the skin
  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec_half);
  }
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added \n";

  // The particle has moved more than half a skin and enters the buffer
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays within the simulation boundary in iteration 2 \n";

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1) << "Fast particle is added into the buffer.\n";
  EXPECT_EQ(container.getNumberOfParticles(), 0) << "Container is empty as the particle has moved to the buffer. \n";

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

  // 3 ITERATION

  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0)
      << "The particle stays in the container as it has not moved since previous iteration. \n";
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << "Dynamic rebuild is false (either by default, or after rebuild due to fast particle in previous step) \n";
  // Tuning sample count is 3, so tuner should trigger a rebuild here.
  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Invalid because tuning triggers rebuild \n";
}

/**
 * Tests dynamic rebuild functionality for one particle moving across the periodic boundary and being added back to the
 * buffer. Additionally, the estimation method for predicting particles in the buffer due to migrating particles is
 * tested.
 */
TEST_F(LogicHandlerTest, testParticleInContainerMoveAcrossPeriodicBoundaryForDynamicRebuild) {
  initLogicHandler();
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

#endif

#ifdef AUTOPAS_DYNAMIC_CONTAINERS_STANDARD
/**
 * Tests dynamic rebuild functionalities for one particle case.
 * Dynamic rebuild should be triggered only when particle moves more than skin/2
 */
TEST_F(LogicHandlerTest, testOneParticleForDynamicRebuild) {
  initLogicHandler();
  Molecule p1({0.5, 1., 1.}, {0., 0., 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  auto &container = _logicHandler->getContainer();
  std::array<double, 3> moveVec{0, container.getVerletSkin() * 0.3, 0};

  // In the beginning, dynamic build is not required, so we expect false
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << " Particle has not moved yet, so no dynamic rebuild required. \n";

  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "Only one particle has been added to the container. \n";

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }

  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "Only one particle should be in the container. \n";
  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << " Particle has moved but not more than half the skin, so dynamic rebuild is not required. \n";

  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }

  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_TRUE(_logicHandler->getDoDynamicRebuild())
      << " Particle has moved more than half the skin, so dynamic rebuild is required. \n";
}
/**
 * Tests that the Velocity Method correctly estimates the rebuild frequency using the formula skin/vmax/deltaT/2 + 1.
 * First we test if the method works with one particle. We then add two more particles and test it again.
 */
TEST_F(LogicHandlerTest, testVelocityMethod) {
  initLogicHandler();
  Molecule p1({0.5, 1., 1.}, {0., 0.2, 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  double skin = 2.;
  double deltaT = 1.;
  double velocityMethodRFEstimate = _logicHandler->getVelocityMethodRFEstimate(skin, deltaT);
  // only one particle in the container --> vmax = 0.2
  // set vmax according to the maximum velocity of the particles that were added to the container before
  double vmax = 0.2;
  EXPECT_NEAR(velocityMethodRFEstimate, skin / vmax / deltaT / 2, 5);

  // test for multiple particles case
  Molecule p2({1.5, 1., 1.}, {0., 0.4, 0.3}, 0, 0);
  Molecule p3({1.5, 2., 2.}, {0., 0., 0.3}, 0, 0);
  _logicHandler->addParticle(p2);
  _logicHandler->addParticle(p3);
  // test the method with different skin length
  skin = 10.;
  velocityMethodRFEstimate = _logicHandler->getVelocityMethodRFEstimate(skin, deltaT);
  // 3 particles in the container --> vmax = (0.3*0.3 + 0.4*0.4)^0.5 = 0.5
  // set vmax according to the maximum velocity of the particles that were added to the container before
  vmax = 0.5;
  EXPECT_NEAR(velocityMethodRFEstimate, skin / vmax / deltaT / 2, 10);
}

/**
 * Tests dynamic rebuild functionalities for one particle moving across the periodic boundary and added to the container
 * when entering from other side. Particle is added to the container when rebuild is expected and hence should not
 * affect simulation pipeline anyway.
 */
TEST_F(LogicHandlerTest, testParticleInContainerMoveAcrossPeriodicBoundaryForDynamicRebuild) {
  initLogicHandler();
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
 * Tests dynamic rebuild functionalities for one particle moving across the periodic boundary and added to the buffer
 * when entering from other side. Particle is added to the buffer as the neighbor lists for the container is valid.
 * Dynamic rebuild should not be triggered even when particle moves more than skin/2 as it is now in buffer and hence
 * should not affect simulation pipeline anyway.
 */
TEST_F(LogicHandlerTest, testParticleInBufferMoveAcrossPeriodicBoundaryForDynamicRebuild) {
  initLogicHandler();
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
#endif