/**
 * @file LogicHandlerTest.cpp
 * @author Manish
 * @date 13.05.24
 */

#include "LogicHandlerTest.h"

#include "autopas/LogicHandler.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
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
      {{autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c08,
        autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled,
        autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec},
       {autopas::ContainerOption::linkedCells, cellSizeFactor, autopas::TraversalOption::lc_c18,
        autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled,
        autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec}});
  _tuningManager = std::make_shared<autopas::TuningManager>(autoTunerInfo);
  _tuningManager->addAutoTuner(
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""),
      autopas::InteractionTypeOption::pairwise);
  _logicHandler =
      std::make_unique<autopas::LogicHandler<Molecule>>(_tuningManager, logicHandlerInfo, verletRebuildFrequency, "");
}

#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
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