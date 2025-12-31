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
        autopas::InteractionTypeOption::pairwise},
       {autopas::ContainerOption::verletLists, cellSizeFactor, autopas::TraversalOption::lc_c08,
        autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled,
        autopas::InteractionTypeOption::pairwise}});
  _tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  _logicHandler =
      std::make_unique<autopas::LogicHandler<Molecule>>(_tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
}

#if defined(AUTOPAS_ENABLE_DYNAMIC_CONTAINERS) && !defined(AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER_LIN)
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
  ASSERT_FALSE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " Particle has not moved yet, so no dynamic rebuild required. \n";

  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added \n";

  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }
  _logicHandler->resetNeighborListsInvalidDoDynamicRebuild();
  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_FALSE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " Particle has moved but not more than half the skin, so dynamic rebuild is not required. \n";

  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }
  _logicHandler->resetNeighborListsInvalidDoDynamicRebuild();
  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_TRUE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " Particle has moved more than half the skin, so dynamic rebuild is required. \n";
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

  Molecule p1({0.5, boxMaxY - skin * 0.15, 1.}, {0., 0., 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added \n";

  // In the beginning, dynamic build is not required, so we expect false
  ASSERT_FALSE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " Particle has not moved yet, so no dynamic rebuild required. \n";

  // In the first step, neighborListsAreValid is false
  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "In the first iteration, neighbor lists are invalid.";

  // 0.3 skin + (boxMaxY - 0.15 skin) = boxMaxY + 0.15 skim -> particle is outside the boundary
  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }
  auto leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 1) << "Exactly one particle has left the container \n";
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 0) << "No particle left in the container \n";

  // shifting particle position to replicate periodic boundary effect
  for (auto particle : leavingParticles) {
    particle.addR(shiftVecPeriodicY);
    _logicHandler->addParticle(particle);
  }
  // As neighbor lists are invalid, particle is added to the container directly, as it will be rebuilt soon
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "One particle added on the other side of periodic boundary \n";

  _logicHandler->resetNeighborListsInvalidDoDynamicRebuild();
  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_TRUE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " Particle has moved across the periodic boundary but not more than half the skin, so dynamic rebuild is not "
         "required. But because we have not shifted the rAtRebuild, we expect otherwise. \n";

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);
  functor.setParticleProperties(24.0, 1);

  // iterate once so that neighbor lists are valid
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "No particle has the container \n";
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // After one iteration, neighbor lists are rebuilt, and neighborListsAreValid is false
  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "After one iteration, neighbor lists are valid.";

  _logicHandler->resetNeighborListsInvalidDoDynamicRebuild();
  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_FALSE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " The neighbor list is rebuilt as previously neighbor lists were invalid. rAtRebuild is reset at the same "
         "time.\n";
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
  ASSERT_FALSE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " Particle has not moved yet, so no dynamic rebuild required. \n";

  // In the first step, neighborListsAreValid is false
  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "In the first iteration, neighbor lists are invalid.";

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);
  functor.setParticleProperties(24.0, 1);

  // iterate once so that neighbor lists are valid
  auto leavingParticles = _logicHandler->updateContainer();
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // After one iteration, neighbor lists are rebuilt, and neighborListsAreValid is false
  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "After one iteration, neighbor lists are valid.";

  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_FALSE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
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

  EXPECT_EQ(_logicHandler->getNumberOfParticlesOwned(), 1)
      << "One particle added on the other side of periodic boundary \n";

  _logicHandler->resetNeighborListsInvalidDoDynamicRebuild();
  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_FALSE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " Particle has moved across the periodic boundary and more than half the skin, but as it is in the buffer, it "
         "doesn't affect the container. \n";
}
#endif

#ifdef AUTOPAS_ENABLE_FAST_PARTICLE_BUFFER_LIN
/**
 * Tests dynamic rebuild functionalities for one particle case.
 * Dynamic rebuild should be triggered only when particle moves more than skin/2
 */
TEST_F(LogicHandlerTest, testOneParticleForDynamicRebuild) {
  initLogicHandler();
  Molecule p1({0.5, 1., 1.}, {0., 0., 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  auto &container = _logicHandler->getContainer();
  std::array<double, 3> moveVec{0, container.getVerletSkin() * 0.1, 0};
  std::array<double, 3> moveVec_half{0, container.getVerletSkin() * 0.5, 0};

  // 0 ITERATION
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added \n";

  // sample tuning is 3
  ASSERT_FALSE(_logicHandler->isTuningInNeedOfRebuild());

  auto leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0);

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0);

  // In the beginning, dynamic build is not required, so we expect false
  ASSERT_FALSE(_logicHandler->neighborListsAreValidBuf()) << " First Iteration rebuild is required. \n";
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild()) << " No decision was triggered initial value false. \n";

  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec);
  }

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);
  // rebuild
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // 1 ITERATION
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1);

  // sample tuning is 3
  ASSERT_FALSE(_logicHandler->isTuningInNeedOfRebuild());

  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0);

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0);

  // In the beginning, dynamic build is not required, so we expect false
  ASSERT_TRUE(_logicHandler->neighborListsAreValidBuf());
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << " No particles in buffer => decision wasn't triggered initial value false. \n";

  // no rebuild
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // 2 ITERATION
  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addR(moveVec_half);
  }
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1);

  // sample tuning is 3
  ASSERT_FALSE(_logicHandler->isTuningInNeedOfRebuild());

  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "Exactly one particle has left the container \n";

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0);
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1);

  ASSERT_TRUE(_logicHandler->getDoDynamicRebuild()) << " Need for first rebuild time estimate triggers rebuild. \n";
  ASSERT_FALSE(_logicHandler->neighborListsAreValidBuf());

  // rebuild
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // 3 ITERATION

  // sample tuning is 3, so here starts tuning iteration
  ASSERT_TRUE(_logicHandler->neighborListsAreValidBuf());
  ASSERT_TRUE(_logicHandler->isTuningInNeedOfRebuild());
  ASSERT_FALSE(_logicHandler->neighborListsAreValidBuf());
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild());
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

  Molecule p1({0.5, boxMaxY - skin * 0.15, 1.}, {0., 0., 0.}, 0, 0);
  _logicHandler->addParticle(p1);
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added \n";

  // In the first step, neighborListsAreValid is false
  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "In the first iteration, neighbor lists are invalid.";
  auto leavingParticles = _logicHandler->updateContainer();

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);
  functor.setParticleProperties(24.0, 1);

  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  ASSERT_TRUE(_logicHandler->neighborListsAreValid());

  // 0.3 skin + (boxMaxY - 0.15 skin) = boxMaxY + 0.15 skim -> particle is outside the boundary
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

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1);
  ASSERT_TRUE(_logicHandler->neighborListsAreValid());
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // iterate once so that neighbor lists are valid
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0);
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 0);
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1);
  // no particles in the buffer
  ASSERT_FALSE(_logicHandler->neighborListsAreValid());

  EXPECT_EQ(_logicHandler->getNumParticlesBufferEstimate(), 2);
}

#endif