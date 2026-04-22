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
       {autopas::ContainerOption::verletLists, cellSizeFactor, autopas::TraversalOption::vl_list_iteration,
        autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled,
        autopas::InteractionTypeOption::pairwise}});
  _tunerMap.emplace(
      autopas::InteractionTypeOption::pairwise,
      std::make_unique<autopas::AutoTuner>(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, ""));
  _logicHandler =
      std::make_unique<autopas::LogicHandler<Molecule>>(_tunerMap, logicHandlerInfo, verletRebuildFrequency, "");
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
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added \n";

  // Test method isTuningInNeedOfRebuild, which should return false in iteration 0
  ASSERT_FALSE(_logicHandler->isTuningInNeedOfRebuild()) << "No tuning rebuild in iteration 0 \n";

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
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "The particle stays in the container in iteration 1 \n";

  // Sample tuning count is 3
  ASSERT_FALSE(_logicHandler->isTuningInNeedOfRebuild())
      << "No tuning rebuild in iteration 1 because 3 tuning samples are required \n";

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

  // Sample tuning count is 3
  ASSERT_FALSE(_logicHandler->isTuningInNeedOfRebuild())
      << "No tuning rebuild in iteration 2 because 3 tuning samples are required \n";

  // The particle has moved more than half a skin and enters the buffer
  leavingParticles = _logicHandler->updateContainer();
  EXPECT_EQ(leavingParticles.size(), 0) << "The particle stays in the container in iteration 2 \n";

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1) << "Fast particle is added into the buffer.\n";
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 0)
      << "Container is empty as no rebuild is triggered. \n";

  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Neighbor lists remain valid. \n";

  // No rebuild
  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  // 3 ITERATION

  // Sample tuning count is 3, so tuning starts here
  ASSERT_TRUE(_logicHandler->isTuningInNeedOfRebuild()) << "Rebuild because of tuning after 3 samples \n";
  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Invalid because tuning triggers rebuild \n";
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild())
      << "Dynamic rebuild is false by default, because rebuild happens because of tuning \n";
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

  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Iteration 0 requires a rebuild. \n";
  auto leavingParticles = _logicHandler->updateContainer();

  constexpr double cutoff = 1.1;
  LJFunctorGlobals functor(cutoff);
  functor.setParticleProperties(24.0, 1);

  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  ASSERT_TRUE(_logicHandler->neighborListsAreValid()) << "No rebuild in iteration 1 \n";

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

  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1) << "Migrating particle in the buffer";

  _logicHandler->computeInteractionsPipeline(&functor, autopas::options::InteractionTypeOption::pairwise);

  leavingParticles = _logicHandler->updateContainer();

  // Buffer contains one migrating particle and neighbor lists remain valid.
  ASSERT_FALSE(_logicHandler->getDoDynamicRebuild()) << "Particle stays in buffer without triggering rebuild. \n";

  ASSERT_FALSE(_logicHandler->neighborListsAreValid()) << "Neighbor lists remain valid. \n";

  EXPECT_EQ(leavingParticles.size(), 0);
  EXPECT_EQ(_logicHandler->getNumberOfParticlesBuffer(), 1)
      << "Particle stays in buffer without triggering rebuild. \n";
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 0)
      << "Rebuild is not triggered, so container remains empty. \n";

  // The estimated number of particles in the buffer is 2 because, at the point of decision-making, the buffer contains
  // 1 particle, and we expect the buffer filling to increase by one additional particle due to migration, as observed
  // in the previous iteration.
  EXPECT_EQ(_logicHandler->getNumParticlesBufferEstimate(), 2)
      << "The expected estimate of particles in the buffer is 2 \n";
}

#endif