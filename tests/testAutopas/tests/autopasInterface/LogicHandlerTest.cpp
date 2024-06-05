/**
 * @file LogicHandlerTest.cpp
 * @author Manish
 * @date 13.05.24
 */

#include "LogicHandlerTest.h"

#include "autopas/LogicHandler.h"
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
        autopas::LoadEstimatorOption::none, autopas::DataLayoutOption::aos, autopas::Newton3Option::enabled}});
  autopas::AutoTuner autoTuner(tuningStrategies, searchSpace, autoTunerInfo, verletRebuildFrequency, "");
  _logicHandler =
      std::make_shared<autopas::LogicHandler<Molecule>>(autoTuner, logicHandlerInfo, verletRebuildFrequency, "");

  /*
  auto leavingParticles = _logicHandler.updateContainer();

  std::cout <<"After update\n";
  std::cout << leavingParticles.size() << "\n";
  std::cout << _logicHandler.getNeighborListsInvalidDoDynamicRebuild() << "\n";
  //std::cout << _logicHandler.neighborListsAreValid() << "\n";
  std::cout << _logicHandler.getContainer().getNumberOfParticles() << "\n";

  //std::cout << _logicHandler.getNeighborListsInvalidDoDynamicRebuild() << "\n";
  ASSERT_EQ(_logicHandler.getNeighborListsInvalidDoDynamicRebuild(), true) << "The neighbor list should require rebuild
  now as the particle has moved more than half the skin \n";
*/
}

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

TEST_F(LogicHandlerTest, testParticleMoveAcrossPeriodicBoundaryForDynamicRebuild) {
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

  // In the beginning, dynamic build is not required, so we expect false
  ASSERT_FALSE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " Particle has not moved yet, so no dynamic rebuild required. \n";

  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1) << "Only one particle has been added \n";
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
  EXPECT_EQ(_logicHandler->getContainer().getNumberOfParticles(), 1)
      << "One particle added on the other side of periodic boundary \n";

  _logicHandler->resetNeighborListsInvalidDoDynamicRebuild();
  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_TRUE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " Particle has moved across the periodic boundary but not more than half the skin, so dynamic rebuild is not "
         "required. But because we have not shifted the rAtRebuild, we expect otherwise. \n";

  // adjusting rAtRebuild for periodic boundary shift
  for (auto iter = _logicHandler->begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    iter->addRAtRebuild(shiftVecPeriodicY);
  }

  _logicHandler->resetNeighborListsInvalidDoDynamicRebuild();
  _logicHandler->checkNeighborListsInvalidDoDynamicRebuild();
  ASSERT_FALSE(_logicHandler->getNeighborListsInvalidDoDynamicRebuild())
      << " Particle has moved across the periodic boundary but not more than half the skin, and rAtRebuild is also "
         "shifted, so dynamic rebuild is not required. \n";
}