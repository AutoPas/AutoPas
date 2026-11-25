/**
 * @file RegressionModelsTest.cpp
 * @author Natallia Padalinskaya
 * @date 28.12.25
 */

#include "DynamicRebuildDecisionWithFPBTest.h"

using namespace autopas::utils::DynamicRebuildDecisionWithFPB;

TEST(DynamicRebuildDecisionWithFPBTest, testRebuildEstimateManager) {
  RebuildEstimateManager rebuildPredictor{};

  rebuildPredictor.addNewPoint(10);
  EXPECT_EQ(rebuildPredictor.predict(), 10);
  rebuildPredictor.addNewPoint(10);
  EXPECT_EQ(rebuildPredictor.predict(), 10);
  rebuildPredictor.addNewPoint(16);
  EXPECT_EQ(rebuildPredictor.predict(), 12);
  rebuildPredictor.addNewPoint(12);
  EXPECT_EQ(rebuildPredictor.predict(), 12);
  rebuildPredictor.addNewPoint(7);
  EXPECT_EQ(rebuildPredictor.predict(), 11);

  // Adding 6th sample, prediction will be average of last 5 samples only.
  rebuildPredictor.addNewPoint(15);
  EXPECT_EQ(rebuildPredictor.predict(), 12);

  // Testing reset function, e.g. when a new configuration is used.
  rebuildPredictor.reset();
  rebuildPredictor.addNewPoint(40);
  EXPECT_EQ(rebuildPredictor.predict(), 40);
}

TEST(DynamicRebuildDecisionWithFPBTest, testRemainderTraversalEstimateManager) {
  RemainderTraversalEstimateManager remTravPredictor{};
  EXPECT_EQ(remTravPredictor.predict(10), 0);

  // First sample
  remTravPredictor.addNewPoint(10, 150);
  EXPECT_EQ(remTravPredictor.predict(15), 0);

  // another sample (x,y) with same x value is ignored and one sample point is not enough for a linear regression, so
  // expected prediction is still zero.
  remTravPredictor.addNewPoint(10, 160);
  EXPECT_EQ(remTravPredictor.predict(15), 0);

  // Second sample
  remTravPredictor.addNewPoint(20, 300);
  EXPECT_EQ(remTravPredictor.predict(30), 450);
  EXPECT_EQ(remTravPredictor.getSumRemainderTraversalTime(), 610);

  // a new sample must have a larger x value than previously seen or it will be ignored.
  remTravPredictor.addNewPoint(18, 300);
  EXPECT_EQ(remTravPredictor.predict(30), 450);

  // Third sample
  remTravPredictor.addNewPoint(30, 450);
  EXPECT_EQ(remTravPredictor.predict(40), 600);

  // Fourth sample
  remTravPredictor.addNewPoint(60, 900);
  EXPECT_EQ(remTravPredictor.predict(80), 1200);

  // Fifth sample
  remTravPredictor.addNewPoint(90, 1350);
  EXPECT_EQ(remTravPredictor.predict(100), 1500);

  // Sample 6-10, to show that only last five samples are used for prediction
  remTravPredictor.addNewPoint(150, 3000);
  remTravPredictor.addNewPoint(300, 6000);
  remTravPredictor.addNewPoint(600, 12000);
  remTravPredictor.addNewPoint(1200, 24000);
  remTravPredictor.addNewPoint(2400, 48000);

  EXPECT_EQ(remTravPredictor.predict(4000), 80000);

  // but sum contains all the samples added so far
  EXPECT_EQ(remTravPredictor.getSumRemainderTraversalTime(), 96610);

  // Testing reset
  remTravPredictor.reset();
  EXPECT_EQ(remTravPredictor.predict(10), 0);
  EXPECT_EQ(remTravPredictor.getSumRemainderTraversalTime(), 0);
}

TEST(DynamicRebuildDecisionWithFPBTest, testRebuildDecisionContext) {
  RebuildDecisionContext rebuildDecisionContext{};

  // Steps since rebuild: 1
  unsigned int rf = 1;
  // First rebuild
  EXPECT_TRUE(std::isnan(rebuildDecisionContext.getRebuildNeighborTimeEstimate()));
  rebuildDecisionContext.afterRebuild(100);
  EXPECT_TRUE(std::isnan(rebuildDecisionContext.getRebuildNeighborTimeEstimate()));

  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 0.);
  rebuildDecisionContext.afterRemainderTraversal(0, 0);
  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 0.);

  rebuildDecisionContext.updateNumParticlesBufferEstimate(3);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 3);
  // Not enough remainder traversal times collected
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(rf));

  // Rebuild estimate is now updated.
  EXPECT_FALSE(std::isnan(rebuildDecisionContext.getRebuildNeighborTimeEstimate()));
  EXPECT_EQ(rebuildDecisionContext.getRebuildNeighborTimeEstimate(), 100);

  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 0.);

  // Steps since rebuild: 2
  rf++;
  rebuildDecisionContext.afterRemainderTraversal(20, 4);
  rebuildDecisionContext.updateNumParticlesBufferEstimate(4);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 5);
  // Not enough remainder traversal times collected
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(rf));

  EXPECT_EQ(rebuildDecisionContext.getRebuildNeighborTimeEstimate(), 100);
  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 25);

  // Steps since rebuild: 3
  rf++;
  // second real sample (as we ignore the samples where number of buffer particles is zero)
  rebuildDecisionContext.afterRemainderTraversal(25, 5);

  rebuildDecisionContext.updateNumParticlesBufferEstimate(6);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 7);

  // Another remainder traversal with a full buffer is cheaper than a rebuild
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(rf));

  // Updated after the last decision
  EXPECT_EQ(rebuildDecisionContext.getRebuildNeighborTimeEstimate(), 100);
  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 35);

  // Steps since rebuild: 4
  rf++;
  rebuildDecisionContext.afterRemainderTraversal(35, 7);
  rebuildDecisionContext.updateNumParticlesBufferEstimate(20);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 21);
  // Rebuild is cheaper than another remainder traversal with a full buffer
  EXPECT_TRUE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(rf));
  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 105);

  // Rebuild
  rebuildDecisionContext.afterRebuild(120);

  // Steps since rebuild: 1
  // Check reset behavior and updated estimate of rebuild cost
  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 0);

  rebuildDecisionContext.updateNumParticlesBufferEstimate(1);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 2);
  // Not enough remainder traversal times collected
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(1));
  EXPECT_EQ(rebuildDecisionContext.getRebuildNeighborTimeEstimate(), 110);
  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 0);

  // Reset internal rebuild time estimator when configuration changes.
  rebuildDecisionContext.afterConfigurationChange();
  rebuildDecisionContext.afterRebuild(10);
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(1));
  EXPECT_EQ(rebuildDecisionContext.getRebuildNeighborTimeEstimate(), 10);
}
