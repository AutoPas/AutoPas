/**
 * @file RegressionModelsTest.cpp
 * @author Natallia Padalinskaya
 * @date 28.12.25
 */

#include "RegressionModelsTest.h"

#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS

using namespace autopas::utils::RegressionModels;

TEST(RegressionModelsTest, testMean) {
  RunningMean mean;
  Mean mean_2(2);  // Mean over last two samples
  Mean mean_3(3);  // Mean over last three samples
  EXPECT_EQ(mean.getN(), 0);
  EXPECT_EQ(mean_2.getN(), 0);
  EXPECT_EQ(mean_3.getN(), 0);

  // Add points
  EXPECT_EQ(mean.addNewPoint(10), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(mean.addNewPoint(20), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(mean_2.addNewPoint(10), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(mean_2.addNewPoint(20), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(mean_3.addNewPoint(10), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(mean_3.addNewPoint(20), RegressionBase::ReturnCode::OK_REG);

  // check predict
  auto result = mean.predict();
  EXPECT_TRUE(result._isOk);
  EXPECT_DOUBLE_EQ(result._value, 15);
  EXPECT_DOUBLE_EQ(mean_2.predict()._value, 15);
  EXPECT_DOUBLE_EQ(mean_3.predict()._value, 15);

  EXPECT_EQ(mean.addNewPoint(30), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(mean_2.addNewPoint(30), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(mean_3.addNewPoint(30), RegressionBase::ReturnCode::OK_REG);

  result = mean.predict();
  EXPECT_TRUE(result._isOk);
  EXPECT_DOUBLE_EQ(result._value, 20);
  EXPECT_DOUBLE_EQ(mean_2.predict()._value, 25);
  EXPECT_DOUBLE_EQ(mean_3.predict()._value, 20);

  // Exceed maximum number of points
  EXPECT_EQ(mean.addNewPoint(40), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(mean_2.addNewPoint(40), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(mean_3.addNewPoint(40), RegressionBase::ReturnCode::OK_REG);

  result = mean.predict();
  EXPECT_TRUE(result._isOk);
  EXPECT_DOUBLE_EQ(result._value, 25);
  EXPECT_DOUBLE_EQ(mean_2.predict()._value, 35);
  EXPECT_DOUBLE_EQ(mean_3.predict()._value, 30);

  // Reset
  mean.reset();
  mean_2.reset();
  mean_3.reset();

  EXPECT_EQ(mean.getN(), 0);

  // Correct behavior if there are not enough points for prediction
  EXPECT_EQ(mean.predict()._returnCode, RegressionBase::ReturnCode::NOT_ENOUGH_POINTS_REG);
  EXPECT_EQ(mean_2.predict()._returnCode, RegressionBase::ReturnCode::NOT_ENOUGH_POINTS_REG);
  EXPECT_EQ(mean_3.predict()._returnCode, RegressionBase::ReturnCode::NOT_ENOUGH_POINTS_REG);

  EXPECT_EQ(mean.addNewPoint(20), RegressionBase::ReturnCode::OK_REG);
  result = mean.predict();
  EXPECT_TRUE(result._isOk);
  EXPECT_DOUBLE_EQ(result._value, 20);
}

TEST(RegressionModelsTest, testSimpleLinearRegressionBoost) {
  SimpleLinearRegressionBoost reg(2, 3);

  EXPECT_EQ(reg.addNewPoint(1, 2), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(reg.addNewPoint(2, 3), RegressionBase::ReturnCode::OK_REG);

  auto result = reg.predict(3);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 4.0, 1e-6);  // y = x + 1 = (x=3) 4

  EXPECT_EQ(reg.addNewPoint(3, 7), RegressionBase::ReturnCode::OK_REG);

  result = reg.predict(6);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 14.0, 1e-6);  // y = 2.5*x - 1 = (x=6) 14

  EXPECT_EQ(reg.addNewPoint(6, 32), RegressionBase::ReturnCode::OK_REG);

  result = reg.predict(7);
  EXPECT_TRUE(result._isOk);
  // Ring buffer behavior
  // Using all points: y = 6.357*x - 8.071 = (x=7) 36.428, but without the oldest point: y = 7.5*x - 13.5 = (x=7) 39
  EXPECT_NEAR(result._value, 39.0, 1e-6);

  // The first added point (1,2) is no longer relevant for prediction but still contributes to the sum of y values
  EXPECT_EQ(reg.getSumY(), 44);

  // Reset
  reg.reset();
  EXPECT_EQ(reg.getN(), 0);
  EXPECT_EQ(reg.getSumY(), 0);

  EXPECT_EQ(reg.predict(6)._returnCode, RegressionBase::ReturnCode::NOT_ENOUGH_POINTS_REG);

  // numDifferentXConsidered
  EXPECT_EQ(reg.getNumDifferentXConsidered(), 0);
  EXPECT_EQ(reg.addNewPoint(1, 1), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(reg.getNumDifferentXConsidered(), 1);
  // adding another point with x<=1 would be ignored.
  EXPECT_EQ(reg.addNewPoint(1, 2), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(reg.getNumDifferentXConsidered(), 1);

  EXPECT_EQ(reg.addNewPoint(2, 2), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(reg.addNewPoint(3, 3), RegressionBase::ReturnCode::OK_REG);
  EXPECT_EQ(reg.getNumDifferentXConsidered(), 3);

  result = reg.predict(5);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 5, 1e-6);  // y = x

  // should replace (1,1)
  EXPECT_EQ(reg.addNewPoint(4, 5), RegressionBase::ReturnCode::OK_REG);
  // should replace (2,2)
  EXPECT_EQ(reg.addNewPoint(6, 9), RegressionBase::ReturnCode::OK_REG);

  EXPECT_EQ(reg.getNumDifferentXConsidered(), 3);

  result = reg.predict(7);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 11, 1e-6);  // y = 2x
}

TEST(RegressionModelsTest, testRebuildDecisionContext) {
  RebuildDecisionContext rebuildDecisionContext{};

  // Steps since rebuild: 1
  unsigned int rf = 1;
  // First rebuild
  EXPECT_TRUE(std::isnan(rebuildDecisionContext.getRebuildNeighborTimeEstimate()));
  rebuildDecisionContext.afterRebuild(100, false);

  EXPECT_FALSE(rebuildDecisionContext.afterRemainderTraversal(10, 0));
  rebuildDecisionContext.updateNumParticlesBufferEstimate(3);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 3);
  // Not enough remainder traversal times collected
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(rf));

  // Steps since rebuild: 2
  rf++;
  EXPECT_FALSE(rebuildDecisionContext.afterRemainderTraversal(20, 4));
  rebuildDecisionContext.updateNumParticlesBufferEstimate(4);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 5);
  // Not enough remainder traversal times collected
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(rf));

  // Steps since rebuild: 3
  rf++;
  EXPECT_FALSE(rebuildDecisionContext.afterRemainderTraversal(30, 5));
  EXPECT_TRUE(std::isnan(rebuildDecisionContext.getRemainderTraversalTimeEstimate()));
  rebuildDecisionContext.updateNumParticlesBufferEstimate(6);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 7);
  // Another remainder traversal with a full buffer is cheaper than a rebuild
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(rf));
  // Updated after the last decision
  EXPECT_EQ(rebuildDecisionContext.getRebuildNeighborTimeEstimate(), 100);
  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 50);

  // Steps since rebuild: 4
  rf++;
  EXPECT_FALSE(rebuildDecisionContext.afterRemainderTraversal(50, 7));
  rebuildDecisionContext.updateNumParticlesBufferEstimate(20);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 21);
  // Rebuild is cheaper than another remainder traversal with a full buffer
  EXPECT_TRUE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(rf));
  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 190);
  rebuildDecisionContext.afterRebuild(120, false);

  // Steps since rebuild: 1
  // Check reset behavior
  EXPECT_TRUE(std::isnan(rebuildDecisionContext.getRemainderTraversalTimeEstimate()));

  rebuildDecisionContext.updateNumParticlesBufferEstimate(1);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 2);
  // Not enough remainder traversal times collected
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(1));
}
#endif