/**
 * @file RegressionModelsTest.cpp
 * @author Natallia Padalinskaya
 * @date 28.12.25
 */

#include "RegressionModelsTest.h"

using namespace autopas::utils;

TEST(RegressionModelsTest, MeanTest) {
  Mean mean(1, 3);
  EXPECT_EQ(mean.getN(), 0u);

  // Add points
  EXPECT_EQ(mean.addNewPoint(10), RegressionBase::ReturnCode::OK);
  EXPECT_EQ(mean.addNewPoint(20), RegressionBase::ReturnCode::OK);

  // check predict
  auto result = mean.predict();
  EXPECT_TRUE(result._isOk);
  EXPECT_DOUBLE_EQ(result._value, 15.0);

  EXPECT_EQ(mean.addNewPoint(30), RegressionBase::ReturnCode::OK);

  result = mean.predict();
  EXPECT_TRUE(result._isOk);
  EXPECT_DOUBLE_EQ(result._value, 20.0);

  // exceed max amount of points
  EXPECT_EQ(mean.addNewPoint(40), RegressionBase::ReturnCode::EXCEEDED_MAX_POINTS);
  result = mean.predict();
  EXPECT_TRUE(result._isOk);
  EXPECT_DOUBLE_EQ(result._value, 20.0);

  // Reset
  mean.reset();
  EXPECT_EQ(mean.getN(), 0u);

  // correct behavior if there's not enough points for predict
  EXPECT_EQ(mean.predict()._returnCode, RegressionBase::ReturnCode::NOT_ENOUGH_POINTS);

  EXPECT_EQ(mean.addNewPoint(20), RegressionBase::ReturnCode::OK);
  result = mean.predict();
  EXPECT_TRUE(result._isOk);
  EXPECT_DOUBLE_EQ(result._value, 20.0);
}

TEST(RegressionModelsTest, SimpleLinearRegressionTest) {
  SimpleLinearRegression reg(2, 3);

  EXPECT_EQ(reg.addNewPoint(1, 2), RegressionBase::ReturnCode::OK);
  EXPECT_EQ(reg.addNewPoint(2, 3), RegressionBase::ReturnCode::OK);

  auto result = reg.predict(3);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 4.0, 1e-6);  // y = x + 1

  EXPECT_EQ(reg.addNewPoint(3, 7), RegressionBase::ReturnCode::OK);

  result = reg.predict(6);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 14.0, 1e-6);  // y = 2.5*x - 1

  EXPECT_EQ(reg.addNewPoint(6, 21), RegressionBase::ReturnCode::EXCEEDED_MAX_POINTS);  // would be 4*x - 3.75

  result = reg.predict(6);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 14.0, 1e-6);  // y = 2.5*x - 1

  // (point not relevant for prediction but for the sum)
  EXPECT_EQ(reg.getSumY(), 33);

  // Reset
  reg.reset();
  EXPECT_EQ(reg.getN(), 0u);
  EXPECT_EQ(reg.getSumY(), 0);

  EXPECT_EQ(reg.predict(6)._returnCode, RegressionBase::ReturnCode::NOT_ENOUGH_POINTS);
}

TEST(RegressionModelsTest, SimpleLinearRegressionBoostTest) {
  SimpleLinearRegressionBoost reg(2, 3);

  EXPECT_EQ(reg.addNewPoint(1, 2), RegressionBase::ReturnCode::OK);
  EXPECT_EQ(reg.addNewPoint(2, 3), RegressionBase::ReturnCode::OK);

  auto result = reg.predict(3);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 4.0, 1e-6);  // y = x + 1

  EXPECT_EQ(reg.addNewPoint(3, 7), RegressionBase::ReturnCode::OK);

  result = reg.predict(6);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 14.0, 1e-6);  // y = 2.5*x - 1

  // ring buffer
  // would be with all y = 6.357*x - 8.071 but without last y = 7.5*x - 13.5
  EXPECT_EQ(reg.addNewPoint(6, 32), RegressionBase::ReturnCode::OK);

  result = reg.predict(7);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 39.0, 1e-6);  // y = 7.5*x - 13.5

  // (point not relevant for prediction but for the sum)
  EXPECT_EQ(reg.getSumY(), 44);

  // Reset testen
  reg.reset();
  EXPECT_EQ(reg.getN(), 0u);
  EXPECT_EQ(reg.getSumY(), 0u);

  EXPECT_EQ(reg.predict(6)._returnCode, RegressionBase::ReturnCode::NOT_ENOUGH_POINTS);

  // numDifferentXConsidered
  EXPECT_EQ(reg.getNumDifferentXConsidered(), 0u);
  EXPECT_EQ(reg.addNewPoint(1, 3), RegressionBase::ReturnCode::OK);
  EXPECT_EQ(reg.getNumDifferentXConsidered(), 1u);
  EXPECT_EQ(reg.addNewPoint(1, 2), RegressionBase::ReturnCode::OK);
  EXPECT_EQ(reg.getNumDifferentXConsidered(), 1u);

  EXPECT_EQ(reg.addNewPoint(2, 3), RegressionBase::ReturnCode::OK);
  EXPECT_EQ(reg.addNewPoint(3, 7), RegressionBase::ReturnCode::OK);

  // 4 points are considered even though _maxN = 3, because duplicate x values are not removed.
  EXPECT_EQ(reg.addNewPoint(6, 21), RegressionBase::ReturnCode::OK);  // 4*x - 3.75

  EXPECT_EQ(reg.getNumDifferentXConsidered(), 3u);

  result = reg.predict(7);
  EXPECT_TRUE(result._isOk);
  EXPECT_NEAR(result._value, 24.25, 1e-6);  // 4*x - 3.75
}

TEST(RegressionModelsTest, RebuildDecisionContext) {
  RebuildDecisionContext rebuildDecisionContext{};

  // steps since rebuild 1
  // first rebuild
  EXPECT_TRUE(std::isnan(rebuildDecisionContext.getRebuildNeighborTimeEstimate()));
  rebuildDecisionContext.afterRebuild(100, false, false);

  EXPECT_FALSE(rebuildDecisionContext.afterRemainderTraversal(10, 0));
  rebuildDecisionContext.updateNumParticlesBufferEstimate(3);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 3u);
  // not enough remainder Traversal times
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(1));

  // steps since rebuild 2
  EXPECT_FALSE(rebuildDecisionContext.afterRemainderTraversal(20, 4));
  rebuildDecisionContext.updateNumParticlesBufferEstimate(4);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 5u);
  // not enough remainder Traversal times
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(2));

  // steps since rebuild 3
  EXPECT_FALSE(rebuildDecisionContext.afterRemainderTraversal(30, 5));
  EXPECT_TRUE(std::isnan(rebuildDecisionContext.getRemainderTraversalTimeEstimate()));
  rebuildDecisionContext.updateNumParticlesBufferEstimate(6);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 7u);
  // doesn't fulfill the criterion
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(3));
  // updated after last decision
  EXPECT_EQ(rebuildDecisionContext.getRebuildNeighborTimeEstimate(), 100);
  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 50);

  // steps since rebuild 4
  EXPECT_FALSE(rebuildDecisionContext.afterRemainderTraversal(50, 7));
  rebuildDecisionContext.updateNumParticlesBufferEstimate(20);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 21u);
  // fulfill the criterion
  EXPECT_TRUE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(4));
  EXPECT_EQ(rebuildDecisionContext.getRemainderTraversalTimeEstimate(), 190);
  rebuildDecisionContext.afterRebuild(120, false, false);

  // steps since rebuild 1
  // check reset
  EXPECT_TRUE(std::isnan(rebuildDecisionContext.getRemainderTraversalTimeEstimate()));

  rebuildDecisionContext.updateNumParticlesBufferEstimate(1);
  EXPECT_EQ(rebuildDecisionContext.getNumParticlesBufferEstimate(), 2u);
  // not enough remainder Traversal times
  EXPECT_FALSE(rebuildDecisionContext.decideToRebuildOnParticleBufferFullness(1));
}