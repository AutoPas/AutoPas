/**
 * @file MPIParallelizedStrategyTest.h
 * @author W. Thieme
 * @date 11.06.2020
 */

#pragma once

#include <gtest/gtest.h>

#include <array>

#include "AutoPasMPITestBase.h"
#include "autopas/tuning/Configuration.h"

class MPIParallelizedStrategyTest : public AutoPasMPITestBase {
 protected:
  static constexpr int numRanksExpected = 3;

  /**
   * Test function that does the following steps:
   *   - set up a AutoTuner with the MPI divide and conquer strategy
   *   - trigger the strategy to resort / filter the queues
   *   - check that all queues have the correct size and add up correctly
   * @param particleDependentBinDensityStdDevs
   * @param expectedNumLocalConfigs
   * @param searchSpace
   */
  void testBucketDistribution(const std::array<double, numRanksExpected> &particleDependentBinDensityStdDevs,
                              const std::array<size_t, numRanksExpected> &expectedNumLocalConfigs,
                              const std::set<autopas::Configuration> &searchSpace);

  autopas::Configuration lc_c01_aos{
      autopas::ContainerOption::linkedCells,    1.,
      autopas::TraversalOption::lc_c01,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  autopas::Configuration lc_c04_aos{
      autopas::ContainerOption::linkedCells,    1.,
      autopas::TraversalOption::lc_c04,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  autopas::Configuration lc_c08_aos{
      autopas::ContainerOption::linkedCells,    1.,
      autopas::TraversalOption::lc_c08,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::aos,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  autopas::Configuration lc_c01_soa{
      autopas::ContainerOption::linkedCells,    1.,
      autopas::TraversalOption::lc_c01,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  autopas::Configuration lc_c04_soa{
      autopas::ContainerOption::linkedCells,    1.,
      autopas::TraversalOption::lc_c04,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p1xVec};
  autopas::Configuration lc_c01_soa_p2xVecDiv2{
      autopas::ContainerOption::linkedCells,    1.,
      autopas::TraversalOption::lc_c01,         autopas::LoadEstimatorOption::none,
      autopas::DataLayoutOption::soa,           autopas::Newton3Option::disabled,
      autopas::InteractionTypeOption::pairwise, autopas::VectorizationPatternOption::p2xVecDiv2};
};