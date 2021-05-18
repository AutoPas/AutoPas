/**
 * @file AutoTunerTest.h
 * @author F. Gratl
 * @date 8/10/18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopas/selectors/AutoTuner.h"
#include "testingHelpers/commonTypedefs.h"

class AutoTunerTest : public AutoPasTestBase {
 public:
  AutoTunerTest() = default;
  ~AutoTunerTest() = default;

  /**
   * Map multiple runtimes to one configuration each.
   */
  using mapConfigTime = std::map<autopas::Configuration, std::vector<long>>;

  /**
   * Create a tuner with
   * @param strategy
   * @param configAndTimes
   * @param expectedBest
   * @param ignoredConfigAndTimes
   */
  void testFastest(autopas::SelectorStrategyOption strategy, mapConfigTime configAndTimes,
                   autopas::Configuration expectedBest, mapConfigTime ignoredConfigAndTimes = {});

  const double _cellSizeFactor{1.};

  const autopas::Configuration _confLc_c01{autopas::ContainerOption::linkedCells, _cellSizeFactor,
                                    autopas::TraversalOption::lc_c01, autopas::LoadEstimatorOption::none,
                                    autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled};
  const autopas::Configuration _confLc_c04{autopas::ContainerOption::linkedCells, _cellSizeFactor,
                                    autopas::TraversalOption::lc_c04, autopas::LoadEstimatorOption::none,
                                    autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled};
  const autopas::Configuration _confLc_c08{autopas::ContainerOption::linkedCells, _cellSizeFactor,
                                    autopas::TraversalOption::lc_c08, autopas::LoadEstimatorOption::none,
                                    autopas::DataLayoutOption::aos, autopas::Newton3Option::disabled};
};
