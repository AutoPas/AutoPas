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

namespace AutoTunerTest {

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
};

}  // end namespace AutoTunerTest
