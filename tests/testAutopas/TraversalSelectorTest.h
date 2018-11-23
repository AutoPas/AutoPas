/**
 * @file TraversalSelectorTest.h
 * @author F. Gratl
 * @date 21.06.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/particles/Particle.h"
#include "autopas/selectors/TraversalSelector.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class TraversalSelectorTest : public AutoPasTestBase {
 public:
  TraversalSelectorTest() = default;
  ~TraversalSelectorTest() override = default;

  void testFastest(autopas::SelectorStrategy strategy,
                   std::unordered_map<autopas::TraversalOptions, std::vector<long>> measurements,
                   autopas::TraversalOptions expectedBest,
                   std::unordered_map<autopas::TraversalOptions, std::vector<long>> ignoredMeasurements = {});
};
