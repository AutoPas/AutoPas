/**
 * @file TraversalSelectorTest.h
 * @author F. Gratl
 * @date 21.06.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/particles/ParticleDefinitions.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "mocks/MockPairwiseFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class TraversalSelectorTest : public AutoPasTestBase {
 public:
  TraversalSelectorTest() = default;
  ~TraversalSelectorTest() override = default;
};
