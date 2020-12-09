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

namespace TraversalSelectorTest {

class TraversalSelectorTest : public AutoPasTestBase {
 public:
  TraversalSelectorTest() = default;
  ~TraversalSelectorTest() override = default;
};

}  // end namespace TraversalSelectorTest
