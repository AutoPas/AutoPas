/**
 * @file TraversalSelectorTest.h
 * @author F. Gratl
 * @date 21.06.18
 */

#pragma once

#include <gtest/gtest.h>
#include <mocks/MockFunctor.h>
#include <particles/Particle.h>
#include <selectors/TraversalSelector.h>
#include <testingHelpers/commonTypedefs.h>
#include "AutoPasTestBase.h"

class TraversalSelectorTest : public AutoPasTestBase {
 public:
  TraversalSelectorTest();
  ~TraversalSelectorTest();
};
