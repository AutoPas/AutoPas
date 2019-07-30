/**
 * @file C01TraversalAdaptiveTest.h
 * @author C. Menges
 * @date 11.07.2019
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/GridGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class C01TraversalAdaptiveTest : public AutoPasTestBase {
 public:
  C01TraversalAdaptiveTest() = default;

  ~C01TraversalAdaptiveTest() override = default;
};
