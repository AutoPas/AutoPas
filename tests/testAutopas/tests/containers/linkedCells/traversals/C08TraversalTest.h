/**
 * @file C08TraversalTest.h
 * @author F. Gratl
 * @date 24.05.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/GridGenerator.h"
#include "testingHelpers/commonTypedefs.h"

class C08TraversalTest : public AutoPasTestBase {
 public:
  C08TraversalTest() = default;

  ~C08TraversalTest() override = default;
};
