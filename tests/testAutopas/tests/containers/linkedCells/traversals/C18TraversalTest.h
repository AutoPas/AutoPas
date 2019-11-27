/**
 * @file C18TraversalTest.h
 * @author S. Seckler
 * @date 10.01.2019
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas-tools/generators/GridGenerator.h"
#include "autopas/AutoPas.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class C18TraversalTest : public AutoPasTestBase {
 public:
  C18TraversalTest() = default;

  ~C18TraversalTest() override = default;
};
