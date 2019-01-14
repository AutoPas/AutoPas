/**
 * @file C18TraversalTest.h
 * @author S. Seckler
 * @date 10.01.2019
 */

#pragma once

#include <gtest/gtest.h>
#include "testingHelpers/commonTypedefs.h"
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/GridGenerator.h"

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

class C18TraversalTest : public AutoPasTestBase {
 public:
  C18TraversalTest() = default;

  ~C18TraversalTest() override = default;
};
