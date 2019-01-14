/**
 * @file C01TraversalTest.h
 * @author S. Seckler
 * @date 10.01.2019
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/GridGenerator.h"
#include "testingHelpers/commonTypedefs.h"

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

class C01TraversalTest : public AutoPasTestBase {
 public:
  C01TraversalTest() = default;

  ~C01TraversalTest() override = default;
};
