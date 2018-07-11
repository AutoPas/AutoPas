/**
 * @file C08TraversalTest.h
 * @author F. Gratl
 * @date 24.05.18
 */

#pragma once

#include <gtest/gtest.h>
#include <testingHelpers/commonTypedefs.h>
#include "autopas/AutoPas.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/GridGenerator.h"

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

class C08TraversalTest : public testing::Test {
 public:
  C08TraversalTest() = default;

  ~C08TraversalTest() = default;
};
