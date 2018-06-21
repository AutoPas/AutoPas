/**
 * C08TraversalTest.h
 *
 *  Created on: 5/24/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPas.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/GridGenerator.h"
#include <testingHelpers/commonTypedefs.h>

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

class C08TraversalTest : public testing::Test {
 public:
  C08TraversalTest() = default;

  ~C08TraversalTest() = default;
};
