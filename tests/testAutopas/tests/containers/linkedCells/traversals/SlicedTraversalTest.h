/**
 * @file SlicedTraversalTest.h
 * @author F. Gratl
 * @date 01.05.18
 */

#pragma once

#include <gtest/gtest.h>
#include "AutoPasTestBase.h"
#include "autopas/AutoPas.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

class SlicedTraversalTest : public AutoPasTestBase {
 public:
  SlicedTraversalTest() = default;

  ~SlicedTraversalTest() override = default;
};
