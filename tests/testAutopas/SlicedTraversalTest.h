/**
 * @file SlicedTraversalTest.h
 * @author F. Gratl
 * @date 01.05.18
 */

#pragma once

#include <gtest/gtest.h>
#include "autopas/AutoPas.h"
#include "mocks/MockFunctor.h"
#include "testingHelpers/commonTypedefs.h"

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

class SlicedTraversalTest : public testing::Test {
 public:
  SlicedTraversalTest() = default;

  ~SlicedTraversalTest() override = default;

  void fillWithParticles(std::vector<FPCell> &cells, std::array<size_t, 3> particlesPerDim);
};
