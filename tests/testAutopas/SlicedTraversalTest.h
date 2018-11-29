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
#include "AutoPasTestBase.h"

#ifdef AUTOPAS_OPENMP
#include <omp.h>
#endif

class SlicedTraversalTest : public AutoPasTestBase {
 public:
  SlicedTraversalTest() = default;

  ~SlicedTraversalTest() override = default;

  void fillWithParticles(std::vector<FPCell> &cells, std::array<size_t, 3> particlesPerDim);
};
