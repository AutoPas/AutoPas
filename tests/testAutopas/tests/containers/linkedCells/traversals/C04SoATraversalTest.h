/**
 * @file C04SoATraversalTest.h
 * @author C.Menges
 * @date 06.07.19
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

class C04SoATraversalTest : public AutoPasTestBase {
 public:
  C04SoATraversalTest() = default;

  ~C04SoATraversalTest() override = default;
};
