/**
 * @file CellFunctorTest.h
 * @author D. Martin
 * @date 29.08.23
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/LJPotential.h"
#include "testingHelpers/commonTypedefs.h"

template <typename T>
class CellFunctorTest : public AutoPasTestBase {
 public:
  CellFunctorTest() = default;

  ~CellFunctorTest() override = default;
};