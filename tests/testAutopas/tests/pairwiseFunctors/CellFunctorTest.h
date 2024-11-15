/**
 * @file CellFunctorTest.h
 * @author D. Martin
 * @date 29.08.23
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/baseFunctors/Functor.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/commonTypedefs.h"

template <typename T>
class CellFunctorTest : public AutoPasTestBase {
 public:
  CellFunctorTest() = default;

  ~CellFunctorTest() override = default;
};