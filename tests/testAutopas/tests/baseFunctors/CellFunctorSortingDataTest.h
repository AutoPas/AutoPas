/**
 * @file CellFunctorSortingDataTest.h
 * @author hmeyran
 * @date 30.06.2026
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/baseFunctors/CellFunctor.h"
#include "molecularDynamicsLibrary/LJFunctor.h"
#include "testingHelpers/commonTypedefs.h"

class CellFunctorSortingDataTest : public AutoPasTestBase {
 protected:
  static constexpr double cutoff = 3.0;
  // Functor declared before cf so it is constructed first (member init order = declaration order).
  LJFunctorType<> _functor{cutoff};
  autopas::internal::CellFunctor<FMCell, LJFunctorType<>, /*bidirectional=*/false> _cf{
      _functor, cutoff, autopas::DataLayoutOption::soa, /*useNewton3=*/false};
};
