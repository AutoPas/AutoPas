/**
 * @file CellFunctorTest3B.h
 * @author muehlhaeusser
 * @date 05.12.23
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/baseFunctors/CellFunctor3B.h"
#include "autopas/baseFunctors/Functor.h"
#include "testingHelpers/commonTypedefs.h"

template <typename T>
class CellFunctorTest3B : public AutoPasTestBase {
 public:
  CellFunctorTest3B() = default;

  ~CellFunctorTest3B() override = default;

  static constexpr double cutoff = 1.0;
  static constexpr double nu = 1.0;
};
