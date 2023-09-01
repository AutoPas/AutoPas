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

/**
 * Helper struct to get access to template parameters
 */
template <autopas::DataLayoutOption::Value dataLayout, bool useNewton3, bool bidirectional>
struct CellFunctorWrapper {
  using CellFT = autopas::internal::CellFunctor<Molecule, autopas::FullParticleCell<Molecule>,
                                                mdLib::LJFunctor<Molecule>, dataLayout, useNewton3, bidirectional>;

  static const autopas::DataLayoutOption::Value dataLayoutV = dataLayout;
  static const bool useNewton3V = useNewton3;
  static const bool bidirectionalV = bidirectional;
};

template <typename T>
class CellFunctorTest : public AutoPasTestBase {
 public:
  CellFunctorTest() = default;

  ~CellFunctorTest() override = default;

 private:
};