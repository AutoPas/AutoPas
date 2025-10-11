/**
 * @file ATFunctorTestGlobals.h
 * @author muehlhaeusser
 * @date 29.08.23
 */

#pragma once

#include <gtest/gtest.h>

#include "ATFunctorTest.h"
#include "testingHelpers/commonTypedefs.h"

template <class FuncType>
class ATFunctorTestGlobals : public ATFunctorTest {
 public:
  struct OwnershipConfig {
    double factor{};
    std::string where_str;
    bool owned1{}, owned2{}, owned3{};
  };

  ATFunctorTestGlobals() : ATFunctorTest() {}

  static void ATFunctorTestGlobalsNoMixingAoS(where_type where, bool newton3);

  void runATFunctorGlobalsTest(where_type where, SoAFunctorType soaFunctorType, bool newton3);

  constexpr static double cutoff{5.};
  constexpr static double nu{0.7};

  constexpr static double absDelta{1e-8};
};
