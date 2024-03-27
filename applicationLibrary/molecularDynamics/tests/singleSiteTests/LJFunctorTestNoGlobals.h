/**
 * @file LJFunctorTestNoGlobals.h
 * @author F. Gratl
 * @date 20.03.20
 */

#pragma once

#include <gtest/gtest.h>

#include "LJFunctorTest.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

template <class FuncType>
class LJFunctorTestNoGlobals : public LJFunctorTest {
 public:
  LJFunctorTestNoGlobals() : LJFunctorTest() {}

  constexpr static double cutoff{1.};
  constexpr static double epsilon{1.};
  constexpr static double sigma{1.};
  constexpr static double epsilon2{2.};
  constexpr static double sigma2{2.};

  // if we make this static constexpr icpc sets them to 0
  const std::array<double, 3> expectedForce{-4547248.8989645941, -9094497.7979291882, -13641746.696893783};
  const std::array<double, 3> expectedForceMixing{-835415983.7676939964294, -1670831967.5353879928588,
                                                  -2506247951.3030819892883};

  constexpr static double absDelta{1e-7};
};
