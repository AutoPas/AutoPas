/**
 * @file LJFunctorTestGlobals.cpp
 * @author F. Gratl
 * @date 20.03.20
 */

#pragma once

#include <gtest/gtest.h>

#include "LJFunctorTest.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopasTools/generators/RandomGenerator.h"

template <class FuncType>
class LJFunctorTestNoGlobals : public LJFunctorTest {
 public:
  LJFunctorTestNoGlobals() : LJFunctorTest() {}

  enum InteractionType { own, pair, verlet };
  enum where_type { inside, boundary, outside };

  template <bool mixing>
  static void testSoANoGlobals(bool newton3, InteractionType interactionType);
  template <bool mixing>
  static std::string testAoSNoGlobals(bool newton3);

  constexpr static double cutoff{1.};
  constexpr static double epsilon{1.};
  constexpr static double sigma{1.};
  constexpr static double epsilon2{2.};
  constexpr static double sigma2{2.};

  constexpr static std::array<double, 3> expectedForce{-4547248.8989645941, -9094497.7979291882, -13641746.696893783};
  constexpr static std::array<double, 3> expectedForceMixing{-835415983.7676939964294, -1670831967.5353879928588,
                                                             -2506247951.3030819892883};

  constexpr static double absDelta{1e-7};
};
