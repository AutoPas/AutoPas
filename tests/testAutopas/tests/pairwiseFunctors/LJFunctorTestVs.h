/**
 * @file LJFunctorTestVs.cpp
 * @author F. Gratl
 * @date 20.03.20
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "LJFunctorTest.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopasTools/generators/RandomGenerator.h"

template <class FuncType>
class LJFunctorTestVs : public LJFunctorTest {
 public:
  LJFunctorTestVs() : LJFunctorTest() {}

  enum InteractionType { own, pair, verlet };

  enum where_type { inside, boundary, outside };

  constexpr static double cutoff{1.};
  constexpr static double epsilon{1.};
  constexpr static double sigma{1.};

  constexpr static std::array<double, 3> expectedForce{-4547248.8989645941, -9094497.7979291882, -13641746.696893783};
  constexpr static std::array<double, 3> expectedForceMixing{-835415983.7676939964294, -1670831967.5353879928588,
                                                             -2506247951.3030819892883};

  constexpr static double expectedVirial{6366148.4585504318};
  constexpr static double expectedEnergy{529783.50857210846};
  constexpr static double absDelta{1e-7};
};
