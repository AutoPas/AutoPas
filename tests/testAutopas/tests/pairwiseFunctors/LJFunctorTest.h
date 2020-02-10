/**
 * @file LJFunctorTest.h
 * @author seckler
 * @date 06.11.18
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopasTools/generators/RandomGenerator.h"

class LJFunctorTest : public AutoPasTestBase {
 public:
  LJFunctorTest() : AutoPasTestBase() {}

  void SetUp() override{};

  void TearDown() override{};

 protected:
  template <bool Mixing>
  void testAoSNoGlobals(bool newton3);

  enum InteractionType { own, pair, verlet };

  template <bool Mixing>
  void testSoANoGlobals(bool newton3, InteractionType interactionType);

  enum where_type { inside, boundary, outside };
  void testAoSGlobals(where_type where, bool newton3, bool duplicatedCalculation);
  void testSoAGlobals(where_type where, bool newton3, bool duplicatedCalculation, InteractionType interactionType,
                      size_t additionalParticlesToVerletNumber);

  constexpr static double cutoff{1.};
  constexpr static double epsilon{1.};
  constexpr static double sigma{1.};
  constexpr static double epsilon2{2.};
  constexpr static double sigma2{2.};

  constexpr static std::array<double, 3> expectedForce{-4547248.8989645941, -9094497.7979291882, -13641746.696893783};
  constexpr static std::array<double, 3> expectedForceMixing{-835415983.7676939964294, -1670831967.5353879928588,
                                                             -2506247951.3030819892883};

  constexpr static double expectedVirial{6366148.4585504318};
  constexpr static double expectedEnergy{529783.50857210846};
  constexpr static double absDelta{1e-7};
};
