/**
 * @file LJFunctorTestGlobals.h
 * @author F. Gratl
 * @date 20.03.20
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopasTools/generators/RandomGenerator.h"

template <class FuncType>
class LJFunctorTestGlobals : public AutoPasTestBase {
 public:
  LJFunctorTestGlobals() : AutoPasTestBase() {}

  /**
   * Checks if the given function throws an exception containing "not implemented".
   * @tparam FunType Type of the given function.
   * @param f Code to be checked as a lambda.
   * @return Empty string if nothing was caught, the exception string if a matching exception was found.
   * If the exception does not match it is rethrown.
   */
  template <class FunType>
  static std::string shouldSkipIfNotImplemented(FunType &&f);

  enum InteractionType { own, pair, verlet };
  enum where_type { inside, boundary, outside };

  static void testAoSGlobals(where_type where, bool newton3, bool duplicatedCalculation);
  static void testSoAGlobals(where_type where, bool newton3, bool duplicatedCalculation,
                             InteractionType interactionType, size_t additionalParticlesToVerletNumber,
                             bool cellWiseOwnedState, uint64_t numParticleReplicas);

  constexpr static double cutoff{1.};
  constexpr static double epsilon{1.};
  constexpr static double sigma{1.};

  constexpr static double expectedVirial{6366148.4585504318};
  constexpr static double expectedEnergy{529783.50857210846};
  constexpr static double absDelta{1e-7};
};
