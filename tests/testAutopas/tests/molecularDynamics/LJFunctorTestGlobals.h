/**
 * @file LJFunctorTestGlobals.h
 * @author F. Gratl
 * @date 20.03.20
 */

#pragma once

#include <gtest/gtest.h>

#include "LJFunctorTest.h"
#include "autopas/molecularDynamics/ParticlePropertiesLibrary.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "testingHelpers/LJPotential.h"

template <class FuncType>
class LJFunctorTestGlobals : public LJFunctorTest {
 public:
  LJFunctorTestGlobals() : LJFunctorTest() {}

  static void testAoSGlobals(where_type where, bool newton3);
  static void testSoAGlobals(where_type where, bool newton3, InteractionType interactionType,
                             size_t additionalParticlesToVerletNumber, uint64_t numParticleReplicas,
                             bool mixedNewton3FunctorCalls);
  static void testAoSGlobalsMixedN3(LJFunctorTestGlobals<FuncType>::where_type where);

  constexpr static double cutoff{1.};
  constexpr static double epsilon{1.};
  constexpr static double sigma{1.};

  constexpr static double absDelta{1e-7};
};
