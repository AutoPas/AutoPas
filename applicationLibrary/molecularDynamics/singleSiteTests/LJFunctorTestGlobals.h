/**
 * @file LJFunctorTestGlobals.h
 * @author F. Gratl
 * @date 20.03.20
 */

#pragma once

#include <gtest/gtest.h>

#include "LJFunctorTest.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "autopasTools/generators/RandomGenerator.h"

template <class FuncType>
class LJFunctorTestGlobals : public LJFunctorTest {
 public:
  LJFunctorTestGlobals() : LJFunctorTest() {}

  static void testAoSGlobals(where_type where, bool newton3);
  static void testSoAGlobals(where_type where, bool newton3, InteractionType interactionType,
                             size_t additionalParticlesToVerletNumber, uint64_t numParticleReplicas);

  constexpr static double cutoff{1.};
  constexpr static double epsilon{1.};
  constexpr static double sigma{1.};

  constexpr static double expectedVirial{6366148.4585504318};
  constexpr static double expectedEnergy{529783.50857210846};
  constexpr static double absDelta{1e-7};
};
