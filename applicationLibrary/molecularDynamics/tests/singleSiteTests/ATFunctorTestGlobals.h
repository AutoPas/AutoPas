/**
 * @file ATFunctorTestGlobals.h
 * @author muehlhaeusser
 * @date 29.08.23
 */

#pragma once

#include <gtest/gtest.h>

#include "ATFunctorTest.h"
#include "AutoPasTestBase.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "testingHelpers/ATPotential.h"
#include "testingHelpers/commonTypedefs.h"

template <class FuncType>
class ATFunctorTestGlobals : public ATFunctorTest {
 public:
  ATFunctorTestGlobals() : ATFunctorTest() {}

  static void ATFunctorTestGlobalsNoMixing(where_type where, bool newton3);
  static void testSoAGlobalsAT(where_type where, bool newton3, InteractionType interactionType,
                               size_t additionalParticlesToVerletNumber, uint64_t numParticleReplicas,
                               bool mixedNewton3FunctorCalls);
  //  static void testAoSGlobalsMixedN3(ATFunctorTestGlobals<FuncType>::where_type where);

  constexpr static double cutoff{5.};
  constexpr static double nu{0.7};

  constexpr static double absDelta{1e-8};
};
