/**
 * @file ATMFunctorTestGlobals.h
 * @author muehlhaeusser
 * @date 29.08.23
 */

#pragma once

#include <gtest/gtest.h>

#include "ATMFunctorTest.h"
#include "AutoPasTestBase.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "testingHelpers/ATMPotential.h"
#include "testingHelpers/commonTypedefs.h"

template <class FuncType>
class ATMFunctorTestGlobals : public ATMFunctorTest {
 public:
  ATMFunctorTestGlobals() : ATMFunctorTest() {}

  static void ATMFunctorTestGlobalsNoMixing(where_type where, bool newton3);
  static void testSoAGlobalsATM(where_type where, bool newton3, InteractionType interactionType,
                                size_t additionalParticlesToVerletNumber, uint64_t numParticleReplicas,
                                bool mixedNewton3FunctorCalls);

  constexpr static double cutoff{5.};
  constexpr static double nu{0.7};

  constexpr static double absDelta{1e-8};
};
