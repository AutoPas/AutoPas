/**
 * @file ATFunctorTestGlobals.h
 * @author muehlhaeusser
 * @date 29.08.23
 */

#pragma once

#include <gtest/gtest.h>

#include "AutoPasTestBase.h"
#include "autopasTools/generators/RandomGenerator.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "testingHelpers/ATPotential.h"
#include "testingHelpers/commonTypedefs.h"

using ATFunctorGlobalsTestingTuple = std::tuple<bool /*mixing*/, bool /*newton3*/>;

class ATFunctorTestGlobals : public AutoPasTestBase, public ::testing::WithParamInterface<ATFunctorGlobalsTestingTuple> {

 public:

  enum InteractionType { own, pair, verlet };
  enum where_type { inside, boundary, outside };

  ATFunctorTestGlobals() : AutoPasTestBase() {}



  static void ATFunctorTestGlobalsNoMixing(bool newton3);
//  static void testSoAGlobals(where_type where, bool newton3, InteractionType interactionType,
//                             size_t additionalParticlesToVerletNumber, uint64_t numParticleReplicas,
//                             bool mixedNewton3FunctorCalls);
//  static void testAoSGlobalsMixedN3(ATFunctorTestGlobals<FuncType>::where_type where);

  constexpr static double cutoff{5.};
  constexpr static double nu{0.7};

  constexpr static double absDelta{1e-10};
};
