/**
 * @file ATMFunctorFlopCounterTest.h
 * @author muehlhaeusser
 * @date 30.07.24
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/options/DataLayoutOption.h"

class ATMFunctorFlopCounterTest : public AutoPasTestBase, public ::testing::WithParamInterface<std::tuple<bool, bool>> {
 public:
  enum FunctorFunction { soaSingle, soaPair, soaTriple, verlet };
  static std::string to_string(FunctorFunction type) {
    switch (type) {
      case FunctorFunction::soaSingle:
        return "SoAFunctorSingle";
      case FunctorFunction::soaPair:
        return "SoAFunctorPair";
      case FunctorFunction::soaTriple:
        return "SoAFunctorTriple";
      case FunctorFunction::verlet:
        return "SoAFunctorVerlet";
      default:
        return "unknown";
    }
  }
  ATMFunctorFlopCounterTest() = default;

  ~ATMFunctorFlopCounterTest() override = default;

  template <bool calculateGlobals>
  void testFLOPCounterAoS(bool newton3);

  template <bool calculateGlobals>
  void testFLOPCounterSoA(bool newton3, FunctorFunction functorFunctionType);

  template <bool calculateGlobals>
  void testFLOPCounterAoSOMP(bool newton3);

  std::array<std::array<double, 3>, 4> getParticlePositions(FunctorFunction functorType);

  template <bool calculateGlobals>
  void testFLOPCounterSoAOMP(bool newton3);

  template <bool calculateGlobals>
  void testFLOPCounterSoAVerletOMP(bool newton3){};
};
