/**
 * @file ATMFunctorTestNoGlobals.h
 * @author muehlhaeusser
 * @date 26.09.23
 */

#pragma once

#include <gtest/gtest.h>

#include "ATMFunctorTest.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

template <class FuncType>
class ATMFunctorTestNoGlobals : public ATMFunctorTest {
 public:
  ATMFunctorTestNoGlobals() : ATMFunctorTest() {}

  constexpr static double cutoff{1.};
  constexpr static double nu{1.};
  constexpr static double nu2{1.25};
  constexpr static double nu3{0.1};

  // These values are obtained from the lammps implementation
  const std::array<double, 3> expectedForceP1{-188761.250385809, -188761.250385809, -188761.250385809};
  const std::array<double, 3> expectedForceP2{-95789.2912405596, 94380.6251929043, 284550.541626368};
  const std::array<double, 3> expectedForceP3{284550.541626368, 94380.6251929043, -95789.2912405596};

  const std::array<double, 3> expectedForceMixingP1{-94380.6251929043, -94380.6251929043, -94380.6251929043};
  const std::array<double, 3> expectedForceMixingP2{-47894.6456202798, 47190.3125964521, 142275.270813184};
  const std::array<double, 3> expectedForceMixingP3{142275.270813184, 47190.3125964521, -47894.6456202798};
  constexpr static double absDelta{1e-7};
};
