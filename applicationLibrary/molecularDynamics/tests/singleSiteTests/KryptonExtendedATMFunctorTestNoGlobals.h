/**
 * @file KryptonExtendedATMFunctorTestNoGlobals.h
 * @author D. Martin
 * @date 15.12.25
 */

#pragma once

#include <gtest/gtest.h>

#include "KryptonExtendedATMFunctorTest.h"

template <class FuncType>
class KryptonExtendedATMFunctorTestNoGlobals : public KryptonExtendedATMFunctorTest {
 public:
  KryptonExtendedATMFunctorTestNoGlobals() : KryptonExtendedATMFunctorTest() {}

  constexpr static double cutoff{2.5};

  // These values are obtained from the Krypton AoS implementation
  const std::array<double, 3> expectedForceP1{-8119.8147229097667, -129917.03556655749, -105557.59139782796};
  const std::array<double, 3> expectedForceP2{41123.065379935862, 64958.517783278701, 49927.783543646176};
  const std::array<double, 3> expectedForceP3{-33003.250657026096, 64958.517783278789, 55629.807854181781};
  constexpr static double absDelta{1e-7};
};
