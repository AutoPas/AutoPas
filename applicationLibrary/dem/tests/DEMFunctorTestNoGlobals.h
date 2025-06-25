/**
 * @file DEMFunctorTestNoGlobals.h
 * @author Joon Kim
 * @date 25.06.25
 */

#pragma once

#include <gtest/gtest.h>

#include "DEMFunctorTest.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

template <class FuncType>
class DEMFunctorTestNoGlobals : public DEMFunctorTest {
 public:
  DEMFunctorTestNoGlobals() : DEMFunctorTest() {}

  constexpr static double epsilon{1e-6};
  constexpr static double cutoff{3.};
  constexpr static double mass{1.};
  constexpr static double radius1{1.};
  constexpr static double radius2{1.5};
  constexpr static double specificHeat{1.};

  constexpr static double elasticStiffness{150.};

  constexpr static double normalViscosity{1e-3};
  constexpr static double frictionViscosity{1e-1};
  constexpr static double rollingViscosity{5e-2};
  constexpr static double torsionViscosity{5e-2};

  constexpr static double staticFrictionCoeff{10.};
  constexpr static double dynamicFrictionCoeff{10.};
  constexpr static double rollingFrictionCoeff{10.};
  constexpr static double torsionFrictionCoeff{10.};

  constexpr static double heatConductivity{1.0};
  constexpr static double heatGenerationFactor{1e-1};

  const std::array<std::array<double, 3>, 2> startingPosOverlap = {{{0., 0., 0.}, {0.1, 0.2, 0.3}}};
  const std::array<std::array<double, 3>, 2> startingVel = {{{3., 3., 3.}, {-3., -3., -3.}}};
  const std::array<std::array<double, 3>, 2> startingAngVel = {{{1., 2., 3.}, {-4., -5., -6.}}};

  constexpr static std::array<double, 3> expectedForce{-65.5088011, -130.4776022, -195.3564034};
  constexpr static std::array<double, 3> expectedTorque{0.0338124, -0.0615624, 0.0230625};
  constexpr static std::array<double, 3> expectedTorqueNewton3{0.0381874, -0.0554374, 0.0309375};
  constexpr static double expectedHeatFlux{0.1556357};

  constexpr static std::array<double, 3> expectedForceMixing{-85.5333496, -170.6068777, -255.4701382};
  constexpr static std::array<double, 3> expectedTorqueMixing{-0.0161539, 0.0171091, -0.0101628};
  constexpr static std::array<double, 3> expectedTorqueMixingNewton3{0.1041895, -0.1301002, 0.0561450};
  constexpr static double expectedHeatFluxMixing{0.1616544};

  constexpr static double absDelta{1e-7};
};