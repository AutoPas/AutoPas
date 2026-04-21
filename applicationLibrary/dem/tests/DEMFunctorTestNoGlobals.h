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

  constexpr static std::array<double, 3> expectedForce{-65.508801145, -130.477602290, -195.356403435};
  constexpr static std::array<double, 3> expectedTorque{0.033812499, -0.061562499, 0.0230625};
  constexpr static std::array<double, 3> expectedTorqueNewton3{0.038187499, -0.055437499, 0.0309375};
  constexpr static double expectedHeatFlux{0.155635714};

  constexpr static std::array<double, 3> expectedForceMixing{-85.533349695, -170.606877763, -255.470138272};
  constexpr static std::array<double, 3> expectedTorqueMixing{-0.016153918, 0.017109176, -0.010162847};
  constexpr static std::array<double, 3> expectedTorqueMixingNewton3{0.104189592, -0.130100257, 0.056145009};
  constexpr static double expectedHeatFluxMixing{0.161654455};

  constexpr static double absDelta{1e-7};
};