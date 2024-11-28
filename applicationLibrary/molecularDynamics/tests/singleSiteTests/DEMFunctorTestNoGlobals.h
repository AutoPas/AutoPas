/**
 * @file DEMFunctorTestNoGlobals.h
 * @author Joon Kim
 * @date 15.11.2024
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

  constexpr static double radius{1.};
  constexpr static double cutoff{3.};
  constexpr static double elasticStiffness{5.};
  constexpr static double adhesiveStiffness{2.5};
  constexpr static double frictionStiffness{1.};
  constexpr static double normalViscosity{5e-1};
  constexpr static double frictionViscosity{1e-1};
  constexpr static double rollingViscosity{5e-2};
  constexpr static double torsionViscosity{5e-2};
  constexpr static double staticFrictionCoeff{10.};
  constexpr static double dynamicFrictionCoeff{10.};
  constexpr static double rollingFrictionCoeff{10.};
  constexpr static double torsionFrictionCoeff{10.};

  constexpr static double epsilon{1.};
  constexpr static double sigma{1.};
  constexpr static double epsilon2{2.};
  constexpr static double sigma2{2.};

  const std::array<std::array<double, 3>, 2> startingPosNoOverlap = {{{0., 0., 0.}, {2.3, 0, 0}}};
  const std::array<std::array<double, 3>, 2> startingPosOverlap = {{{0., 0., 0.}, {0.1, 0.2, 0.3}}};
  const std::array<std::array<double, 3>, 2> startingVel = {{{3., 3., 3.}, {-3., -3., -3.}}};
  const std::array<std::array<double, 3>, 2> startingAngVel = {{{1., 2., 3.}, {-4., -5., -6.}}};

  constexpr static std::array<double, 3> expectedNormalContactForceOverlap{-3.458327, -6.916653, -10.37498};
  constexpr static std::array<double, 3> expectedNormalContactForceMixingOverlap{-3.458327, -6.916653, -10.37498};

  constexpr static std::array<double, 3> expectedNormalContactForceNoOverlap{0., 0., 0.};
  constexpr static std::array<double, 3> expectedNormalContactForceMixingNoOverlap{0., 0., 0.};

  constexpr static std::array<double, 3> expectedNormalVdWForceOverlap{0.0, 0.0, 0.0};
  constexpr static std::array<double, 3> expectedNormalVdWForceMixingOverlap{0.0, 0.0, 0.0};

  // Deactivate VdW forces
  // constexpr static std::array<double, 3> expectedNormalVdWForceNoOverlap{0.059514258085717045, 0.000000, 0.000000};
  // constexpr static std::array<double, 3> expectedNormalVdWForceMixingNoOverlap{0.23967546841725773, 0.000000,
  // 0.000000};
  constexpr static std::array<double, 3> expectedNormalVdWForceNoOverlap{0., 0.000000, 0.000000};
  constexpr static std::array<double, 3> expectedNormalVdWForceMixingNoOverlap{0., 0.000000, 0.000000};

  constexpr static std::array<double, 3> expectedFrictionForceOverlap{-0.327857142, -0.1157143, 0.1864286};
  constexpr static std::array<double, 3> expectedFrictionForceMixingOverlap{-0.327857, -0.1157143, 0.1864286};

  constexpr static std::array<double, 3> expectedFrictionForceNoOverlap{0., 0., 0.};
  constexpr static std::array<double, 3> expectedFrictionForceMixingNoOverlap{0., 0., 0.};

  constexpr static std::array<double, 3> expectedFrictionTorqueIOverlap{0.035999999, -0.05849999, 0.027000000};
  constexpr static std::array<double, 3> expectedFrictionTorqueIMixingOverlap{0.03599999, -0.0584999, 0.0270000};

  constexpr static std::array<double, 3> expectedFrictionTorqueINoOverlap{0.0, 0.0, 0.0};
  constexpr static std::array<double, 3> expectedFrictionTorqueIMixingNoOverlap{0.0, 0.0, 0.0};

  constexpr static std::array<double, 3> expectedRollingTorqueIOverlap{-0.00075000, -0.00018750, 0.00037500};
  constexpr static std::array<double, 3> expectedRollingTorqueIMixingOverlap{-0.0007500, -0.00018750, 0.00037500};

  constexpr static std::array<double, 3> expectedRollingTorqueINoOverlap{0.0, 0.0, 0.0};
  constexpr static std::array<double, 3> expectedRollingTorqueIMixingNoOverlap{0.0, 0.0, 0.0};

  constexpr static std::array<double, 3> expectedTorsionTorqueIOverlap{-0.00143750, -0.00287500, -0.00431250};
  constexpr static std::array<double, 3> expectedTorsionTorqueIMixingOverlap{-0.00143750, -0.00287500, -0.00431250};

  constexpr static std::array<double, 3> expectedTorsionTorqueINoOverlap{0.0, 0.0, 0.0};
  constexpr static std::array<double, 3> expectedTorsionTorqueIMixingNoOverlap{0.0, 0.0, 0.0};

  constexpr static double absDelta{1e-6};
};