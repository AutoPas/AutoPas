/**
 * @file ParticleMatcher.h
 * @author J. Schuhmacher
 * @date 11.02.26
 */

#pragma once

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "autopas/utils/FloatComparison.h"

/**
 * Check if two particle-like objects are strictly equal.
 * Strict equality means that position, velocity, force, and ID are compared using `operator==`
 * on the respective return types.
 *
 * @tparam L Particle-like type of the left-hand side.
 * @tparam R Particle-like type of the right-hand side.
 * @param lhs Left-hand side particle.
 * @param rhs Right-hand side particle.
 * @return True if `getR()`, `getV()`, `getF()`, and `getID()` are all strictly equal; false otherwise.
 */
bool equalParticles(const auto &lhs, const auto &rhs) {
  return lhs.getR() == rhs.getR() && lhs.getV() == rhs.getV() && lhs.getF() == rhs.getF() && lhs.getID() == rhs.getID();
}

/**
 * GoogleTest Matcher for strict particle equality.
 * Compares `getR()`, `getV()`, `getF()`, and `getID()` via `operator==`.
 *
 * Intended usage:
 * @code{.cpp}
 * EXPECT_THAT(std::tie(actual, expected), ParticleEq());
 * @endcode
 */
MATCHER(ParticleEq, "Comparing if two particles are strictly equal to each other") {
  const auto &lhs = std::get<0>(arg);
  const auto &rhs = std::get<1>(arg);

  return lhs.getR() == rhs.getR() && lhs.getV() == rhs.getV() && lhs.getF() == rhs.getF() && lhs.getID() == rhs.getID();
}

/**
 * Check if two particle-like objects are almost equal using relative comparisons.
 *
 * Position, velocity, and force are compared using `autopas::utils::almostEqualRelative(...)`
 * with a configurable relative tolerance. Additionally, `getTypeId()` must match exactly.
 *
 * @tparam L Particle-like type of the left-hand side.
 * @tparam R Particle-like type of the right-hand side.
 * @param lhs Left-hand side particle.
 * @param rhs Right-hand side particle.
 * @param epsilon Relative tolerance used for the floating-point comparisons.
 * @return True if `getR()`, `getV()`, and `getF()` are almost equal (relative) and `getTypeId()` matches.
 */
bool almostEqualParticles(const auto &lhs, const auto &rhs, double epsilon = autopas::utils::EPSILON_ALMOST_EQUAL) {
  return autopas::utils::almostEqualRelative(lhs.getR(), rhs.getR(), epsilon) &&
         autopas::utils::almostEqualRelative(lhs.getV(), rhs.getV(), epsilon) &&
         autopas::utils::almostEqualRelative(lhs.getF(), rhs.getF(), epsilon) && lhs.getTypeId() == rhs.getTypeId();
}

/**
 * GoogleTest Matcher for approximate particle equality with configurable relative tolerance.
 * Uses `almostEqualParticles(lhs, rhs, epsilon)` internally.
 *
 * Intended usage:
 * @code{.cpp}
 * EXPECT_THAT(std::tie(actual, expected), ParticleAlmostEq(1e-12));
 * @endcode
 *
 * @param epsilon Relative tolerance forwarded to the underlying comparisons.
 */
MATCHER_P(ParticleAlmostEq, epsilon,
          "Comparing if two particles are almost equal to each other given a relative epsilon") {
  const auto &lhs = std::get<0>(arg);
  const auto &rhs = std::get<1>(arg);
  return almostEqualParticles(lhs, rhs, epsilon);
}

/**
 * Check if two particle-like objects are almost equal using an ULP distance.
 * Position, velocity, and force are compared using `autopas::utils::almostEqualUlps(...)`
 * with a configurable maximum ULP distance. Additionally, `getTypeId()` must match exactly.
 *
 * @tparam L Particle-like type of the left-hand side.
 * @tparam R Particle-like type of the right-hand side.
 * @param lhs Left-hand side particle.
 * @param rhs Right-hand side particle.
 * @param ulpDistance Maximum allowed distance in ULPs.
 * @return True if `getR()`, `getV()`, and `getF()` are almost equal (ULPs) and `getTypeId()` matches.
 */
bool almostEqualParticlesUlps(const auto &lhs, const auto &rhs,
                              unsigned int ulpDistance = autopas::utils::MAX_ULP_DISTANCE) {
  return autopas::utils::almostEqualUlps(lhs.getR(), rhs.getR(), ulpDistance) &&
         autopas::utils::almostEqualUlps(lhs.getV(), rhs.getV(), ulpDistance) &&
         autopas::utils::almostEqualUlps(lhs.getF(), rhs.getF(), ulpDistance) && lhs.getTypeId() == rhs.getTypeId();
}

/**
 * GoogleTest Matcher for approximate particle equality with configurable ULP distance.
 * Uses `almostEqualParticlesUlps(lhs, rhs, ulpDistance)` internally.
 *
 * Intended usage:
 * @code{.cpp}
 * EXPECT_THAT(std::tie(actual, expected), ParticleAlmostEqUlps(8));
 * @endcode
 *
 * @param ulpDistance Maximum allowed distance in ULPs forwarded to the underlying comparisons.
 */
MATCHER_P(ParticleUlpsEq, ulpDistance,
          "Comparing if two particles are almost equal to each other given an ULP distance") {
  const auto &lhs = std::get<0>(arg);
  const auto &rhs = std::get<1>(arg);
  return almostEqualParticlesUlps(lhs, rhs, ulpDistance);
}