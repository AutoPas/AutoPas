/**
 * @file ForceCalculationTest.h
 * @author F. Gratl
 * @date 13.04.18
 */
#ifdef __AVX__
#pragma once


#include <array>

#include "AutoPasTestBase.h"
#include "autopas/utils/SoA.h"
#include "autopas/options/DataLayoutOption.h"
#include "testingHelpers/commonTypedefs.h"

#include <cstdint>
class MieForceCalculationTest : public AutoPasTestBase {
 public:
  MieForceCalculationTest() = default;

  ~MieForceCalculationTest() override = default;

  /**
   * Creates a test with four particles placed on the corners of a square,
   * executes a LJFunctor and compares the results within a given tolerance.
   *
   * @param n first exponent  of the Mie Potential
   * @param m second exponent of the Mie Potential
   * @param particleSpacing distance between two particles
   * @param cutoff cutoff radius for the LJ
   * @param dataLayoutOption aos or soa
   * @param expectedForces 2D array containing the forces which are expected for
   * each molecule in the end
   * @param tolerance error tolerance
   */
   void testMie(uint16_t n, uint16_t m, double particleSpacing, double cutoff, autopas::DataLayoutOption dataLayoutOption,
              std::array<std::array<double, 3>, 4> expectedForces, double tolerance);
};
#endif