/**
 * @file LJFunctorTestHWY.h
 * @author Luis Gall
 * @date 23/04/24
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/utils/SoA.h"
#include "molecularDynamicsLibrary/LJFunctorHWY.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"

using VectorizationPattern = autopas::VectorizationPatternOption::Value;

using LJFunctorHWYTestingTuple =
    std::tuple<bool /*mixing*/, bool /*newton3*/, bool /*doDeleteSomeParticles*/, VectorizationPattern>;

class LJFunctorTestHWY : public AutoPasTestBase, public ::testing::WithParamInterface<LJFunctorHWYTestingTuple> {
 public:
  LJFunctorTestHWY() = default;

  constexpr static double _maxError = 1e-8;

  /**
   * Is called before each test case
   */
  void SetUp() override;

  /**
   * Called from SetUp to initialize the particle properties library with common values.
   */
  void setupPPL();

  /**
   * Checks equality of SoALoader, SoAFunctorPair and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is HWY, second Autovec
   *
   * Checks SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam mixing
   * @param newton3
   * @param doDeleteSomeParticles
   * @param useUnalignedViews
   * @param pattern
   */
  template <bool mixing>
  void testLJFunctorvsLJFunctorHWYTwoCells(const bool newton3, const bool doDeleteSomeParticles,
                                           const bool useUnalignedViews, const VectorizationPattern pattern);

  /**
   * Checks equality of SoALoader, SoAFunctorSingle and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is HWY, second Autovec
   *
   * Checks SoAFunctorSingle(soa, newton3)
   *
   * @tparam mixing
   * @param newton3
   * @param doDeleteSomeParticles
   * @param useUnalignedViews
   * @param pattern
   */
  template <bool mixing>
  void testLJFunctorvsLJFunctorHWYOneCell(const bool newton3, const bool doDeleteSomeParticles,
                                          const bool useUnalignedViews, const VectorizationPattern pattern);

  /**
   * Creates two cells, generates neighbor lists manually and then compares the SoAFunctorVerlet calls.
   *
   * @tparam mixing
   * @param newton3
   * @param doDeleteSomeParticles
   */
  template <bool mixing>
  void testLJFunctorvsLJFunctorHWYVerlet(const bool newton3, const bool doDeleteSomeParticles);

  /**
   * Create two cells and compare the HWY AoSFunctor to the Autovec functor.
   *
   * @tparam mixing
   * @param newton3
   * @param doDeleteSomeParticles
   */
  template <bool mixing>
  void testLJFunctorvsLJFunctorHWYAoS(const bool newton3, const bool doDeleteSomeParticles);

  /**
   * Checks that two non empty SoAs' particles are equal
   *
   * @tparam SoAType
   * @param soa1
   * @param soa2
   * @return true if the two SoAs' particles are equal
   */
  template <class SoAType>
  bool checkSoAParticlesAreEqual(const autopas::SoA<SoAType> &soa1, const autopas::SoA<SoAType> &soa2);

  /**
   * Checks if the two cells cell 1 and cell 2 are equal. They are equal, if the size is the same and all contained
   * particles are equal.
   *
   * @param cell1
   * @param cell2
   * @return if the two cells cell 1 and cell 2 are equal.
   */
  bool checkAoSParticlesAreEqual(const FMCell &cell1, const FMCell &cell2);

  /**
   * Checks if the two particles p1 and p2 are equal.
   *
   * @param p1 particle 1
   * @param p2 particle 2
   * @return if the two particles p1 and p2 are equal.
   */
  bool checkParticlesAreEqual(const Molecule &p1, const Molecule &p2);

  // Values for the test scenario.
  constexpr static double _cutoff{6.};
  constexpr static double _skin{2.};
  constexpr static double _interactionLengthSquare{(_cutoff + _skin) * (_cutoff + _skin)};
  constexpr static double _epsilon{1.};
  constexpr static double _sigma{1.};
  const std::array<double, 3> _lowCorner{0., 0., 0.};
  const std::array<double, 3> _highCorner{8., 8., 8.};

  // Global PPL for all tests.
  ParticlePropertiesLibrary<double, size_t> _PPL{_cutoff};
};
