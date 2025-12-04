/**
 * @file ATMFunctorHWYTest.h
 * @author D. Martin
 * @date 11/17/25
 */

#pragma once

#include "AutoPasTestBase.h"
#include "autopas/utils/SoA.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"

using ATMFunctorAVXTestingTuple = std::tuple<bool /*mixing*/, bool /*newton3*/, bool /*doDeleteSomeParticles*/>;

class ATMFunctorHWYTest : public AutoPasTestBase, public ::testing::WithParamInterface<ATMFunctorAVXTestingTuple> {
 public:
  ATMFunctorHWYTest() : AutoPasTestBase() {}

  /**
   *  Maximum error allowed for comparisons.
   */
  constexpr static double _maxError = 1e-12;

  /**
   * Checks equality of SoALoader, SoAFunctorTriple and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is HWY, second non-HWY
   *
   * Checks SoAFunctorTriple(soa1, soa2, soa3, newton3)
   *
   * @tparam mixing
   * @param newton3
   * @param doDeleteSomeParticles
   * @param useUnalignedViews
   */
  template <bool mixing>
  void testATMFunctorVSATMFunctorHWYThreeCells(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews);

  /**
   * Checks equality of SoALoader, SoAFunctorPair and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is HWY, second non-HWY
   *
   * Checks SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam mixing
   * @param newton3
   * @param doDeleteSomeParticles
   * @param useUnalignedViews
   */
  template <bool mixing>
  void testATMFunctorVSATMFunctorHWYTwoCells(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews);

  /**
   * Checks equality of SoALoader, SoAFunctorSingle and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is HWY, second non-HWY
   *
   * Checks SoAFunctorSingle(soa, newton3)
   *
   * @tparam mixing
   * @param newton3
   * @param doDeleteSomeParticles
   * @param useUnalignedViews
   */
  template <bool mixing>
  void testATMFunctorVSATMFunctorHWYOneCell(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews);

  /**
   * Create two cells and compare AoSFunctor
   * @tparam mixing
   * @param newton3
   * @param doDeleteSomeParticles
   */
  template <bool mixing>
  void testATMFunctorVSATMFunctorHWYAoS(bool newton3, bool doDeleteSomeParticles);

  /**
   * Checks that two non empty SoAs' particles are equal
   * @tparam SoAType
   * @param soa1
   * @param soa2
   * @return
   */
  template <class SoAType>
  bool SoAParticlesEqual(autopas::SoA<SoAType> &soa1, autopas::SoA<SoAType> &soa2);

  /**
   * Check that two non empty AoSs' (=Cells) particles are equal.
   * @param cell1
   * @param cell2
   * @return
   */
  bool AoSParticlesEqual(FMCell &cell1, FMCell &cell2);

  /**
   * Check that two particles are equal.
   * @param p1
   * @param p2
   * @return
   */
  bool particleEqual(Molecule &p1, Molecule &p2);

  constexpr static double _cutoff{6.};
  constexpr static double _skin{2.};
  constexpr static unsigned int _rebuildFrequency{20};
  constexpr static double _interactionLengthSquare{(_cutoff + _skin) * (_cutoff + _skin)};
  // Parameters for mixing = false
  constexpr static double _nu{1.};

  const std::array<double, 3> _lowCorner{0., 0., 0.};
  const std::array<double, 3> _highCorner{6., 6., 6.};
};