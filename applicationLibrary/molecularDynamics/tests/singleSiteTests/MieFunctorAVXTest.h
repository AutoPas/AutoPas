//
// Created by kay on 03.10.23.
//
/**
 * @file LJFunctorAVXTest.h
 * @author F. Gratl
 * @date 12/17/18
 */

#ifdef __AVX__
#pragma once

#include "AutoPasTestBase.h"
#include "autopas/utils/SoA.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "testingHelpers/commonTypedefs.h"

using MieFunctorAVXTestingTuple = std::tuple<bool /*newton3*/, bool /*doDeleteSomeParticles*/>;

class MieFunctorAVXTest : public AutoPasTestBase, public ::testing::WithParamInterface<MieFunctorAVXTestingTuple> {
 public:
  MieFunctorAVXTest() : AutoPasTestBase() {}

  /**
   *  Maximum error allowed for comparisons.
   */
  constexpr static double _maxError = 1e-12;



  /**
   * Checks equality of SoALoader, SoAFunctorPair and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is MieAVX, second LJ-AVX
   *
   * Checks SoAFunctorPair(soa1, soa2, newton3)
   *
   * @param newton3
   * @param doDeleteSomeParticles
   * @param useUnalignedViews
   */
  void testMieFunctorAVXVSLJFunctorAVXTwoCells(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews);

  /**
   * Checks equality of SoALoader, SoAFunctorSingle and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
* In all comparisons first is MieAVX, second LJ-AVX
   *
   * Checks SoAFunctorSingle(soa, newton3)
   *
   * @param newton3
   * @param doDeleteSomeParticles
   * @param useUnalignedViews
   */
  void testMieFunctorAVXVSLJFunctorAVXOneCell(bool newton3, bool doDeleteSomeParticles, bool useUnalignedViews);

  /**
   * Creates two cells, generates neighbor lists manually and then compares the SoAFunctorVerlet calls.
   * @param newton3
   * @param doDeleteSomeParticles
   */
  void testMieFunctorAVXVSLJFunctorAVXVerlet(bool newton3, bool doDeleteSomeParticles);

  /**
   * Create two cells and compare AoSFunctor
   * @param newton3
   * @param doDeleteSomeParticles
   */
  void testMieFunctorAVXAoS(bool newton3, bool doDeleteSomeParticles);

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
  bool particleEqual(Particle &p1, Particle &p2);

  constexpr static double _cutoff{6.};
  constexpr static double _skinPerTimestep{0.1};
  constexpr static unsigned int _rebuildFrequency{20};
  constexpr static double _interactionLengthSquare{(_cutoff + _skinPerTimestep * _rebuildFrequency) *
                                                   (_cutoff + _skinPerTimestep * _rebuildFrequency)};
  constexpr static double _epsilon{1.};
  constexpr static double _sigma{1.};
  const std::array<double, 3> _lowCorner{0., 0., 0.};
  const std::array<double, 3> _highCorner{6., 6., 6.};
};

#endif  // __AVX__
