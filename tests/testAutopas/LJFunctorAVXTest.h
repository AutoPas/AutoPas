/**
 * @file LJFunctorAVX2Test.h
 * @author F. Gratl
 * @date 12/17/18
 */

#ifdef __AVX__
#pragma once

#include "AutoPasTestBase.h"
#include "autopas/utils/SoA.h"
#include "testingHelpers/commonTypedefs.h"

class LJFunctorAVXTest : public AutoPasTestBase {
 public:
  LJFunctorAVXTest()
      : AutoPasTestBase(), _cutoff{1.}, _epsilon{1.}, _sigma{1.}, _lowCorner{0, 0, 0}, _highCorner{2, 1, 1} {}

  /**
   *  Maximum error allowed for comparisons.
   */
  constexpr static double _maxError = 1e-14;

  /**
   * Checks equality of SoALoader, SoAFunctor and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is AVX2, second non-AVX2
   *
   * Checks SoAFunctor(soa1, soa2, newton3)
   *
   * @param newton3
   */
  void testLJFunctorVSLJFunctorAVXTwoCells(bool newton3);

  /**
   * Checks equality of SoALoader, SoAFunctor and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is AVX2, second non-AVX2
   *
   * Checks SoAFunctor(soa, newton3)
   *
   * @param newton3
   */
  void testLJFunctorVSLJFunctorAVXOneCell(bool newton3);

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
  bool AoSParticlesEqual(FPCell &cell1, FPCell &cell2);

  /**
   * Check that two particles are equal.
   * @param p1
   * @param p2
   * @return
   */
  bool particleEqual(Particle &p1, Particle &p2);

  const double _cutoff;
  const double _epsilon;
  const double _sigma;
  const std::array<double, 3> _lowCorner;
  const std::array<double, 3> _highCorner;
};
#endif  // __AVX__