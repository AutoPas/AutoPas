/**
 * @file LJFunctorCudaTest.h
 * @author jspahl
 * @date 2/11/18
 */

#pragma once
#if defined(AUTOPAS_CUDA)

#include "AutoPasTestBase.h"
#include "autopas/utils/SoA.h"
#include "testingHelpers/commonTypedefs.h"

class LJFunctorCudaTest : public AutoPasTestBase,
                          public ::testing::WithParamInterface<std::tuple<bool, bool, int, int>> {
 public:
  LJFunctorCudaTest()
      : AutoPasTestBase(), _cutoff{1.}, _epsilon{2}, _sigma{0.05}, _lowCorner{0, 0, 0}, _highCorner{2, 1, 1} {}

  /**
   *  Maximum error allowed for comparisons.
   */
  constexpr static double _maxError = 1e-12;

  /**
   * Checks equality of SoALoader, SoAFunctor and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is AVX2, second non-AVX2
   *
   * Checks CudaFunctor(soa1, soa2, newton3)
   *
   * @param newton3
   */
  template <typename ParticleType, bool useNewton3, bool calculateGlobals>
  void testLJFunctorVSLJFunctorCudaTwoCells(size_t numParticles, size_t numParticles2);

  /**
   * Checks equality of SoALoader, SoAFunctor and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is AVX2, second non-AVX2
   *
   * Checks CudaFunctor(soa, newton3)
   *
   * @param newton3
   */
  template <typename ParticleType, bool useNewton3, bool calculateGlobals>
  void testLJFunctorVSLJFunctorCudaOneCell(size_t numParticles);

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

#endif  // AUTOPAS_CUDA
