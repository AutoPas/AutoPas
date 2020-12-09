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

namespace LJFunctorCudaTest {

class LJFunctorCudaTest : public AutoPasTestBase,
                          public ::testing::WithParamInterface<
                              std::tuple<bool /*newton3 */, bool /*calculateGlobals*/, bool /*withDeletions*/,
                                         int /*ParticlesCell1*/, int /*ParticlesCell2*/>> {
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
   * LJFunctor works only with MoleculeLJ
   * @param newton3
   */
  template <typename ParticleType, bool calculateGlobals>
  void testLJFunctorVSLJFunctorCudaTwoCells(size_t numParticles, size_t numParticles2, bool useNewton3,
                                            bool withDeletions);

  /**
   * Checks equality of SoALoader, SoAFunctor and SoAExtractor.
   * Expects that particles are loaded and extracted in the same order.
   * In all comparisons first is AVX2, second non-AVX2
   *
   * Checks CudaFunctor(soa, newton3)
   * LJFunctor only works with MoleculeLJ
   * @param newton3
   */
  template <typename ParticleType, bool calculateGlobals>
  void testLJFunctorVSLJFunctorCudaOneCell(size_t numParticles, bool useNewton3, bool withDeletions);

  /**
   * Checks that two non empty SoAs' particles are equal
   * @tparam Particle
   * @param soa1
   * @param soa2
   * @return
   */
  template <class Particle>
  bool SoAParticlesEqual(autopas::SoA<typename Particle::SoAArraysType> &soa1,
                         autopas::SoA<typename Particle::SoAArraysType> &soa2);

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

  const double _cutoff;
  double _epsilon;
  double _sigma;
  const std::array<double, 3> _lowCorner;
  const std::array<double, 3> _highCorner;
};

}  // end namespace LJFunctorCudaTest
#endif  // AUTOPAS_CUDA
