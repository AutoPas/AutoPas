/**
 * @file LJFunctorAVX.h
 *
 * @date 17 Jan 2018
 * @author F. Gratl
 */
#pragma once

#include <immintrin.h>
#include <array>
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This Version is implemented using AVX intrinsics.
 * @tparam Particle The type of particle.
 * @tparam ParticleCell The type of particlecell.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, class ParticleCell, FunctorN3Modes useNewton3 = FunctorN3Modes::Both,
          bool calculateGlobals = false, bool relevantForTuning = true>
class LJFunctorAVX : public Functor<Particle, ParticleCell, typename Particle::SoAArraysType> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorAVX() = delete;

  /**
   * Constructor, which sets the global values, i.e. cutoff, epsilon, sigma and shift.
   * @param cutoff
   * @param epsilon
   * @param sigma
   * @param shift
   * @param lowCorner Lower corner of the local simulation domain.
   * @param highCorner Upper corner of the local simulation domain.
   * @param duplicatedCalculation Defines whether duplicated calculations are happening across processes / over the
   * simulation boundary. e.g. eightShell: false, fullShell: true.
   */
  explicit LJFunctorAVX(double cutoff, double epsilon, double sigma, double shift,
                        std::array<double, 3> lowCorner = {0., 0., 0.}, std::array<double, 3> highCorner = {0., 0., 0.},
                        bool duplicatedCalculation = false)
#ifdef __AVX__
      : _one{_mm256_set1_pd(1.)},
        _masks{
            _mm256_set_epi64x(0, 0, 0, -1),
            _mm256_set_epi64x(0, 0, -1, -1),
            _mm256_set_epi64x(0, -1, -1, -1),
        },
        _cutoffsquare{_mm256_set1_pd(cutoff * cutoff)},
        _epsilon24{_mm256_set1_pd(epsilon * 24.0)},
        _sigmasquare{_mm256_set1_pd(sigma * sigma)},
        _shift6{shift * 6.0},
        _upotSum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _duplicatedCalculations{duplicatedCalculation},
        _lowCorner(lowCorner),
        _highCorner(highCorner),
        _postProcessed{false} {

    if (calculateGlobals and duplicatedCalculation) {
      if (lowCorner == highCorner) {
        throw utils::ExceptionHandler::AutoPasException(
            "Please specify the lowCorner and highCorner properly if calculateGlobals and duplicatedCalculation are "
            "set to true.");
      }
    }
    if (calculateGlobals) {
      _aosThreadData.resize(autopas_get_max_threads());
    }
  }
#else
      : _one{0}, _masks{0, 0, 0}, _cutoffsquare{0}, _epsilon24{0}, _sigmasquare{0} {
    utils::ExceptionHandler::exception("AutoPas was compiled without AVX support!");
  }
#endif

  bool isRelevantForTuning() override { return relevantForTuning; }

  bool allowsNewton3() override {
    return useNewton3 == FunctorN3Modes::Newton3Only or useNewton3 == FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() override {
    return useNewton3 == FunctorN3Modes::Newton3Off or useNewton3 == FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
    utils::ExceptionHandler::exception("LJFunctorAVX.AoSFunctor() not implemented!");
  }

  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa, bool newton3)
   * This functor ignores the newton3 value, as we do not expect any benefit from disabling newton3.
   */
  void SoAFunctor(SoA<SoAArraysType> &soa, bool newton3) override {
#ifdef __AVX__
    if (soa.getNumParticles() == 0) return;

    double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    double *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    double *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d upotSum = _mm256_setzero_pd();

    if (calculateGlobals) {
      // Checks if the cell is a halo cell, if it is, we skip it.
      // We cannot do this in normal cases (where we do not calculate globals), as _lowCorner and _highCorner are not
      // set. (as of 23.11.2018)
      bool isHaloCell = false;
      isHaloCell |= xptr[0] < _lowCorner[0] || xptr[0] >= _highCorner[0];
      isHaloCell |= yptr[0] < _lowCorner[1] || yptr[0] >= _highCorner[1];
      isHaloCell |= zptr[0] < _lowCorner[2] || zptr[0] >= _highCorner[2];
      if (isHaloCell) {
        return;
      }
    }

    // reverse outer loop s.th. inner loop always beginns at aligned array start
    // typecast to detect underflow
    for (size_t i = soa.getNumParticles() - 1; (long)i >= 0; --i) {
      __m256d fxacc = _mm256_setzero_pd();
      __m256d fyacc = _mm256_setzero_pd();
      __m256d fzacc = _mm256_setzero_pd();

      const __m256d x1 = _mm256_broadcast_sd(&xptr[i]);
      const __m256d y1 = _mm256_broadcast_sd(&yptr[i]);
      const __m256d z1 = _mm256_broadcast_sd(&zptr[i]);

      // floor soa numParticles to multiple of vecLength
      unsigned int j = 0;
      for (; j < (i & ~(vecLength - 1)); j += 4) {
        SoAKernel<true, false>(j, x1, y1, z1, xptr, yptr, zptr, fxptr, fyptr, fzptr, fxacc, fyacc, fzacc, &virialSumX,
                               &virialSumY, &virialSumZ, &upotSum);
      }
      const int rest = (int)(i & (vecLength - 1));
      if (rest > 0)
        SoAKernel<true, true>(j, x1, y1, z1, xptr, yptr, zptr, fxptr, fyptr, fzptr, fxacc, fyacc, fzacc, &virialSumX,
                              &virialSumY, &virialSumZ, &upotSum, rest);

      // horizontally reduce fDacc to sumfD
      const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
      const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

      const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
      const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

      const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
      const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

      const union {
        __m128d reg;
        double arr[2];
      } sumfxfyVEC = {.reg = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh)};
      const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

      const double sumfx = sumfxfyVEC.arr[0];
      const double sumfy = sumfxfyVEC.arr[1];
      const double sumfz = _mm_cvtsd_f64(sumfzVEC);

      fxptr[i] += sumfx;
      fyptr[i] += sumfy;
      fzptr[i] += sumfz;
    }

    if (calculateGlobals) {
      const int threadnum = autopas_get_thread_num();

      // horizontally reduce virialSumX and virialSumY
      const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
      const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
      const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
      const __m128d hSumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      // horizontally reduce virialSumZ and upotSum
      const __m256d hSumVirialzUpot = _mm256_hadd_pd(virialSumZ, upotSum);
      const __m128d hSumVirialzUpotLow = _mm256_extractf128_pd(hSumVirialzUpot, 0);
      const __m128d hSumVirialzUpotHigh = _mm256_extractf128_pd(hSumVirialzUpot, 1);
      const __m128d hSumVirialzUpotVec = _mm_add_pd(hSumVirialzUpotHigh, hSumVirialzUpotLow);

      // globals = {virialX, virialY, virialZ, uPot}
      double globals[4];
      _mm_store_pd(&globals[0], hSumVirialxyVec);
      _mm_store_pd(&globals[2], hSumVirialzUpotVec);

      // if newton3 is false, then we divide by 2 later on, so we multiply by two here (very hacky, but needed for
      // AoS)
      _aosThreadData[threadnum].virialSum[0] += globals[0] * (newton3 ? 1 : 2);
      _aosThreadData[threadnum].virialSum[1] += globals[1] * (newton3 ? 1 : 2);
      _aosThreadData[threadnum].virialSum[2] += globals[2] * (newton3 ? 1 : 2);
      _aosThreadData[threadnum].upotSum += globals[3] * (newton3 ? 1 : 2);
    }
#endif
  }

 private:
  template <bool newton3, bool masked>
  inline void SoAKernel(size_t j, const __m256d &x1, const __m256d &y1, const __m256d &z1,
                        double *const __restrict__ &x2ptr, double *const __restrict__ &y2ptr,
                        double *const __restrict__ &z2ptr, double *const __restrict__ &fx2ptr,
                        double *const __restrict__ &fy2ptr, double *const __restrict__ &fz2ptr, __m256d &fxacc,
                        __m256d &fyacc, __m256d &fzacc, __m256d *virialSumX, __m256d *virialSumY, __m256d *virialSumZ,
                        __m256d *upotSum, const unsigned int rest = 0) {
#ifdef __AVX__
    const __m256d x2 = masked ? _mm256_maskload_pd(&x2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&x2ptr[j]);
    const __m256d y2 = masked ? _mm256_maskload_pd(&y2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&y2ptr[j]);
    const __m256d z2 = masked ? _mm256_maskload_pd(&z2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&z2ptr[j]);

    const __m256d drx = _mm256_sub_pd(x1, x2);
    const __m256d dry = _mm256_sub_pd(y1, y2);
    const __m256d drz = _mm256_sub_pd(z1, z2);

    const __m256d drx2 = _mm256_mul_pd(drx, drx);
    const __m256d dry2 = _mm256_mul_pd(dry, dry);
    const __m256d drz2 = _mm256_mul_pd(drz, drz);

    const __m256d dr2PART = _mm256_add_pd(drx2, dry2);
    const __m256d dr2 = _mm256_add_pd(dr2PART, drz2);

    // _CMP_LE_OS == Less-Equal-then (ordered, signaling)
    // signaling = throw error if NaN is encountered
    // dr2 <= _cutoffsquare ? 0xFFFFFFFFFFFFFFFF : 0
    const __m256d cutoffMask = _mm256_cmp_pd(dr2, _cutoffsquare, _CMP_LE_OS);

    const __m256d invdr2 = _mm256_div_pd(_one, dr2);
    const __m256d lj2 = _mm256_mul_pd(_sigmasquare, invdr2);
    const __m256d lj4 = _mm256_mul_pd(lj2, lj2);
    const __m256d lj6 = _mm256_mul_pd(lj2, lj4);
    const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
    const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);
    const __m256d lj12m6alj12 = _mm256_add_pd(lj12m6, lj12);
    const __m256d lj12m6alj12e = _mm256_mul_pd(lj12m6alj12, _epsilon24);
    const __m256d fac = _mm256_mul_pd(lj12m6alj12e, invdr2);

    const __m256d facMasked = masked
                                  ? _mm256_and_pd(fac, _mm256_and_pd(cutoffMask, _mm256_castsi256_pd(_masks[rest - 1])))
                                  : _mm256_and_pd(fac, cutoffMask);

    const __m256d fx = _mm256_mul_pd(drx, facMasked);
    const __m256d fy = _mm256_mul_pd(dry, facMasked);
    const __m256d fz = _mm256_mul_pd(drz, facMasked);

    fxacc = _mm256_add_pd(fxacc, fx);
    fyacc = _mm256_add_pd(fyacc, fy);
    fzacc = _mm256_add_pd(fzacc, fz);

    // if newton 3 is used subtract fD from particle j
    if (newton3) {
      const __m256d fx2 = masked ? _mm256_maskload_pd(&fx2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&fx2ptr[j]);
      const __m256d fy2 = masked ? _mm256_maskload_pd(&fy2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&fy2ptr[j]);
      const __m256d fz2 = masked ? _mm256_maskload_pd(&fz2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&fz2ptr[j]);

      const __m256d fx2new = _mm256_sub_pd(fx2, fx);
      const __m256d fy2new = _mm256_sub_pd(fy2, fy);
      const __m256d fz2new = _mm256_sub_pd(fz2, fz);

      masked ? _mm256_maskstore_pd(&fx2ptr[j], _masks[rest - 1], fx2new) : _mm256_store_pd(&fx2ptr[j], fx2new);
      masked ? _mm256_maskstore_pd(&fy2ptr[j], _masks[rest - 1], fy2new) : _mm256_store_pd(&fy2ptr[j], fy2new);
      masked ? _mm256_maskstore_pd(&fz2ptr[j], _masks[rest - 1], fz2new) : _mm256_store_pd(&fz2ptr[j], fz2new);
    }

    if (calculateGlobals) {
      // Global Virial
      const __m256d virialx = _mm256_mul_pd(fx, drx);
      const __m256d virialy = _mm256_mul_pd(fy, dry);
      const __m256d virialz = _mm256_mul_pd(fz, drz);

      // Global Potential
      const __m256d upotPART1 = _mm256_mul_pd(_epsilon24, lj12m6);
      const __m256d shift6V = _mm256_set1_pd(_shift6);
      const __m256d upotPART2 = _mm256_add_pd(shift6V, upotPART1);
      const __m256d upot =
          masked ? _mm256_and_pd(upotPART2, _mm256_and_pd(cutoffMask, _mm256_castsi256_pd(_masks[rest - 1])))
                 : _mm256_and_pd(upotPART2, cutoffMask);

      *virialSumX = _mm256_add_pd(*virialSumX, virialx);
      *virialSumY = _mm256_add_pd(*virialSumY, virialy);
      *virialSumZ = _mm256_add_pd(*virialSumZ, virialz);
      *upotSum = _mm256_add_pd(*upotSum, upot);
    }
#endif
  }


 private:
  template <bool newton3, bool masked>
  inline void SoAKernelNew(const __m256d &x1, const __m256d &y1, const __m256d &z1,
                        const __m256d &x2, const __m256d &y2, const __m256d &z2,
                        __m256d *fx1, __m256d *fy1, __m256d *fz1,
                        __m256d *fx2, __m256d *fy2, __m256d *fz2,
                        __m256d *virialSumX, __m256d *virialSumY, __m256d *virialSumZ,
                        __m256d *upotSum, __m256d &mask) {
#ifdef __AVX__

    const __m256d drx = _mm256_sub_pd(x1, x2);
    const __m256d dry = _mm256_sub_pd(y1, y2);
    const __m256d drz = _mm256_sub_pd(z1, z2);

    const __m256d drx2 = _mm256_mul_pd(drx, drx);
    const __m256d dry2 = _mm256_mul_pd(dry, dry);
    const __m256d drz2 = _mm256_mul_pd(drz, drz);

    const __m256d dr2PART = _mm256_add_pd(drx2, dry2);
    const __m256d dr2 = _mm256_add_pd(dr2PART, drz2);

    // _CMP_LE_OS == Less-Equal-then (ordered, signaling)
    // signaling = throw error if NaN is encountered
    // dr2 <= _cutoffsquare ? 0xFFFFFFFFFFFFFFFF : 0
    const __m256d cutoffMask = _mm256_cmp_pd(dr2, _cutoffsquare, _CMP_LE_OS);

    const __m256d invdr2 = _mm256_div_pd(_one, dr2);
    const __m256d lj2 = _mm256_mul_pd(_sigmasquare, invdr2);
    const __m256d lj4 = _mm256_mul_pd(lj2, lj2);
    const __m256d lj6 = _mm256_mul_pd(lj2, lj4);
    const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
    const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);
    const __m256d lj12m6alj12 = _mm256_add_pd(lj12m6, lj12);
    const __m256d lj12m6alj12e = _mm256_mul_pd(lj12m6alj12, _epsilon24);
    const __m256d fac = _mm256_mul_pd(lj12m6alj12e, invdr2);

      //TODO #####################################################                v
    const __m256d facMasked = masked
                                  ? _mm256_and_pd(fac, _mm256_and_pd(cutoffMask, mask))
                                  : _mm256_and_pd(fac, cutoffMask);

    const __m256d fx = _mm256_mul_pd(drx, facMasked);
    const __m256d fy = _mm256_mul_pd(dry, facMasked);
    const __m256d fz = _mm256_mul_pd(drz, facMasked);

    *fx1 = _mm256_add_pd(*fx1, fx);
    *fy1 = _mm256_add_pd(*fy1, fy);
    *fz1 = _mm256_add_pd(*fz1, fz);

    // if newton 3 is used subtract fD from particle j
    if (newton3) {
      *fx2 = _mm256_sub_pd(*fx2, fx);
      *fy2 = _mm256_sub_pd(*fy2, fy);
      *fz2 = _mm256_sub_pd(*fz2, fz);
    }

    if (calculateGlobals) {
      // Global Virial
      const __m256d virialx = _mm256_mul_pd(fx, drx);
      const __m256d virialy = _mm256_mul_pd(fy, dry);
      const __m256d virialz = _mm256_mul_pd(fz, drz);

      // Global Potential
      const __m256d upotPART1 = _mm256_mul_pd(_epsilon24, lj12m6);
      const __m256d shift6V = _mm256_set1_pd(_shift6);
      const __m256d upotPART2 = _mm256_add_pd(shift6V, upotPART1);

      //TODO #####################################################     v
      const __m256d upot =
          masked ? _mm256_and_pd(upotPART2, _mm256_and_pd(cutoffMask, mask))
                 : _mm256_and_pd(upotPART2, cutoffMask);

      *virialSumX = _mm256_add_pd(*virialSumX, virialx);
      *virialSumY = _mm256_add_pd(*virialSumY, virialy);
      *virialSumZ = _mm256_add_pd(*virialSumZ, virialz);
      *upotSum = _mm256_add_pd(*upotSum, upot);
    }
#endif
  }

 private:
  void permute(int* perm, auto* from, auto* to) {
    for (int a = 0; a < 4; a++) {
      to[perm[a]] = from[a];
    }
  }
 private:
  void invpermute(int* perm, auto* from, auto* to, int max = 4) {
    for (int a = 0; a < max; a++) {
      to[a] = from[perm[a]];
    }
  }

 private:
  template <bool masked1, bool masked2>
  inline void SoAPermKernel(double* const x1ptr, double* const y1ptr, double* const z1ptr, double* fx1ptr, double* fy1ptr, double* fz1ptr, int r1, 
                     double* const x2ptr, double* const y2ptr, double* const z2ptr, double* fx2ptr, double* fy2ptr, double* fz2ptr, int r2,
                     __m256d *virialSumX, __m256d *virialSumY, __m256d *virialSumZ, __m256d *upotSum, bool newton3){

    int permArray[4][4] = {{0,1,2,3},{1,2,3,0},{2,3,0,1},{3,0,1,2}};


    //unsigned long masks [4][4] = {{0,0,0,0},{ULONG_MAX,ULONG_MAX,ULONG_MAX,0},{ULONG_MAX,ULONG_MAX,0,0},{ULONG_MAX,0,0,0}};
    unsigned long masks [4][4] = {{0,0,0,0},{ULONG_MAX,0,0,0},{ULONG_MAX,ULONG_MAX,0,0},{ULONG_MAX,ULONG_MAX,ULONG_MAX,0}};
    unsigned long *mask1ptr = masks[r1]; //depends on r1
    unsigned long *mask2ptr = masks[r2]; //depends on r2

    /*
     *if (masked1) {
     *  std::cout << "r1: " << r1 << " " << mask1ptr[0] << " " << mask1ptr[1] << " " << mask1ptr[2] << " " << mask1ptr[3] << " " << std::endl;
     *}
     *if (masked2) {
     *  std::cout << "r2: " << r2 << " " << mask2ptr[0] << " " << mask2ptr[1] << " " << mask2ptr[2] << " " << mask2ptr[3] << " " << std::endl;
     *}
     */

    __m256i mask1 = _mm256_loadu_si256((__m256i *) &mask1ptr[0]);

    __m256d x1 = masked1 ? _mm256_maskload_pd(&x1ptr[0], mask1) : _mm256_load_pd(&x1ptr[0]);
    __m256d y1 = masked1 ? _mm256_maskload_pd(&y1ptr[0], mask1) : _mm256_load_pd(&y1ptr[0]);
    __m256d z1 = masked1 ? _mm256_maskload_pd(&z1ptr[0], mask1) : _mm256_load_pd(&z1ptr[0]);

    __m256d fx1 = masked1 ? _mm256_maskload_pd(&fx1ptr[0], mask1) : _mm256_load_pd(&fx1ptr[0]);
    __m256d fy1 = masked1 ? _mm256_maskload_pd(&fy1ptr[0], mask1) : _mm256_load_pd(&fy1ptr[0]);
    __m256d fz1 = masked1 ? _mm256_maskload_pd(&fz1ptr[0], mask1) : _mm256_load_pd(&fz1ptr[0]);

    for (int a = 0; a < 4; a++) {
      double x2ptrPerm[4];
      double y2ptrPerm[4];
      double z2ptrPerm[4];
      permute(permArray[a], x2ptr, x2ptrPerm);
      permute(permArray[a], y2ptr, y2ptrPerm);
      permute(permArray[a], z2ptr, z2ptrPerm);

      double fx2ptrPerm[4];
      double fy2ptrPerm[4];
      double fz2ptrPerm[4];
      permute(permArray[a], fx2ptr, fx2ptrPerm);
      permute(permArray[a], fy2ptr, fy2ptrPerm);
      permute(permArray[a], fz2ptr, fz2ptrPerm);

      unsigned long mask2ptrPerm[4];
      permute(permArray[a], mask2ptr, mask2ptrPerm);
      __m256i mask2 = _mm256_loadu_si256((__m256i *) &mask2ptrPerm[0]);


      __m256d x2 = masked2 ? _mm256_maskload_pd(&x2ptrPerm[0], mask2) : _mm256_load_pd(&x2ptrPerm[0]);
      __m256d y2 = masked2 ? _mm256_maskload_pd(&y2ptrPerm[0], mask2) : _mm256_load_pd(&y2ptrPerm[0]);
      __m256d z2 = masked2 ? _mm256_maskload_pd(&z2ptrPerm[0], mask2) : _mm256_load_pd(&z2ptrPerm[0]);

      __m256d fx2 = masked2 ? _mm256_maskload_pd(&fx2ptrPerm[0], mask2) : _mm256_load_pd(&fx2ptrPerm[0]);
      __m256d fy2 = masked2 ? _mm256_maskload_pd(&fy2ptrPerm[0], mask2) : _mm256_load_pd(&fy2ptrPerm[0]);
      __m256d fz2 = masked2 ? _mm256_maskload_pd(&fz2ptrPerm[0], mask2) : _mm256_load_pd(&fz2ptrPerm[0]);

      __m256d mask1and2;
      bool masked = false;
      if (masked1 and masked2) {
        mask1and2 = _mm256_and_pd(_mm256_castsi256_pd(mask1), _mm256_castsi256_pd(mask2));
        masked = true;
      } else if (masked1) {
        mask1and2 = _mm256_castsi256_pd(mask1);
        masked = true;
      } else if (masked2) {
        mask1and2 = _mm256_castsi256_pd(mask2);
        masked = true;
      }

      if (masked) {
        if (newton3) {
          SoAKernelNew<true, true>(x1, y1, z1, x2, y2, z2, &fx1, &fy1, &fz1, &fx2, &fy2, &fz2,
                                    virialSumX, virialSumY, virialSumZ, upotSum, mask1and2);
        } else {
          SoAKernelNew<false, true>(x1, y1, z1, x2, y2, z2, &fx1, &fy1, &fz1, &fx2, &fy2, &fz2,
                                    virialSumX, virialSumY, virialSumZ, upotSum, mask1and2);
        }
      } else {
        if (newton3) {
          SoAKernelNew<true, false>(x1, y1, z1, x2, y2, z2, &fx1, &fy1, &fz1, &fx2, &fy2, &fz2,
                                    virialSumX, virialSumY, virialSumZ, upotSum, mask1and2);
        } else {
          SoAKernelNew<false, false>(x1, y1, z1, x2, y2, z2, &fx1, &fy1, &fz1, &fx2, &fy2, &fz2,
                                    virialSumX, virialSumY, virialSumZ, upotSum, mask1and2);
        }
      }
      

      double fx2resultPtrPerm [4];
      double fy2resultPtrPerm [4];
      double fz2resultPtrPerm [4];

      masked2 ? _mm256_maskstore_pd(&fx2resultPtrPerm[0], mask2, fx2) : _mm256_store_pd(&fx2resultPtrPerm[0], fx2);
      masked2 ? _mm256_maskstore_pd(&fy2resultPtrPerm[0], mask2, fy2) : _mm256_store_pd(&fy2resultPtrPerm[0], fy2);
      masked2 ? _mm256_maskstore_pd(&fz2resultPtrPerm[0], mask2, fz2) : _mm256_store_pd(&fz2resultPtrPerm[0], fz2);

      if (!masked2){
        invpermute(permArray[a], &fx2resultPtrPerm[0], &fx2ptr[0]);
        invpermute(permArray[a], &fy2resultPtrPerm[0], &fy2ptr[0]);
        invpermute(permArray[a], &fz2resultPtrPerm[0], &fz2ptr[0]);
      } else {
        invpermute(permArray[a], &fx2resultPtrPerm[0], &fx2ptr[0], r2);
        invpermute(permArray[a], &fy2resultPtrPerm[0], &fy2ptr[0], r2);
        invpermute(permArray[a], &fz2resultPtrPerm[0], &fz2ptr[0], r2);
      }

    }
      masked1 ? _mm256_maskstore_pd(&fx1ptr[0], mask1, fx1) : _mm256_store_pd(&fx1ptr[0], fx1);
      masked1 ? _mm256_maskstore_pd(&fy1ptr[0], mask1, fy1) : _mm256_store_pd(&fy1ptr[0], fy1);
      masked1 ? _mm256_maskstore_pd(&fz1ptr[0], mask1, fz1) : _mm256_store_pd(&fz1ptr[0], fz1);
  }

 public:
  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, bool newton3)
   */
  void SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, const bool newton3) override {
#ifdef __AVX__
    if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

    double *const __restrict__ x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    double *const __restrict__ x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    double *const __restrict__ y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    double *const __restrict__ z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    double *const __restrict__ fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    double *const __restrict__ fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    double *const __restrict__ fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    double *const __restrict__ fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    double *const __restrict__ fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    double *const __restrict__ fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d upotSum = _mm256_setzero_pd();

    // Copied from pairwiseFunctor/LJFunctor.h:264-281
    bool isHaloCell1 = false;
    bool isHaloCell2 = false;
    // Checks whether the cells are halo cells.
    // This check cannot be done if _lowCorner and _highCorner are not set. So we do this only if calculateGlobals is
    // defined. (as of 23.11.2018)
    if (calculateGlobals) {
      isHaloCell1 |= x1ptr[0] < _lowCorner[0] || x1ptr[0] >= _highCorner[0];
      isHaloCell1 |= y1ptr[0] < _lowCorner[1] || y1ptr[0] >= _highCorner[1];
      isHaloCell1 |= z1ptr[0] < _lowCorner[2] || z1ptr[0] >= _highCorner[2];
      isHaloCell2 |= x2ptr[0] < _lowCorner[0] || x2ptr[0] >= _highCorner[0];
      isHaloCell2 |= y2ptr[0] < _lowCorner[1] || y2ptr[0] >= _highCorner[1];
      isHaloCell2 |= z2ptr[0] < _lowCorner[2] || z2ptr[0] >= _highCorner[2];

      // This if is commented out because the AoS vs SoA test would fail otherwise. Even though it is physically
      // correct!
      /*if(_duplicatedCalculations and isHaloCell1 and isHaloCell2){
        return;
      }*/
    }

    unsigned int i = 0;
    for (; i < (soa1.getNumParticles() & ~(vecLength - 1)); i += 4) {
      unsigned int j = 0;
      for (; j < (soa2.getNumParticles() & ~(vecLength - 1)); j += 4) {

          //std::cout << "########## false, false " << i << " " << j << std::endl;
          SoAPermKernel<false, false>(&x1ptr[i], &y1ptr[i], &z1ptr[i], &fx1ptr[i], &fy1ptr[i], &fz1ptr[i], 0,
                               &x2ptr[j], &y2ptr[j], &z2ptr[j], &fx2ptr[j], &fy2ptr[j], &fz2ptr[j], 0,
                               &virialSumX, &virialSumY, &virialSumZ, &upotSum, newton3);
        }
        const int rest2 = (int)(soa2.getNumParticles() & (vecLength - 1));
        if (rest2 > 0) {
          //std::cout << "########## false, true " << i << std::endl;
          SoAPermKernel<false, true>(&x1ptr[i], &y1ptr[i], &z1ptr[i], &fx1ptr[i], &fy1ptr[i], &fz1ptr[i], 0,
                               &x2ptr[j], &y2ptr[j], &z2ptr[j], &fx2ptr[j], &fy2ptr[j], &fz2ptr[j], rest2,
                               &virialSumX, &virialSumY, &virialSumZ, &upotSum, newton3);

        }
      }

      const int rest1 = (int)(soa1.getNumParticles() & (vecLength - 1));
      if (rest1 > 0) {
        
        unsigned int j = 0;
        for (; j < (soa2.getNumParticles() & ~(vecLength - 1)); j += 4) {

            //std::cout << "########## true, false " << j << std::endl;
            SoAPermKernel<true, false>(&x1ptr[i], &y1ptr[i], &z1ptr[i], &fx1ptr[i], &fy1ptr[i], &fz1ptr[i], rest1,
                                 &x2ptr[j], &y2ptr[j], &z2ptr[j], &fx2ptr[j], &fy2ptr[j], &fz2ptr[j], 0,
                                 &virialSumX, &virialSumY, &virialSumZ, &upotSum, newton3);
            }
          const int rest2 = (int)(soa2.getNumParticles() & (vecLength - 1));
          if (rest2 > 0) {
            //std::cout << "########## true, true " << std::endl;
            SoAPermKernel<true, true>(&x1ptr[i], &y1ptr[i], &z1ptr[i], &fx1ptr[i], &fy1ptr[i], &fz1ptr[i], rest1,
                                 &x2ptr[j], &y2ptr[j], &z2ptr[j], &fx2ptr[j], &fy2ptr[j], &fz2ptr[j], rest2,
                                 &virialSumX, &virialSumY, &virialSumZ, &upotSum, newton3);

        }
      }

/*
 *    for (unsigned int i = 0; i < soa1.getNumParticles(); i ++) {
 *      __m256d fxacc = _mm256_setzero_pd();
 *      __m256d fyacc = _mm256_setzero_pd();
 *      __m256d fzacc = _mm256_setzero_pd();
 *
 *      const __m256d x1 = _mm256_broadcast_sd(&x1ptr[i]);
 *      const __m256d y1 = _mm256_broadcast_sd(&y1ptr[i]);
 *      const __m256d z1 = _mm256_broadcast_sd(&z1ptr[i]);
 *
 *      // floor soa2 numParticles to multiple of vecLength
 *      // soon to be deprecated
 *      if (newton3) {
 *        unsigned int j = 0;
 *        for (; j < (soa2.getNumParticles() & ~(vecLength - 1)); j += 4) {
 *          SoAKernel<true, false>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, fxacc, fyacc, fzacc,
 *                                 &virialSumX, &virialSumY, &virialSumZ, &upotSum);
 *        }
 *        const int rest = (int)(soa2.getNumParticles() & (vecLength - 1));
 *        if (rest > 0)
 *          SoAKernel<true, true>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, fxacc, fyacc, fzacc,
 *                                &virialSumX, &virialSumY, &virialSumZ, &upotSum, rest);
 *      } else {
 *        unsigned int j = 0;
 *        for (; j < (soa2.getNumParticles() & ~(vecLength - 1)); j += 4) {
 *          SoAKernel<false, false>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, fxacc, fyacc, fzacc,
 *                                  &virialSumX, &virialSumY, &virialSumZ, &upotSum);
 *        }
 *        const int rest = (int)(soa2.getNumParticles() & (vecLength - 1));
 *        if (rest > 0)
 *          SoAKernel<false, true>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, fxacc, fyacc, fzacc,
 *                                 &virialSumX, &virialSumY, &virialSumZ, &upotSum, rest);
 *      }
 *
 *      // horizontally reduce fDacc to sumfD
 *      // deprecated
 *      const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
 *      const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);
 *
 *      const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
 *      const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);
 *
 *      const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
 *      const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);
 *
 *      const union {
 *        __m128d reg;
 *        double arr[2];
 *      } sumfxfyVEC = {.reg = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh)};
 *      const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);
 *
 *      const double sumfx = sumfxfyVEC.arr[0];
 *      const double sumfy = sumfxfyVEC.arr[1];
 *      const double sumfz = _mm_cvtsd_f64(sumfzVEC);
 *
 *      fx1ptr[i] += sumfx;
 *      fy1ptr[i] += sumfy;
 *      fz1ptr[i] += sumfz;
 *    }
 */

    if (calculateGlobals) {
      const int threadnum = autopas_get_thread_num();

      // horizontally reduce virialSumX and virialSumY
      const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
      const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
      const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
      const __m128d hSumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      // horizontally reduce virialSumZ and upotSum
      const __m256d hSumVirialzUpot = _mm256_hadd_pd(virialSumZ, upotSum);
      const __m128d hSumVirialzUpotLow = _mm256_extractf128_pd(hSumVirialzUpot, 0);
      const __m128d hSumVirialzUpotHigh = _mm256_extractf128_pd(hSumVirialzUpot, 1);
      const __m128d hSumVirialzUpotVec = _mm_add_pd(hSumVirialzUpotHigh, hSumVirialzUpotLow);

      // globals = {virialX, virialY, virialZ, uPot}
      double globals[4];
      _mm_store_pd(&globals[0], hSumVirialxyVec);
      _mm_store_pd(&globals[2], hSumVirialzUpotVec);

      double energyfactor = 1.;
      if (_duplicatedCalculations) {
        // if we have duplicated calculations, i.e., we calculate interactions multiple times, we have to take care
        // that we do not add the energy multiple times!
        energyfactor = isHaloCell1 ? 0. : 1.;
        if (newton3) {
          energyfactor += isHaloCell2 ? 0. : 1.;
          energyfactor *= 0.5;  // we count the energies partly to one of the two cells!
        }
      }

      _aosThreadData[threadnum].virialSum[0] += globals[0] * energyfactor;
      _aosThreadData[threadnum].virialSum[1] += globals[1] * energyfactor;
      _aosThreadData[threadnum].virialSum[2] += globals[2] * energyfactor;
      _aosThreadData[threadnum].upotSum += globals[3] * energyfactor;
    }
#endif
  }

  // clang-format off
  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType> &soa, const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom, size_t iTo, bool newton3)
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
   */
  // clang-format on
  void SoAFunctor(SoA<SoAArraysType> &soa,
                  const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom,
                  size_t iTo, bool newton3) override {
    utils::ExceptionHandler::exception("Verlet SoA functor not implemented!");
  }

  /**
   * SoALoader
   * @param cell
   * @param soa
   * @param offset
   */
  AUTOPAS_FUNCTOR_SOALOADER(
      cell, soa, offset,
      // @todo it is probably better to resize the soa only once, before calling
      // SoALoader (verlet-list only)
      soa.resizeArrays(offset + cell.numParticles());

      if (cell.numParticles() == 0) return;

      unsigned long *const __restrict__ idptr = soa.template begin<Particle::AttributeNames::id>();
      double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();
      double *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
      double *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
      double *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

      auto cellIter = cell.begin();
      // load particles in SoAs
      for (size_t i = offset; cellIter.isValid(); ++cellIter, ++i) {
        idptr[i] = cellIter->getID();
        xptr[i] = cellIter->getR()[0];
        yptr[i] = cellIter->getR()[1];
        zptr[i] = cellIter->getR()[2];
        fxptr[i] = cellIter->getF()[0];
        fyptr[i] = cellIter->getF()[1];
        fzptr[i] = cellIter->getF()[2];
      })
  /**
   * soaextractor
   * @param cell
   * @param soa
   * @param offset
   */
  AUTOPAS_FUNCTOR_SOAEXTRACTOR(
      cell, soa, offset,
      // body start
      if (soa.getNumParticles() == 0) return;

      auto cellIter = cell.begin();

#ifndef NDEBUG
      unsigned long *const __restrict__ idptr = soa.template begin<Particle::AttributeNames::id>();
#endif

      double *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
      double *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
      double *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

      for (size_t i = offset; cellIter.isValid(); ++i, ++cellIter) {
        assert(idptr[i] == cellIter->getID());
        cellIter->setF({fxptr[i], fyptr[i], fzptr[i]});
      })

  /**
   * Get the number of flops used per kernel call. This should count the
   * floating point operations needed for two particles that lie within a cutoff
   * radius.
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall() {
    // Kernel: 12 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply
    // scale) sum Forces: 6 (forces) kernel total = 12 + 6 = 18
    return 18ul;
  }

  /**
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   */
  void initTraversal() override {
    _upotSum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _aosThreadData[i].setZero();
    }
  }

  /**
   * Accumulates global values, e.g. upot and virial.
   * @param newton3
   */
  void endTraversal(bool newton3) override {
    if (_postProcessed) {
      throw utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _upotSum += _aosThreadData[i].upotSum;
      _virialSum = ArrayMath::add(_virialSum, _aosThreadData[i].virialSum);
    }
    if (not newton3) {
      // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2 here.
      _upotSum *= 0.5;
      _virialSum = ArrayMath::mulScalar(_virialSum, 0.5);
    }
    // we have always calculated 6*upot, so we divide by 6 here!
    _upotSum /= 6.;
    _postProcessed = true;
  }

  /**
   * Get the potential Energy
   * @return the potential Energy
   */
  double getUpot() {
    if (not calculateGlobals) {
      throw utils::ExceptionHandler::AutoPasException(
          "Trying to get upot even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw utils::ExceptionHandler::AutoPasException("Cannot get upot, because endTraversal was not called.");
    }
    return _upotSum;
  }

  /**
   * Get the virial
   * @return the virial
   */
  double getVirial() {
    if (not calculateGlobals) {
      throw utils::ExceptionHandler::AutoPasException(
          "Trying to get virial even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw utils::ExceptionHandler::AutoPasException("Cannot get virial, because endTraversal was not called.");
    }
    return _virialSum[0] + _virialSum[1] + _virialSum[2];
  }

 private:
  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, upotSum{0.} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      upotSum = 0.;
    }

    // variables
    std::array<double, 3> virialSum;
    double upotSum;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[4];
  };
  // make sure of the size of AoSThreadData
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

  //  double _cutoffsquare, _epsilon24, _sigmasquare, _shift6;
  const __m256d _one;
  //  const std::vector<__m256i> _masks;
  const __m256i _masks[3];
  const __m256d _cutoffsquare;
  const __m256d _epsilon24;
  const __m256d _sigmasquare;
  double _shift6;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _upotSum;
  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // bool that defines whether duplicate calculations are happening
  bool _duplicatedCalculations;
  // lower and upper corner of the domain of the current process
  std::array<double, 3> _lowCorner, _highCorner;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  // number of double values that fit into a vector register.
  constexpr static size_t vecLength = 4;

};  // namespace autopas
}  // namespace autopas
