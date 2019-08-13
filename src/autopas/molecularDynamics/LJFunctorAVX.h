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
template <class Particle, class ParticleCell, bool useMixing = false, FunctorN3Modes useNewton3 = FunctorN3Modes::Both,
          bool calculateGlobals = false, bool relevantForTuning = true>
class LJFunctorAVX
    : public Functor<Particle, ParticleCell, typename Particle::SoAArraysType, LJFunctorAVX<Particle, ParticleCell>> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorAVX() = delete;

  /**
   * Constructor, which sets the global values, i.e. cutoff, epsilon, sigma and shift.
   * @param cutoff
   * @param shift
   * @param duplicatedCalculation Defines whether duplicated calculations are happening across processes / over the
   * simulation boundary. e.g. eightShell: false, fullShell: true.
   */
  explicit LJFunctorAVX(double cutoff, double shift, bool duplicatedCalculation = false)
#ifdef __AVX__
      : Functor<Particle, ParticleCell, SoAArraysType, LJFunctorAVX<Particle, ParticleCell>>(cutoff),
        _one{_mm256_set1_pd(1.)},
        _masks{
            _mm256_set_epi64x(0, 0, 0, -1),
            _mm256_set_epi64x(0, 0, -1, -1),
            _mm256_set_epi64x(0, -1, -1, -1),
        },
        _cutoffsquare{_mm256_set1_pd(cutoff * cutoff)},
        _shift6{_mm256_set1_pd(shift * 6.0)},
        _upotSum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _duplicatedCalculations{duplicatedCalculation},
        _postProcessed{false} {

    if (calculateGlobals) {
      _aosThreadData.resize(autopas_get_max_threads());
    }
  }
#else
      : Functor<Particle, ParticleCell, SoAArraysType, LJFunctorAVX<Particle, ParticleCell>>(cutoff),
        _one{0},
        _masks{0, 0, 0},
        _cutoffsquare{0},
        _shift6{0},
        _epsilon24{0},
        _sigmaSquare{0} {
    utils::ExceptionHandler::exception("AutoPas was compiled without AVX support!");
  }
#endif

  /**
   * Constructor, which sets the global values, i.e. cutoff, epsilon, sigma and shift.
   * @param cutoff
   * @param shift
   * @param particlePropertiesLibrary
   * @param duplicatedCalculation Defines whether duplicated calculations are happening across processes / over the
   * simulation boundary. e.g. eightShell: false, fullShell: true.
   */
  explicit LJFunctorAVX(double cutoff, double shift,
                        ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary,
                        bool duplicatedCalculation = true)
      : LJFunctorAVX(cutoff, shift, duplicatedCalculation) {
    _PPLibrary = &particlePropertiesLibrary;
    //@TODO: check that this constructor is only called with mixing == true
  }

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
   * @copydoc Functor::SoAFunctor(SoAView<SoAArraysType> soa, bool newton3)
   * This functor ignores the newton3 value, as we do not expect any benefit from disabling newton3.
   */
  void SoAFunctor(SoAView<SoAArraysType> soa, bool newton3) override {
#ifdef __AVX__
    if (soa.getNumParticles() == 0) return;

    const auto *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict__ ownedPtr = soa.template begin<Particle::AttributeNames::owned>();

    auto *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict__ typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d upotSum = _mm256_setzero_pd();

    if (calculateGlobals) {
      // Checks if the cell is a halo cell, if it is, we skip it.
      bool isHaloCell = not ownedPtr[0];
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
      size_t j = 0;
      for (; j < (i & ~(vecLength - 1)); j += 4) {
        SoAKernel<true, false>(j, x1, y1, z1, xptr, yptr, zptr, fxptr, fyptr, fzptr, typeIDptr, typeIDptr, fxacc, fyacc,
                               fzacc, &virialSumX, &virialSumY, &virialSumZ, &upotSum, 0);
      }
      const int rest = (int)(i & (vecLength - 1));
      if (rest > 0)
        SoAKernel<true, true>(j, x1, y1, z1, xptr, yptr, zptr, fxptr, fyptr, fzptr, typeIDptr, typeIDptr, fxacc, fyacc,
                              fzacc, &virialSumX, &virialSumY, &virialSumZ, &upotSum, rest);

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
  inline void SoAKernel(const size_t j, const __m256d &x1, const __m256d &y1, const __m256d &z1,
                        const double *const x2ptr, const double *const y2ptr, const double *const z2ptr,
                        double *const fx2ptr, double *const fy2ptr, double *const fz2ptr,
                        const size_t *const typeID1ptr, const size_t *const type2IDptr, __m256d &fxacc, __m256d &fyacc,
                        __m256d &fzacc, __m256d *virialSumX, __m256d *virialSumY, __m256d *virialSumZ, __m256d *upotSum,
                        const unsigned int rest = 0) {
#ifdef __AVX__
    __m256d epsilon24s = _epsilon24;
    __m256d sigmaSquares = _sigmaSquare;
    if (useMixing) {
      epsilon24s = _mm256_set_pd(_PPLibrary->mixing24Epsilon(*typeID1ptr, *(type2IDptr + 0)),
                                 _PPLibrary->mixing24Epsilon(*typeID1ptr, *(type2IDptr + 1)),
                                 _PPLibrary->mixing24Epsilon(*typeID1ptr, *(type2IDptr + 2)),
                                 _PPLibrary->mixing24Epsilon(*typeID1ptr, *(type2IDptr + 3)));
      sigmaSquares = _mm256_set_pd(_PPLibrary->mixingSigmaSquare(*typeID1ptr, *(type2IDptr + 0)),
                                   _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(type2IDptr + 1)),
                                   _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(type2IDptr + 2)),
                                   _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(type2IDptr + 3)));
    }

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
    // @TODO abort if everything is masked away

    const __m256d invdr2 = _mm256_div_pd(_one, dr2);
    const __m256d lj2 = _mm256_mul_pd(sigmaSquares, invdr2);
    const __m256d lj4 = _mm256_mul_pd(lj2, lj2);
    const __m256d lj6 = _mm256_mul_pd(lj2, lj4);
    const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
    const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);
    const __m256d lj12m6alj12 = _mm256_add_pd(lj12m6, lj12);
    const __m256d lj12m6alj12e = _mm256_mul_pd(lj12m6alj12, epsilon24s);
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
      const __m256d upotPART1 = _mm256_mul_pd(epsilon24s, lj12m6);
      const __m256d upotPART2 = _mm256_add_pd(_shift6, upotPART1);
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

 public:
  /**
   * @copydoc Functor::SoAFunctor(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
   */
  void SoAFunctor(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, const bool newton3) override {
#ifdef __AVX__
    if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

    const auto *const __restrict__ x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict__ y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict__ z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict__ x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict__ y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict__ z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict__ ownedPtr1 = soa1.template begin<Particle::AttributeNames::owned>();
    const auto *const __restrict__ ownedPtr2 = soa2.template begin<Particle::AttributeNames::owned>();

    auto *const __restrict__ fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict__ fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict__ fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict__ fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict__ fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict__ fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict__ typeID1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    const auto *const __restrict__ typeID2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

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
      isHaloCell1 = not ownedPtr1[0];
      isHaloCell2 = not ownedPtr2[0];

      // This if is commented out because the AoS vs SoA test would fail otherwise. Even though it is physically
      // correct!
      /*if(_duplicatedCalculations and isHaloCell1 and isHaloCell2){
        return;
      }*/
    }

    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      __m256d fxacc = _mm256_setzero_pd();
      __m256d fyacc = _mm256_setzero_pd();
      __m256d fzacc = _mm256_setzero_pd();

      const __m256d x1 = _mm256_broadcast_sd(&x1ptr[i]);
      const __m256d y1 = _mm256_broadcast_sd(&y1ptr[i]);
      const __m256d z1 = _mm256_broadcast_sd(&z1ptr[i]);

      // floor soa2 numParticles to multiple of vecLength
      if (newton3) {
        unsigned int j = 0;
        for (; j < (soa2.getNumParticles() & ~(vecLength - 1)); j += 4) {
          SoAKernel<true, false>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr,
                                 fxacc, fyacc, fzacc, &virialSumX, &virialSumY, &virialSumZ, &upotSum, 0);
        }
        const int rest = (int)(soa2.getNumParticles() & (vecLength - 1));
        if (rest > 0)
          SoAKernel<true, true>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr,
                                fxacc, fyacc, fzacc, &virialSumX, &virialSumY, &virialSumZ, &upotSum, rest);
      } else {
        unsigned int j = 0;
        for (; j < (soa2.getNumParticles() & ~(vecLength - 1)); j += 4) {
          SoAKernel<false, false>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr,
                                  fxacc, fyacc, fzacc, &virialSumX, &virialSumY, &virialSumZ, &upotSum, 0);
        }
        const int rest = (int)(soa2.getNumParticles() & (vecLength - 1));
        if (rest > 0)
          SoAKernel<false, true>(j, x1, y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr,
                                 fxacc, fyacc, fzacc, &virialSumX, &virialSumY, &virialSumZ, &upotSum, rest);
      }

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

      fx1ptr[i] += sumfx;
      fy1ptr[i] += sumfy;
      fz1ptr[i] += sumfz;
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
   * @copydoc Functor::SoAFunctor(SoAView<SoAArraysType> soa, const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom, size_t iTo, bool newton3)
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
   */
  // clang-format on
  void SoAFunctor(SoAView<SoAArraysType> soa,
                  const std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &neighborList, size_t iFrom,
                  size_t iTo, bool newton3) override {
    utils::ExceptionHandler::exception("Verlet SoA functor not implemented!");
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static const std::array<typename Particle::AttributeNames, 8> getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 8>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::owned};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static const std::array<typename Particle::AttributeNames, 5> getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 5>{
        Particle::AttributeNames::id, Particle::AttributeNames::posX, Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::owned};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static const std::array<typename Particle::AttributeNames, 3> getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ};
  }

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

  /**
 * Sets the particle properties constants for this functor.
 *
 * This is only necessary if no particlePropertiesLibrary is used.
 *
 * @param epsilon24
 * @param sigmaSquare
 */
  void setParticleProperties(double epsilon24, double sigmaSquare) {
    _epsilon24 = _mm256_set1_pd(epsilon24);
    _sigmaSquare = _mm256_set1_pd(sigmaSquare);
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

  const __m256d _one;
  const __m256i _masks[3];
  const __m256d _cutoffsquare;
  const __m256d _shift6;
  __m256d _epsilon24;
  __m256d _sigmaSquare;

  ParticlePropertiesLibrary<double, size_t> *_PPLibrary = nullptr;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _upotSum;
  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // bool that defines whether duplicate calculations are happening
  bool _duplicatedCalculations;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  // number of double values that fit into a vector register.
  constexpr static size_t vecLength = 4;

};  // namespace autopas
}  // namespace autopas
