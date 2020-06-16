/**
 * @file LJFunctorAVX.h
 *
 * @date 17 Jan 2018
 * @author F. Gratl
 */
#pragma once

#include <immintrin.h>

#include <array>

#include "ParticlePropertiesLibrary.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticSelectorMacros.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * This Version is implemented using AVX intrinsics.
 * @tparam Particle The type of particle.
 * @tparam ParticleCell The type of particlecell.
 * @tparam applyShift Switch for the lj potential to be truncated shifted.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, class ParticleCell, bool applyShift = false, bool useMixing = false,
          FunctorN3Modes useNewton3 = FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true>
class LJFunctorAVX : public Functor<Particle, ParticleCell, typename Particle::SoAArraysType,
                                    LJFunctorAVX<Particle, ParticleCell, applyShift, useMixing, useNewton3,
                                                 calculateGlobals, relevantForTuning>> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorAVX() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @param dummy unused, only there to make the signature different from the public constructor.
   */
  explicit LJFunctorAVX(double cutoff, void * /*dummy*/)
#ifdef __AVX__
      : Functor<Particle, ParticleCell, SoAArraysType,
                LJFunctorAVX<Particle, ParticleCell, applyShift, useMixing, useNewton3, calculateGlobals,
                             relevantForTuning>>(cutoff),
        _cutoffsquare{_mm256_set1_pd(cutoff * cutoff)},
        _upotSum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas_get_max_threads());
    }
  }
#else
      : Functor<Particle, ParticleCell, SoAArraysType,
                LJFunctorAVX<Particle, ParticleCell, applyShift, useMixing, useNewton3, calculateGlobals,
                             relevantForTuning>>(cutoff) {
    utils::ExceptionHandler::exception("AutoPas was compiled without AVX support!");
  }
#endif
 public:
  /**
   * Constructor, which sets the global values, i.e. cutoff, epsilon, sigma and shift.
   *
   * @note Only to be used with mixing == false.
   *
   * @param cutoff
   */
  explicit LJFunctorAVX(double cutoff) : LJFunctorAVX(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor, which sets the global values, i.e. cutoff, epsilon, sigma and shift.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit LJFunctorAVX(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctorAVX(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  bool isRelevantForTuning() override { return relevantForTuning; }

  bool allowsNewton3() override {
    return useNewton3 == FunctorN3Modes::Newton3Only or useNewton3 == FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() override {
    return useNewton3 == FunctorN3Modes::Newton3Off or useNewton3 == FunctorN3Modes::Both;
  }

  bool isAppropriateClusterSize(unsigned int clusterSize, DataLayoutOption::Value dataLayout) const override {
    return dataLayout == DataLayoutOption::aos;  // LJFunctorAVX does only support clusters via aos.
  }

  void AoSFunctor(Particle &i, Particle &j, bool newton3) override {
    utils::ExceptionHandler::exception("LJFunctorAVX.AoSFunctor() not implemented!");
  }

  /**
   * @copydoc Functor::SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3)
   * This functor ignores the newton3 value, as we do not expect any benefit from disabling newton3.
   */
  void SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3) override {
    // using nested withStaticBool is not possible because of bug in gcc < 9 (and the intel compiler)
    /// @todo c++20: gcc < 9 can probably be dropped, replace with nested lambdas.
    utils::withStaticBool(newton3, [&](auto newton3) { SoAFunctorSingleImpl<newton3>(soa); });
  }

  // clang-format off
  /**
   * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
   */
  // clang-format on
  void SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, const bool newton3) override {
    utils::withStaticBool(newton3, [&](auto newton3) { SoAFunctorPairImpl<newton3>(soa1, soa2); });
  }

 private:
  /**
   * Templatized version of SoAFunctorSingle actually doing what the latter should.
   * @tparam newton3
   * @param soa
   */
  template <bool newton3>
  void SoAFunctorSingleImpl(SoAView<SoAArraysType> soa) {
#ifdef __AVX__
    if (soa.getNumParticles() == 0) return;

    const auto *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict__ ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict__ fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict__ fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict__ fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict__ typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d upotSum = _mm256_setzero_pd();

    // reverse outer loop s.th. inner loop always beginns at aligned array start
    // typecast to detect underflow
    for (size_t i = soa.getNumParticles() - 1; (long)i >= 0; --i) {
      if (ownedStatePtr[i] == OwnershipState::dummy) {
        // If the i-th particle is a dummy, skip this loop iteration.
        continue;
      }

      static_assert(std::is_same_v<std::underlying_type_t<OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");
      // ownedStatePtr contains int64_t, so we broadcast these to make an __m256i.
      // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
      __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr[i]));

      __m256d fxacc = _mm256_setzero_pd();
      __m256d fyacc = _mm256_setzero_pd();
      __m256d fzacc = _mm256_setzero_pd();

      const __m256d x1 = _mm256_broadcast_sd(&xptr[i]);
      const __m256d y1 = _mm256_broadcast_sd(&yptr[i]);
      const __m256d z1 = _mm256_broadcast_sd(&zptr[i]);

      // floor soa numParticles to multiple of vecLength
      size_t j = 0;
      for (; j < (i & ~(vecLength - 1)); j += 4) {
        SoAKernel<true, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr,
                               zptr, fxptr, fyptr, fzptr, typeIDptr, typeIDptr, fxacc, fyacc, fzacc, &virialSumX,
                               &virialSumY, &virialSumZ, &upotSum, 0);
      }
      const int rest = (int)(i & (vecLength - 1));
      if (rest > 0) {
        SoAKernel<true, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr,
                              zptr, fxptr, fyptr, fzptr, typeIDptr, typeIDptr, fxacc, fyacc, fzacc, &virialSumX,
                              &virialSumY, &virialSumZ, &upotSum, rest);
      }

      // horizontally reduce fDacc to sumfD
      const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
      const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

      const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
      const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

      const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
      const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

      const __m128d sumfxfyVEC = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
      const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

      const double sumfx = sumfxfyVEC[0];
      const double sumfy = sumfxfyVEC[1];
      const double sumfz = _mm_cvtsd_f64(sumfzVEC);

      fxptr[i] += sumfx;
      fyptr[i] += sumfy;
      fzptr[i] += sumfz;
    }

    if constexpr (calculateGlobals) {
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

      double factor = 1.;
      // we assume newton3 to be enabled in this functor call, thus we multiply by two if the value of newton3 is false,
      // since for newton3 disabled we divide by two later on.
      factor *= newton3 ? .5 : 1.;
      // In case we have a non-cell-wise owned state, we have multiplied everything by two, so we divide it by 2 again.
      _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
      _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
      _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
      _aosThreadData[threadnum].upotSum += globals[3] * factor;
    }
#endif
  }

  template <bool newton3>
  void SoAFunctorPairImpl(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2) {
#ifdef __AVX__
    if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

    const auto *const __restrict__ x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict__ y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict__ z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict__ x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict__ y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict__ z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict__ ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict__ ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

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

    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      if (ownedStatePtr1[i] == OwnershipState::dummy) {
        // If the i-th particle is a dummy, skip this loop iteration.
        continue;
      }

      __m256d fxacc = _mm256_setzero_pd();
      __m256d fyacc = _mm256_setzero_pd();
      __m256d fzacc = _mm256_setzero_pd();

      static_assert(std::is_same_v<std::underlying_type_t<OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");
      // ownedStatePtr1 contains int64_t, so we broadcast these to make an __m256i.
      // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
      __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr1[i]));

      const __m256d x1 = _mm256_broadcast_sd(&x1ptr[i]);
      const __m256d y1 = _mm256_broadcast_sd(&y1ptr[i]);
      const __m256d z1 = _mm256_broadcast_sd(&z1ptr[i]);

      // floor soa2 numParticles to multiple of vecLength
      unsigned int j = 0;
      for (; j < (soa2.getNumParticles() & ~(vecLength - 1)); j += 4) {
        SoAKernel<newton3, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                  y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc,
                                  &virialSumX, &virialSumY, &virialSumZ, &upotSum, 0);
      }
      const int rest = (int)(soa2.getNumParticles() & (vecLength - 1));
      if (rest > 0)
        SoAKernel<newton3, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                 y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc,
                                 &virialSumX, &virialSumY, &virialSumZ, &upotSum, rest);

      // horizontally reduce fDacc to sumfD
      const __m256d hSumfxfy = _mm256_hadd_pd(fxacc, fyacc);
      const __m256d hSumfz = _mm256_hadd_pd(fzacc, fzacc);

      const __m128d hSumfxfyLow = _mm256_extractf128_pd(hSumfxfy, 0);
      const __m128d hSumfzLow = _mm256_extractf128_pd(hSumfz, 0);

      const __m128d hSumfxfyHigh = _mm256_extractf128_pd(hSumfxfy, 1);
      const __m128d hSumfzHigh = _mm256_extractf128_pd(hSumfz, 1);

      const __m128d sumfxfyVEC = _mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
      const __m128d sumfzVEC = _mm_add_pd(hSumfzLow, hSumfzHigh);

      const double sumfx = sumfxfyVEC[0];
      const double sumfy = sumfxfyVEC[1];
      const double sumfz = _mm_cvtsd_f64(sumfzVEC);

      fx1ptr[i] += sumfx;
      fy1ptr[i] += sumfy;
      fz1ptr[i] += sumfz;
    }

    if constexpr (calculateGlobals) {
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

      // we have duplicated calculations, i.e., we calculate interactions multiple times, so we have to take care
      // that we do not add the energy multiple times!
      double energyfactor = 1.;
      if constexpr (newton3) {
        energyfactor *= 0.5;  // we count the energies partly to one of the two cells!
      }

      _aosThreadData[threadnum].virialSum[0] += globals[0] * energyfactor;
      _aosThreadData[threadnum].virialSum[1] += globals[1] * energyfactor;
      _aosThreadData[threadnum].virialSum[2] += globals[2] * energyfactor;
      _aosThreadData[threadnum].upotSum += globals[3] * energyfactor;
    }
#endif
  }

  /**
   * Actual inner kernel of the SoAFunctors.
   *
   * @tparam newton3
   * @tparam remainderIsMasked If false the full vector length is used. Otherwise the last entries are masked away
   * depending on the argument "rest".
   * @param j
   * @param x1
   * @param y1
   * @param z1
   * @param x2ptr
   * @param y2ptr
   * @param z2ptr
   * @param fx2ptr
   * @param fy2ptr
   * @param fz2ptr
   * @param typeID1ptr
   * @param typeID2ptr
   * @param fxacc
   * @param fyacc
   * @param fzacc
   * @param virialSumX
   * @param virialSumY
   * @param virialSumZ
   * @param upotSum
   * @param rest
   */
  template <bool newton3, bool remainderIsMasked>
  inline void SoAKernel(const size_t j, const __m256i ownedStateI, const int64_t *const __restrict__ ownedStatePtr2,
                        const __m256d &x1, const __m256d &y1, const __m256d &z1, const double *const __restrict__ x2ptr,
                        const double *const __restrict__ y2ptr, const double *const __restrict__ z2ptr,
                        double *const __restrict__ fx2ptr, double *const __restrict__ fy2ptr,
                        double *const __restrict__ fz2ptr, const size_t *const typeID1ptr,
                        const size_t *const typeID2ptr, __m256d &fxacc, __m256d &fyacc, __m256d &fzacc,
                        __m256d *virialSumX, __m256d *virialSumY, __m256d *virialSumZ, __m256d *upotSum,
                        const unsigned int rest = 0) {
#ifdef __AVX__
    __m256d epsilon24s = _epsilon24;
    __m256d sigmaSquares = _sigmaSquare;
    __m256d shift6s = _shift6;
    if (useMixing) {
      // the first argument for set lands in the last bits of the register
      epsilon24s = _mm256_set_pd(
          not remainderIsMasked or rest > 3 ? _PPLibrary->mixing24Epsilon(*typeID1ptr, *(typeID2ptr + 3)) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->mixing24Epsilon(*typeID1ptr, *(typeID2ptr + 2)) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->mixing24Epsilon(*typeID1ptr, *(typeID2ptr + 1)) : 0,
          _PPLibrary->mixing24Epsilon(*typeID1ptr, *(typeID2ptr + 0)));
      sigmaSquares = _mm256_set_pd(
          not remainderIsMasked or rest > 3 ? _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(typeID2ptr + 3)) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(typeID2ptr + 2)) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(typeID2ptr + 1)) : 0,
          _PPLibrary->mixingSigmaSquare(*typeID1ptr, *(typeID2ptr + 0)));
      if constexpr (applyShift) {
        shift6s = _mm256_set_pd(
            (not remainderIsMasked or rest > 3) ? _PPLibrary->mixingShift6(*typeID1ptr, *(typeID2ptr + 3)) : 0,
            (not remainderIsMasked or rest > 2) ? _PPLibrary->mixingShift6(*typeID1ptr, *(typeID2ptr + 2)) : 0,
            (not remainderIsMasked or rest > 1) ? _PPLibrary->mixingShift6(*typeID1ptr, *(typeID2ptr + 1)) : 0,
            _PPLibrary->mixingShift6(*typeID1ptr, *(typeID2ptr + 0)));
      }
    }

    const __m256d x2 = remainderIsMasked ? _mm256_maskload_pd(&x2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&x2ptr[j]);
    const __m256d y2 = remainderIsMasked ? _mm256_maskload_pd(&y2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&y2ptr[j]);
    const __m256d z2 = remainderIsMasked ? _mm256_maskload_pd(&z2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&z2ptr[j]);

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

    // This requires that dummy is zero (otherwise when loading using a mask the owned state will not be zero)
    const __m256i ownedStateJ =
        remainderIsMasked
            ? _mm256_maskload_epi64(reinterpret_cast<long long const *>(&ownedStatePtr2[j]), _masks[rest - 1])
            : _mm256_load_si256(reinterpret_cast<const __m256i *>(&ownedStatePtr2[j]));
    // This requires that dummy is the first entry in OwnershipState!
    const __m256i dummyMask = _mm256_cmpgt_epi64(ownedStateJ, _ownedStateDummyMM256i);
    const __m256d cutoffDummyMask = _mm256_and_pd(cutoffMask, _mm256_castsi256_pd(dummyMask));

    // if everything is masked away return from this function.
    if (_mm256_movemask_pd(cutoffDummyMask) == 0) {
      return;
    }

    const __m256d invdr2 = _mm256_div_pd(_one, dr2);
    const __m256d lj2 = _mm256_mul_pd(sigmaSquares, invdr2);
    const __m256d lj4 = _mm256_mul_pd(lj2, lj2);
    const __m256d lj6 = _mm256_mul_pd(lj2, lj4);
    const __m256d lj12 = _mm256_mul_pd(lj6, lj6);
    const __m256d lj12m6 = _mm256_sub_pd(lj12, lj6);
    const __m256d lj12m6alj12 = _mm256_add_pd(lj12m6, lj12);
    const __m256d lj12m6alj12e = _mm256_mul_pd(lj12m6alj12, epsilon24s);
    const __m256d fac = _mm256_mul_pd(lj12m6alj12e, invdr2);

    const __m256d facMasked =
        remainderIsMasked ? _mm256_and_pd(fac, _mm256_and_pd(cutoffDummyMask, _mm256_castsi256_pd(_masks[rest - 1])))
                          : _mm256_and_pd(fac, cutoffDummyMask);

    const __m256d fx = _mm256_mul_pd(drx, facMasked);
    const __m256d fy = _mm256_mul_pd(dry, facMasked);
    const __m256d fz = _mm256_mul_pd(drz, facMasked);

    fxacc = _mm256_add_pd(fxacc, fx);
    fyacc = _mm256_add_pd(fyacc, fy);
    fzacc = _mm256_add_pd(fzacc, fz);

    // if newton 3 is used subtract fD from particle j
    if (newton3) {
      const __m256d fx2 =
          remainderIsMasked ? _mm256_maskload_pd(&fx2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&fx2ptr[j]);
      const __m256d fy2 =
          remainderIsMasked ? _mm256_maskload_pd(&fy2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&fy2ptr[j]);
      const __m256d fz2 =
          remainderIsMasked ? _mm256_maskload_pd(&fz2ptr[j], _masks[rest - 1]) : _mm256_load_pd(&fz2ptr[j]);

      const __m256d fx2new = _mm256_sub_pd(fx2, fx);
      const __m256d fy2new = _mm256_sub_pd(fy2, fy);
      const __m256d fz2new = _mm256_sub_pd(fz2, fz);

      remainderIsMasked ? _mm256_maskstore_pd(&fx2ptr[j], _masks[rest - 1], fx2new)
                        : _mm256_store_pd(&fx2ptr[j], fx2new);
      remainderIsMasked ? _mm256_maskstore_pd(&fy2ptr[j], _masks[rest - 1], fy2new)
                        : _mm256_store_pd(&fy2ptr[j], fy2new);
      remainderIsMasked ? _mm256_maskstore_pd(&fz2ptr[j], _masks[rest - 1], fz2new)
                        : _mm256_store_pd(&fz2ptr[j], fz2new);
    }

    if (calculateGlobals) {
      // Global Virial
      const __m256d virialX = _mm256_mul_pd(fx, drx);
      const __m256d virialY = _mm256_mul_pd(fy, dry);
      const __m256d virialZ = _mm256_mul_pd(fz, drz);

      // Global Potential
      const __m256d upot = wrapperFMA(epsilon24s, lj12m6, shift6s);

      const __m256d upotMasked =
          remainderIsMasked ? _mm256_and_pd(upot, _mm256_and_pd(cutoffDummyMask, _mm256_castsi256_pd(_masks[rest - 1])))
                            : _mm256_and_pd(upot, cutoffDummyMask);

      __m256i ownedMaskI = _mm256_cmpeq_epi64(ownedStateI, _ownedStateOwnedMM256i);
      __m256d energyFactor = _mm256_blendv_pd(_zero, _one, _mm256_castsi256_pd(ownedMaskI));
      if constexpr (newton3) {
        __m256i ownedMaskJ = _mm256_cmpeq_epi64(ownedStateJ, _ownedStateOwnedMM256i);
        energyFactor = _mm256_add_pd(energyFactor, _mm256_blendv_pd(_zero, _one, _mm256_castsi256_pd(ownedMaskJ)));
      }
      *upotSum = wrapperFMA(energyFactor, upotMasked, *upotSum);
      *virialSumX = wrapperFMA(energyFactor, virialX, *virialSumX);
      *virialSumY = wrapperFMA(energyFactor, virialY, *virialSumY);
      *virialSumZ = wrapperFMA(energyFactor, virialZ, *virialSumZ);
    }
#endif
  }

 public:
  // clang-format off
  /**
   * @copydoc Functor::SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
   */
  // clang-format on
  void SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) override {
    utils::ExceptionHandler::exception("Verlet SoA functor not implemented!");
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 9>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 6>{
        Particle::AttributeNames::id,   Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 3>{
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ};
  }

  /**
   *
   * @return useMixing
   */
  constexpr static bool getMixing() { return useMixing; }

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

    if (calculateGlobals) {
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        _upotSum += _aosThreadData[i].upotSum;
        _virialSum = utils::ArrayMath::add(_virialSum, _aosThreadData[i].virialSum);
      }
      if (not newton3) {
        // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
        // here.
        _upotSum *= 0.5;
        _virialSum = utils::ArrayMath::mulScalar(_virialSum, 0.5);
      }
      // we have always calculated 6*upot, so we divide by 6 here!
      _upotSum /= 6.;
      _postProcessed = true;
    }
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
#ifdef __AVX__
    _epsilon24 = _mm256_set1_pd(epsilon24);
    _sigmaSquare = _mm256_set1_pd(sigmaSquare);
    if constexpr (applyShift) {
      _shift6 = _mm256_set1_pd(
          ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, _cutoffsquare[0]));
    } else {
      _shift6 = _mm256_setzero_pd();
    }
#endif
  }

 private:
  /**
   * Wrapper function for FMA. If FMA is not supported it executes first the multiplication then the addition.
   * @param factorA
   * @param factorB
   * @param summandC
   * @return A * B + C
   */
  inline __m256d wrapperFMA(const __m256d &factorA, const __m256d &factorB, const __m256d &summandC) {
#ifdef __FMA__
    return _mm256_fmadd_pd(factorA, factorB, summandC);
#elif __AVX__
    const __m256d tmp = _mm256_mul_pd(factorA, factorB);
    return _mm256_add_pd(summandC, tmp);
#else
    // dummy return. If no vectorization is available this whole class is pointless anyways.
    return __m256d();
#endif
  }

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

#ifdef __AVX__
  const __m256d _zero{_mm256_set1_pd(0.)};
  const __m256d _one{_mm256_set1_pd(1.)};
  const __m256i _masks[3]{
      _mm256_set_epi64x(0, 0, 0, -1),
      _mm256_set_epi64x(0, 0, -1, -1),
      _mm256_set_epi64x(0, -1, -1, -1),
  };
  const __m256i _ownedStateDummyMM256i{0x0};
  const __m256i _ownedStateOwnedMM256i{_mm256_set1_epi64x(static_cast<int64_t>(OwnershipState::owned))};
  const __m256d _cutoffsquare{};
  __m256d _shift6 = _mm256_setzero_pd();
  __m256d _epsilon24{};
  __m256d _sigmaSquare{};
#endif

  ParticlePropertiesLibrary<double, size_t> *_PPLibrary = nullptr;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _upotSum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  // number of double values that fit into a vector register.
  constexpr static size_t vecLength = 4;

};  // namespace autopas
}  // namespace autopas
