/**
 * @file LJFunctorAVX.h
 *
 * @date 17 Jan 2018
 * @author F. Gratl
 */
#pragma once
#ifndef __AVX__
#pragma message "Requested to compile LJFunctorAVX but AVX is not available!"
#else
#include <immintrin.h>
#endif

#include <array>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

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
 * @tparam countFLOPs counts FLOPs and hitrate. Not implemented for this functor. Please use the AutoVec functor.
 */
template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctorAVX
    : public autopas::PairwiseFunctor<Particle, LJFunctorAVX<Particle, applyShift, useMixing, useNewton3,
                                                             calculateGlobals, countFLOPs, relevantForTuning>> {
  /**
   * FloatType used for calculations
   */
  using CalcType = typename Particle::ParticleCalcType;

  /**
   * FloatType used for accumulations or more relevant calculations
   */
  using AccuType = typename Particle::ParticleAccuType;

  using SIMDCalcType = typename std::conditional_t<std::is_same<CalcType, float>::value, __m256, __m256d>;
  using SIMDAccuType = typename std::conditional_t<std::is_same<AccuType, float>::value, __m256, __m256d>;

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
   * @note param dummy unused, only there to make the signature different from the public constructor.
   */
  explicit LJFunctorAVX(CalcType cutoff, void * /*dummy*/)
#ifdef __AVX__
      : autopas::PairwiseFunctor<Particle, LJFunctorAVX<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
                                                        countFLOPs, relevantForTuning>>(cutoff),
#if AUTOPAS_PRECISION_MODE == SPSP || AUTOPAS_PRECISION_MODE == SPDP
        _cutoffSquared{_mm256_set1_ps(cutoff * cutoff)},
#else
        _cutoffSquared{_mm256_set1_pd(cutoff * cutoff)},
#endif
        _cutoffSquaredAoS(cutoff * cutoff),
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      AutoPasLog(DEBUG, "Using LJFunctorAVX with countFLOPs but FLOP counting is not implemented.");
    }
  }
#else
      : autopas::Functor<
            Particle, LJFunctorAVX<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff) {
    autopas::utils::ExceptionHandler::exception("AutoPas was compiled without AVX support!");
  }
#endif
 public:
  /**
   * Constructor for Functor with mixing disabled. When using this functor it is necessary to call
   * setParticleProperties() to set internal constants because it does not use a particle properties library.
   *
   * @note Only to be used with mixing == false.
   *
   * @param cutoff
   */
  explicit LJFunctorAVX(CalcType cutoff) : LJFunctorAVX(cutoff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like sigma, epsilon and shift.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit LJFunctorAVX(CalcType cutoff, ParticlePropertiesLibrary<CalcType, size_t> &particlePropertiesLibrary)
      : LJFunctorAVX(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  std::string getName() final { return "LJFunctorAVX"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  inline void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy()) {
      return;
    }

    CalcType sigmaSquared = _sigmaSquaredAoS;
    CalcType epsilon24 = _epsilon24AoS;
    CalcType shift6 = _shift6AoS;

    if constexpr (useMixing) {
      sigmaSquared = _PPLibrary->getMixingSigmaSquared(i.getTypeId(), j.getTypeId());
      epsilon24 = _PPLibrary->getMixing24Epsilon(i.getTypeId(), j.getTypeId());
      if constexpr (applyShift) {
        shift6 = _PPLibrary->getMixingShift6(i.getTypeId(), j.getTypeId());
      }
    }

    std::array<CalcType, 3> dr = i.getR() - j.getR();
    CalcType dr2 = autopas::utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffSquaredAoS) {
      return;
    }

    CalcType invdr2 = static_cast<CalcType>(1.) / dr2;
    CalcType lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    CalcType lj12 = lj6 * lj6;
    CalcType lj12m6 = lj12 - lj6;
    CalcType fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    std::array<CalcType, 3> f = dr * fac;

    const std::array<AccuType, 3> convertedF = autopas::utils::ArrayUtils::static_cast_copy_array<AccuType>(f);

    i.addF(convertedF);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(convertedF);
    }

    if (calculateGlobals) {
      // We always add the full contribution for each owned particle and divide the sums by 2 in endTraversal().
      // Potential energy has an additional factor of 6, which is also handled in endTraversal().

      std::array<CalcType, 3> virial = dr * f;
      CalcType potentialEnergy6 = epsilon24 * lj12m6 + shift6;

      const int threadnum = autopas::autopas_get_thread_num();
      if (i.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadData[threadnum].virialSum += autopas::utils::ArrayUtils::static_cast_copy_array<AccuType>(virial);
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadData[threadnum].virialSum += autopas::utils::ArrayUtils::static_cast_copy_array<AccuType>(virial);
      }
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always do a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
   */
  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (newton3) {
      SoAFunctorSingleImpl<true>(soa);
    } else {
      SoAFunctorSingleImpl<false>(soa);
    }
  }

  // clang-format off
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
   */
  // clang-format on
  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                             const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

 private:
  /**
   * Templatized version of SoAFunctorSingle actually doing what the latter should.
   * @tparam newton3
   * @param soa
   */
  template <bool newton3>
  inline void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
#ifdef __AVX__
    if (soa.size() == 0) return;

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

#if AUTOPAS_PRECISION_MODE == SPSP
    __m256 virialSumX = _mm256_setzero_ps();
    __m256 virialSumY = _mm256_setzero_ps();
    __m256 virialSumZ = _mm256_setzero_ps();
    __m256 potentialEnergySum = _mm256_setzero_ps();
#else
    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d potentialEnergySum = _mm256_setzero_pd();
#endif

    // reverse outer loop s.th. inner loop always beginns at aligned array start
    // typecast to detect underflow
    for (size_t i = soa.size() - 1; (long)i >= 0; --i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        // If the i-th particle is a dummy, skip this loop iteration.
        continue;
      }

// values for calculations
#if AUTOPAS_PRECISION_MODE == SPSP || AUTOPAS_PRECISION_MODE == SPDP
      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int32_t>,
                    "OwnershipStates underlying type should be int32_t!");
      // ownedStatePtr contains int32_t, so we broadcast these to make an __m256i.
      // _mm256_set1_epi32 broadcasts a 32-bit integer, we use this instruction to have 8 values!
      __m256i ownedStateI = _mm256_set1_epi32(static_cast<int32_t>(ownedStatePtr[i]));

      const __m256 x1 = _mm256_broadcast_ss(&xptr[i]);
      const __m256 y1 = _mm256_broadcast_ss(&yptr[i]);
      const __m256 z1 = _mm256_broadcast_ss(&zptr[i]);
#else
      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");
      // ownedStatePtr contains int64_t, so we broadcast these to make an __m256i.
      // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
      __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr[i]));

      const __m256d x1 = _mm256_broadcast_sd(&xptr[i]);
      const __m256d y1 = _mm256_broadcast_sd(&yptr[i]);
      const __m256d z1 = _mm256_broadcast_sd(&zptr[i]);
#endif

// these are about accumulation, so different check
#if AUTOPAS_PRECISION_MODE == SPSP
      __m256 fxacc = _mm256_setzero_ps();
      __m256 fyacc = _mm256_setzero_ps();
      __m256 fzacc = _mm256_setzero_ps();
#else
      __m256d fxacc = _mm256_setzero_pd();
      __m256d fyacc = _mm256_setzero_pd();
      __m256d fzacc = _mm256_setzero_pd();
#endif

      size_t j = 0;
      // floor soa numParticles to multiple of vecLength
      // If b is a power of 2 the following holds:
      // a & ~(b -1) == a - (a mod b)
      for (; j < (i & ~(vecLength - 1)); j += vecLength) {
        SoAKernel<true, false>(j, ownedStateI, reinterpret_cast<const autopas::OwnershipType *>(ownedStatePtr), x1, y1,
                               z1, xptr, yptr, zptr, fxptr, fyptr, fzptr, &typeIDptr[i], &typeIDptr[j], fxacc, fyacc,
                               fzacc, &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, 0);
      }
      // If b is a power of 2 the following holds:
      // a & (b -1) == a mod b
      const int rest = (int)(i & (vecLength - 1));
      if (rest > 0) {
        SoAKernel<true, true>(j, ownedStateI, reinterpret_cast<const autopas::OwnershipType *>(ownedStatePtr), x1, y1,
                              z1, xptr, yptr, zptr, fxptr, fyptr, fzptr, &typeIDptr[i], &typeIDptr[j], fxacc, fyacc,
                              fzacc, &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, rest);
      }

#if AUTOPAS_PRECISION_MODE == SPSP
      // horizontally reduce fDacc to sumfD
      const __m256 hSumfxfy = _mm256_hadd_ps(fxacc, fyacc);
      const __m256 hSumfz = _mm256_hadd_ps(fzacc, fzacc);

      const __m128 hSumfxfyLow = _mm256_extractf128_ps(hSumfxfy, 0);
      const __m128 hSumfzLow = _mm256_extractf128_ps(hSumfz, 0);

      const __m128 hSumfxfyHigh = _mm256_extractf128_ps(hSumfxfy, 1);
      const __m128 hSumfzHigh = _mm256_extractf128_ps(hSumfz, 1);

      const __m128 sumfxfyVEC = _mm_add_ps(hSumfxfyLow, hSumfxfyHigh);
      const __m128 sumfzVEC = _mm_add_ps(hSumfzLow, hSumfzHigh);

      const __m128 hsumfxfyVEC = _mm_hadd_ps(sumfxfyVEC, sumfxfyVEC);
      const __m128 hsumfzVEC = _mm_hadd_ps(sumfzVEC, sumfzVEC);

      const float sumfx = hsumfxfyVEC[0];
      const float sumfy = hsumfxfyVEC[1];
      const float sumfz = _mm_cvtss_f32(hsumfzVEC);
#else
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
#endif

      fxptr[i] += sumfx;
      fyptr[i] += sumfy;
      fzptr[i] += sumfz;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

#if AUTOPAS_PRECISION_MODE == SPSP
      // horizontally reduce virialSumX and virialSumY
      const __m256 hSumVirialxy = _mm256_hadd_ps(virialSumX, virialSumY);
      const __m128 hSumVirialxyLow = _mm256_extractf128_ps(hSumVirialxy, 0);
      const __m128 hSumVirialxyHigh = _mm256_extractf128_ps(hSumVirialxy, 1);
      const __m128 sumVirialxyVec = _mm_add_ps(hSumVirialxyHigh, hSumVirialxyLow);
      const __m128 hSumVirialxyVec = _mm_hadd_ps(sumVirialxyVec, sumVirialxyVec);

      // horizontally reduce virialSumZ and potentialEnergySum
      const __m256 hSumVirialzPotentialEnergy = _mm256_hadd_ps(virialSumZ, potentialEnergySum);
      const __m128 hSumVirialzPotentialEnergyLow = _mm256_extractf128_ps(hSumVirialzPotentialEnergy, 0);
      const __m128 hSumVirialzPotentialEnergyHigh = _mm256_extractf128_ps(hSumVirialzPotentialEnergy, 1);
      const __m128 sumVirialzPotentialEnergyVec =
          _mm_add_ps(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);
      const __m128 hSumVirialzPotentialEnergyVec =
          _mm_hadd_ps(sumVirialzPotentialEnergyVec, sumVirialzPotentialEnergyVec);

      // globals = {virialX, virialY, virialZ, potentialEnergy}
      float globals[4];
      _mm_store_ps(&globals[0], hSumVirialxyVec);
      _mm_store_ps(&globals[2], hSumVirialzPotentialEnergyVec);
#else
      // horizontally reduce virialSumX and virialSumY
      const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
      const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
      const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
      const __m128d sumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      // horizontally reduce virialSumZ and potentialEnergySum
      const __m256d hSumVirialzPotentialEnergy = _mm256_hadd_pd(virialSumZ, potentialEnergySum);
      const __m128d hSumVirialzPotentialEnergyLow = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 0);
      const __m128d hSumVirialzPotentialEnergyHigh = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 1);
      const __m128d sumVirialzPotentialEnergyVec =
          _mm_add_pd(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);

      // globals = {virialX, virialY, virialZ, potentialEnergy}
      double globals[4];
      _mm_store_pd(&globals[0], sumVirialxyVec);
      _mm_store_pd(&globals[2], sumVirialzPotentialEnergyVec);
#endif
      _aosThreadData[threadnum].virialSum[0] += globals[0];
      _aosThreadData[threadnum].virialSum[1] += globals[1];
      _aosThreadData[threadnum].virialSum[2] += globals[2];
      _aosThreadData[threadnum].potentialEnergySum += globals[3];
    }
#endif
  }

  template <bool newton3>
  inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
#ifdef __AVX__
    if (soa1.size() == 0 || soa2.size() == 0) return;

    const auto *const __restrict x1ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeID1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    const auto *const __restrict typeID2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

#if AUTOPAS_PRECISION_MODE == SPSP
    __m256 virialSumX = _mm256_setzero_ps();
    __m256 virialSumY = _mm256_setzero_ps();
    __m256 virialSumZ = _mm256_setzero_ps();
    __m256 potentialEnergySum = _mm256_setzero_ps();
#else
    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d potentialEnergySum = _mm256_setzero_pd();
#endif

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
        // If the i-th particle is a dummy, skip this loop iteration.
        continue;
      }

#if AUTOPAS_PRECISION_MODE == SPSP || AUTOPAS_PRECISION_MODE == SPDP
      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int32_t>,
                    "OwnershipStates underlying type should be int32_t!");
      // ownedStatePtr contains int32_t, so we broadcast these to make an __m256i.
      // _mm256_set1_epi32 broadcasts a 32-bit integer, we use this instruction to have 8 values!
      __m256i ownedStateI = _mm256_set1_epi32(static_cast<int32_t>(ownedStatePtr1[i]));

      const __m256 x1 = _mm256_broadcast_ss(&x1ptr[i]);
      const __m256 y1 = _mm256_broadcast_ss(&y1ptr[i]);
      const __m256 z1 = _mm256_broadcast_ss(&z1ptr[i]);
#else
      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");
      // ownedStatePtr1 contains int64_t, so we broadcast these to make an __m256i.
      // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
      __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr1[i]));

      const __m256d x1 = _mm256_broadcast_sd(&x1ptr[i]);
      const __m256d y1 = _mm256_broadcast_sd(&y1ptr[i]);
      const __m256d z1 = _mm256_broadcast_sd(&z1ptr[i]);
#endif

#if AUTOPAS_PRECISION_MODE == SPSP
      __m256 fxacc = _mm256_setzero_ps();
      __m256 fyacc = _mm256_setzero_ps();
      __m256 fzacc = _mm256_setzero_ps();
#else
      __m256d fxacc = _mm256_setzero_pd();
      __m256d fyacc = _mm256_setzero_pd();
      __m256d fzacc = _mm256_setzero_pd();
#endif

      // floor soa2 numParticles to multiple of vecLength
      unsigned int j = 0;
      for (; j < (soa2.size() & ~(vecLength - 1)); j += vecLength) {
        SoAKernel<newton3, false>(j, ownedStateI, reinterpret_cast<const autopas::OwnershipType *>(ownedStatePtr2), x1,
                                  y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, &typeID1ptr[i], &typeID2ptr[j],
                                  fxacc, fyacc, fzacc, &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, 0);
      }
      const int rest = (int)(soa2.size() & (vecLength - 1));
      if (rest > 0)
        SoAKernel<newton3, true>(j, ownedStateI, reinterpret_cast<const autopas::OwnershipType *>(ownedStatePtr2), x1,
                                 y1, z1, x2ptr, y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, &typeID1ptr[i], &typeID2ptr[j],
                                 fxacc, fyacc, fzacc, &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, rest);

#if AUTOPAS_PRECISION_MODE == SPSP
      // horizontally reduce fDacc to sumfD
      const __m256 hSumfxfy = _mm256_hadd_ps(fxacc, fyacc);
      const __m256 hSumfz = _mm256_hadd_ps(fzacc, fzacc);

      const __m128 hSumfxfyLow = _mm256_extractf128_ps(hSumfxfy, 0);
      const __m128 hSumfzLow = _mm256_extractf128_ps(hSumfz, 0);

      const __m128 hSumfxfyHigh = _mm256_extractf128_ps(hSumfxfy, 1);
      const __m128 hSumfzHigh = _mm256_extractf128_ps(hSumfz, 1);

      const __m128 sumfxfyVEC = _mm_add_ps(hSumfxfyLow, hSumfxfyHigh);
      const __m128 sumfzVEC = _mm_add_ps(hSumfzLow, hSumfzHigh);

      const __m128 hsumfxfyVEC = _mm_hadd_ps(sumfxfyVEC, sumfxfyVEC);
      const __m128 hsumfzVEC = _mm_hadd_ps(sumfzVEC, sumfzVEC);

      const float sumfx = hsumfxfyVEC[0];
      const float sumfy = hsumfxfyVEC[1];
      const float sumfz = _mm_cvtss_f32(hsumfzVEC);
#else
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
#endif

      fx1ptr[i] += sumfx;
      fy1ptr[i] += sumfy;
      fz1ptr[i] += sumfz;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

#if AUTOPAS_PRECISION_MODE == SPSP
      // horizontally reduce virialSumX and virialSumY
      const __m256 hSumVirialxy = _mm256_hadd_ps(virialSumX, virialSumY);
      const __m128 hSumVirialxyLow = _mm256_extractf128_ps(hSumVirialxy, 0);
      const __m128 hSumVirialxyHigh = _mm256_extractf128_ps(hSumVirialxy, 1);
      const __m128 sumVirialxyVec = _mm_add_ps(hSumVirialxyHigh, hSumVirialxyLow);
      const __m128 hSumVirialxyVec = _mm_hadd_ps(sumVirialxyVec, sumVirialxyVec);

      // horizontally reduce virialSumZ and potentialEnergySum
      const __m256 hSumVirialzPotentialEnergy = _mm256_hadd_ps(virialSumZ, potentialEnergySum);
      const __m128 hSumVirialzPotentialEnergyLow = _mm256_extractf128_ps(hSumVirialzPotentialEnergy, 0);
      const __m128 hSumVirialzPotentialEnergyHigh = _mm256_extractf128_ps(hSumVirialzPotentialEnergy, 1);
      const __m128 sumVirialzPotentialEnergyVec =
          _mm_add_ps(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);
      const __m128 hSumVirialzPotentialEnergyVec =
          _mm_hadd_ps(sumVirialzPotentialEnergyVec, sumVirialzPotentialEnergyVec);

      // globals = {virialX, virialY, virialZ, potentialEnergy}
      float globals[4];
      _mm_store_ps(&globals[0], hSumVirialxyVec);
      _mm_store_ps(&globals[2], hSumVirialzPotentialEnergyVec);
#else
      // horizontally reduce virialSumX and virialSumY
      const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
      const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
      const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
      const __m128d sumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      // horizontally reduce virialSumZ and potentialEnergySum
      const __m256d hSumVirialzPotentialEnergy = _mm256_hadd_pd(virialSumZ, potentialEnergySum);
      const __m128d hSumVirialzPotentialEnergyLow = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 0);
      const __m128d hSumVirialzPotentialEnergyHigh = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 1);
      const __m128d sumVirialzPotentialEnergyVec =
          _mm_add_pd(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);

      // globals = {virialX, virialY, virialZ, potentialEnergy}
      double globals[4];
      _mm_store_pd(&globals[0], sumVirialxyVec);
      _mm_store_pd(&globals[2], sumVirialzPotentialEnergyVec);
#endif

      _aosThreadData[threadnum].virialSum[0] += globals[0];
      _aosThreadData[threadnum].virialSum[1] += globals[1];
      _aosThreadData[threadnum].virialSum[2] += globals[2];
      _aosThreadData[threadnum].potentialEnergySum += globals[3];
    }
#endif
  }

#ifdef __AVX__
  /**
   * Actual inner kernel of the SoAFunctors.
   *
   * @tparam newton3
   * @tparam remainderIsMasked If false the full vector length is used. Otherwise the last entries are masked away
   * depending on the argument "rest".
   * @param j
   * @param ownedStateI
   * @param ownedStatePtr2
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
   * @param potentialEnergySum
   * @param rest
   */

  template <bool newton3, bool remainderIsMasked>
  inline void SoAKernel(const size_t j, const __m256i ownedStateI,
                        const autopas::OwnershipType *const __restrict ownedStatePtr2, const SIMDCalcType &x1,
                        const SIMDCalcType &y1, const SIMDCalcType &z1, const CalcType *const __restrict x2ptr,
                        const CalcType *const __restrict y2ptr, const CalcType *const __restrict z2ptr,
                        AccuType *const __restrict fx2ptr, AccuType *const __restrict fy2ptr,
                        AccuType *const __restrict fz2ptr, const size_t *const typeID1ptr,
                        const size_t *const typeID2ptr, SIMDAccuType &fxacc, SIMDAccuType &fyacc, SIMDAccuType &fzacc,
                        SIMDAccuType *virialSumX, SIMDAccuType *virialSumY, SIMDAccuType *virialSumZ,
                        SIMDAccuType *potentialEnergySum, const unsigned int rest = 0) {
    SIMDCalcType epsilon24s = _epsilon24;
    SIMDCalcType sigmaSquareds = _sigmaSquared;
    SIMDCalcType shift6s = _shift6;

#if AUTOPAS_PRECISION_MODE == SPSP || AUTOPAS_PRECISION_MODE == SPDP
    // code for calculations in single precision

    if constexpr (useMixing) {
      // the first argument for set lands in the last bits of the register
      epsilon24s = _mm256_set_ps(
          not remainderIsMasked or rest > 7 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 7)) : 0,
          not remainderIsMasked or rest > 6 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 6)) : 0,
          not remainderIsMasked or rest > 5 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 5)) : 0,
          not remainderIsMasked or rest > 4 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 4)) : 0,
          not remainderIsMasked or rest > 3 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 3)) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 2)) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 1)) : 0,
          _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 0)));
      sigmaSquareds = _mm256_set_ps(
          not remainderIsMasked or rest > 7 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 7)) : 0,
          not remainderIsMasked or rest > 6 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 6)) : 0,
          not remainderIsMasked or rest > 5 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 5)) : 0,
          not remainderIsMasked or rest > 4 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 4)) : 0,
          not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 3)) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 2)) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 1)) : 0,
          _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 0)));
      if constexpr (applyShift) {
        shift6s = _mm256_set_ps(
            (not remainderIsMasked or rest > 7) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 7)) : 0,
            (not remainderIsMasked or rest > 6) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 6)) : 0,
            (not remainderIsMasked or rest > 5) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 5)) : 0,
            (not remainderIsMasked or rest > 4) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 4)) : 0,
            (not remainderIsMasked or rest > 3) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 3)) : 0,
            (not remainderIsMasked or rest > 2) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 2)) : 0,
            (not remainderIsMasked or rest > 1) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 1)) : 0,
            _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 0)));
      }
    }

    const __m256 x2 = remainderIsMasked ? _mm256_maskload_ps(&x2ptr[j], _masks[rest - 1]) : _mm256_loadu_ps(&x2ptr[j]);
    const __m256 y2 = remainderIsMasked ? _mm256_maskload_ps(&y2ptr[j], _masks[rest - 1]) : _mm256_loadu_ps(&y2ptr[j]);
    const __m256 z2 = remainderIsMasked ? _mm256_maskload_ps(&z2ptr[j], _masks[rest - 1]) : _mm256_loadu_ps(&z2ptr[j]);

    const __m256 drx = _mm256_sub_ps(x1, x2);
    const __m256 dry = _mm256_sub_ps(y1, y2);
    const __m256 drz = _mm256_sub_ps(z1, z2);

    const __m256 drx2 = _mm256_mul_ps(drx, drx);
    const __m256 dry2 = _mm256_mul_ps(dry, dry);
    const __m256 drz2 = _mm256_mul_ps(drz, drz);

    const __m256 dr2PART = _mm256_add_ps(drx2, dry2);
    const __m256 dr2 = _mm256_add_ps(dr2PART, drz2);

    // _CMP_LE_OS == Less-Equal-then (ordered, signaling)
    // signaling = throw error if NaN is encountered
    // dr2 <= _cutoffSquared ? 0xFFFFFFFFFFFFFFFF : 0
    const __m256 cutoffMask = _mm256_cmp_ps(dr2, _cutoffSquared, _CMP_LE_OS);

    // This requires that dummy is zero (otherwise when loading using a mask the owned state will not be zero)
    // Uses _mm256_maskload_ps because _mm256_maskload_epi32 requires AVX2, this is a workaround
    const __m256i ownedStateJ = remainderIsMasked
                                    ? _mm256_castps_si256(_mm256_maskload_ps(
                                          reinterpret_cast<float const *>(&ownedStatePtr2[j]), _masks[rest - 1]))
                                    : _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&ownedStatePtr2[j]));

    // This requires that dummy is the first entry in OwnershipState!
    const __m256 dummyMask = _mm256_cmp_ps(_mm256_castsi256_ps(ownedStateJ), _zero, _CMP_NEQ_UQ);
    const __m256 cutoffDummyMask = _mm256_and_ps(cutoffMask, dummyMask);

    // if everything is masked away return from this function.
    if (_mm256_movemask_ps(cutoffDummyMask) == 0) {
      return;
    }

    const __m256 invdr2 = _mm256_div_ps(_one, dr2);
    const __m256 lj2 = _mm256_mul_ps(sigmaSquareds, invdr2);
    const __m256 lj4 = _mm256_mul_ps(lj2, lj2);
    const __m256 lj6 = _mm256_mul_ps(lj2, lj4);
    const __m256 lj12 = _mm256_mul_ps(lj6, lj6);
    const __m256 lj12m6 = _mm256_sub_ps(lj12, lj6);
    const __m256 lj12m6alj12 = _mm256_add_ps(lj12m6, lj12);
    const __m256 lj12m6alj12e = _mm256_mul_ps(lj12m6alj12, epsilon24s);
    const __m256 fac = _mm256_mul_ps(lj12m6alj12e, invdr2);

    const __m256 facMasked =
        remainderIsMasked ? _mm256_and_ps(fac, _mm256_and_ps(cutoffDummyMask, _mm256_castsi256_ps(_masks[rest - 1])))
                          : _mm256_and_ps(fac, cutoffDummyMask);

    const __m256 fx = _mm256_mul_ps(drx, facMasked);
    const __m256 fy = _mm256_mul_ps(dry, facMasked);
    const __m256 fz = _mm256_mul_ps(drz, facMasked);
#else
    // code for calculations in double precision

    if (useMixing) {
      // the first argument for set lands in the last bits of the register
      epsilon24s = _mm256_set_pd(
          not remainderIsMasked or rest > 3 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 3)) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 2)) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 1)) : 0,
          _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 0)));
      sigmaSquareds = _mm256_set_pd(
          not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 3)) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 2)) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 1)) : 0,
          _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 0)));
      if constexpr (applyShift) {
        shift6s = _mm256_set_pd(
            (not remainderIsMasked or rest > 3) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 3)) : 0,
            (not remainderIsMasked or rest > 2) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 2)) : 0,
            (not remainderIsMasked or rest > 1) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 1)) : 0,
            _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 0)));
      }
    }

    const __m256d x2 = remainderIsMasked ? _mm256_maskload_pd(&x2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&x2ptr[j]);
    const __m256d y2 = remainderIsMasked ? _mm256_maskload_pd(&y2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&y2ptr[j]);
    const __m256d z2 = remainderIsMasked ? _mm256_maskload_pd(&z2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&z2ptr[j]);

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
    // dr2 <= _cutoffSquared ? 0xFFFFFFFFFFFFFFFF : 0
    const __m256d cutoffMask = _mm256_cmp_pd(dr2, _cutoffSquared, _CMP_LE_OS);

    // This requires that dummy is zero (otherwise when loading using a mask the owned state will not be zero)
    const __m256i ownedStateJ = remainderIsMasked
                                    ? _mm256_castpd_si256(_mm256_maskload_pd(
                                          reinterpret_cast<double const *>(&ownedStatePtr2[j]), _masks[rest - 1]))
                                    : _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&ownedStatePtr2[j]));
    // This requires that dummy is the first entry in OwnershipState!
    const __m256d dummyMask = _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateJ), _zero, _CMP_NEQ_UQ);
    const __m256d cutoffDummyMask = _mm256_and_pd(cutoffMask, dummyMask);

    // if everything is masked away return from this function.
    if (_mm256_movemask_pd(cutoffDummyMask) == 0) {
      return;
    }

    const __m256d invdr2 = _mm256_div_pd(_one, dr2);
    const __m256d lj2 = _mm256_mul_pd(sigmaSquareds, invdr2);
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
#endif

#if AUTOPAS_PRECISION_MODE == SPSP
    // single precision

    fxacc = _mm256_add_ps(fxacc, fx);
    fyacc = _mm256_add_ps(fyacc, fy);
    fzacc = _mm256_add_ps(fzacc, fz);

    // if newton 3 is used subtract fD from particle j
    if constexpr (newton3) {
      const __m256 fx2 =
          remainderIsMasked ? _mm256_maskload_ps(&fx2ptr[j], _masks[rest - 1]) : _mm256_loadu_ps(&fx2ptr[j]);
      const __m256 fy2 =
          remainderIsMasked ? _mm256_maskload_ps(&fy2ptr[j], _masks[rest - 1]) : _mm256_loadu_ps(&fy2ptr[j]);
      const __m256 fz2 =
          remainderIsMasked ? _mm256_maskload_ps(&fz2ptr[j], _masks[rest - 1]) : _mm256_loadu_ps(&fz2ptr[j]);

      const __m256 fx2new = _mm256_sub_ps(fx2, fx);
      const __m256 fy2new = _mm256_sub_ps(fy2, fy);
      const __m256 fz2new = _mm256_sub_ps(fz2, fz);

      if (remainderIsMasked) {
        _mm256_maskstore_ps(&fx2ptr[j], _masks[rest - 1], fx2new);
        _mm256_maskstore_ps(&fy2ptr[j], _masks[rest - 1], fy2new);
        _mm256_maskstore_ps(&fz2ptr[j], _masks[rest - 1], fz2new);
      } else {
        _mm256_storeu_ps(&fx2ptr[j], fx2new);
        _mm256_storeu_ps(&fy2ptr[j], fy2new);
        _mm256_storeu_ps(&fz2ptr[j], fz2new);
      }
    }

    if constexpr (calculateGlobals) {
      // Global Virial
      const __m256 virialX = _mm256_mul_ps(fx, drx);
      const __m256 virialY = _mm256_mul_ps(fy, dry);
      const __m256 virialZ = _mm256_mul_ps(fz, drz);

      // Global Potential
      const __m256 potentialEnergy = wrapperFMA(epsilon24s, lj12m6, shift6s);

      const __m256 potentialEnergyMasked =
          remainderIsMasked
              ? _mm256_and_ps(potentialEnergy, _mm256_and_ps(cutoffDummyMask, _mm256_castsi256_ps(_masks[rest - 1])))
              : _mm256_and_ps(potentialEnergy, cutoffDummyMask);

      __m256 ownedMaskI =
          _mm256_cmp_ps(_mm256_castsi256_ps(ownedStateI), _mm256_castsi256_ps(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
      __m256 energyFactor = _mm256_blendv_ps(_zero, _one, ownedMaskI);

      if constexpr (newton3) {
        __m256 ownedMaskJ =
            _mm256_cmp_ps(_mm256_castsi256_ps(ownedStateJ), _mm256_castsi256_ps(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
        energyFactor = _mm256_add_ps(energyFactor, _mm256_blendv_ps(_zero, _one, ownedMaskJ));
      }
      *potentialEnergySum = wrapperFMA(energyFactor, potentialEnergyMasked, *potentialEnergySum);
      *virialSumX = wrapperFMA(energyFactor, virialX, *virialSumX);
      *virialSumY = wrapperFMA(energyFactor, virialY, *virialSumY);
      *virialSumZ = wrapperFMA(energyFactor, virialZ, *virialSumZ);
    }
#elif AUTOPAS_PRECISION_MODE == SPDP
    // special mixed precision case, we have to convert the 8 floats to doubles first before adding them

    // split into 2*4 floats and convert to 2*4 doubles, _mm256_cvtps_pd are expensive! _mm512_cvtps_pd introduced by
    // AVX512 could reduce this to one instruction each
    const __m256d fxCvtLow = _mm256_cvtps_pd(_mm256_extractf128_ps(fx, 0));
    const __m256d fyCvtLow = _mm256_cvtps_pd(_mm256_extractf128_ps(fy, 0));
    const __m256d fzCvtLow = _mm256_cvtps_pd(_mm256_extractf128_ps(fz, 0));

    // add the converted values
    fxacc = _mm256_add_pd(fxacc, fxCvtLow);
    fyacc = _mm256_add_pd(fyacc, fyCvtLow);
    fzacc = _mm256_add_pd(fzacc, fzCvtLow);

    __m256d fxCvtHigh;
    __m256d fyCvtHigh;
    __m256d fzCvtHigh;
    // execute the second half if needed
    if (not remainderIsMasked or (remainderIsMasked and rest > 4)) {
      fxCvtHigh = _mm256_cvtps_pd(_mm256_extractf128_ps(fx, 1));
      fyCvtHigh = _mm256_cvtps_pd(_mm256_extractf128_ps(fy, 1));
      fzCvtHigh = _mm256_cvtps_pd(_mm256_extractf128_ps(fz, 1));
      fxacc = _mm256_add_pd(fxacc, fxCvtHigh);
      fyacc = _mm256_add_pd(fyacc, fyCvtHigh);
      fzacc = _mm256_add_pd(fzacc, fzCvtHigh);
    }

    // if newton 3 is used subtract fD from particle j
    if constexpr (newton3) {
      // we need to potentially load and store double the number of particles in this case
      const __m256d fx2first = _mm256_loadu_pd(&fx2ptr[j]);
      const __m256d fy2first = _mm256_loadu_pd(&fy2ptr[j]);
      const __m256d fz2first = _mm256_loadu_pd(&fz2ptr[j]);

      const __m256d fx2newfirst = _mm256_sub_pd(fx2first, fxCvtLow);
      const __m256d fy2newfirst = _mm256_sub_pd(fy2first, fyCvtLow);
      const __m256d fz2newfirst = _mm256_sub_pd(fz2first, fzCvtLow);

      _mm256_storeu_pd(&fx2ptr[j], fx2newfirst);
      _mm256_storeu_pd(&fy2ptr[j], fy2newfirst);
      _mm256_storeu_pd(&fz2ptr[j], fz2newfirst);

      // execute the second half if needed
      if (not remainderIsMasked or (remainderIsMasked and rest > 4)) {
        const __m256d fx2second = _mm256_loadu_pd(&fx2ptr[j + 4]);
        const __m256d fy2second = _mm256_loadu_pd(&fy2ptr[j + 4]);
        const __m256d fz2second = _mm256_loadu_pd(&fz2ptr[j + 4]);
        const __m256d fx2newsecond = _mm256_sub_pd(fx2second, fxCvtHigh);
        const __m256d fy2newsecond = _mm256_sub_pd(fy2second, fyCvtHigh);
        const __m256d fz2newsecond = _mm256_sub_pd(fz2second, fzCvtHigh);
        _mm256_storeu_pd(&fx2ptr[j + 4], fx2newsecond);
        _mm256_storeu_pd(&fy2ptr[j + 4], fy2newsecond);
        _mm256_storeu_pd(&fz2ptr[j + 4], fz2newsecond);
      }
    }

    if constexpr (calculateGlobals) {
      // Global Virial
      const __m256 virialX = _mm256_mul_ps(fx, drx);
      const __m256 virialY = _mm256_mul_ps(fy, dry);
      const __m256 virialZ = _mm256_mul_ps(fz, drz);

      // Global Potential
      const __m256 potentialEnergy = wrapperFMA(epsilon24s, lj12m6, shift6s);

      const __m256 potentialEnergyMasked =
          remainderIsMasked
              ? _mm256_and_ps(potentialEnergy, _mm256_and_ps(cutoffDummyMask, _mm256_castsi256_ps(_masks[rest - 1])))
              : _mm256_and_ps(potentialEnergy, cutoffDummyMask);

      const __m256 ownedMaskI =
          _mm256_cmp_ps(_mm256_castsi256_ps(ownedStateI), _mm256_castsi256_ps(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
      __m256 energyFactor = _mm256_blendv_ps(_zero, _one, ownedMaskI);
      if constexpr (newton3) {
        __m256 ownedMaskJ =
            _mm256_cmp_ps(_mm256_castsi256_ps(ownedStateJ), _mm256_castsi256_ps(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
        energyFactor = _mm256_add_ps(energyFactor, _mm256_blendv_ps(_zero, _one, ownedMaskJ));
      }
      const __m256d energyFactorLow = _mm256_cvtps_pd(_mm256_extractf128_ps(energyFactor, 0));
      const __m256d energyFactorHigh = _mm256_cvtps_pd(_mm256_extractf128_ps(energyFactor, 1));
      *potentialEnergySum = wrapperFMA(
          energyFactorLow, _mm256_cvtps_pd(_mm256_extractf128_ps(potentialEnergyMasked, 0)), *potentialEnergySum);
      *potentialEnergySum = wrapperFMA(
          energyFactorHigh, _mm256_cvtps_pd(_mm256_extractf128_ps(potentialEnergyMasked, 1)), *potentialEnergySum);
      *virialSumX = wrapperFMA(energyFactorLow, _mm256_cvtps_pd(_mm256_extractf128_ps(virialX, 0)), *virialSumX);
      *virialSumX = wrapperFMA(energyFactorHigh, _mm256_cvtps_pd(_mm256_extractf128_ps(virialX, 1)), *virialSumX);
      *virialSumY = wrapperFMA(energyFactorLow, _mm256_cvtps_pd(_mm256_extractf128_ps(virialY, 0)), *virialSumY);
      *virialSumY = wrapperFMA(energyFactorHigh, _mm256_cvtps_pd(_mm256_extractf128_ps(virialY, 1)), *virialSumY);
      *virialSumZ = wrapperFMA(energyFactorLow, _mm256_cvtps_pd(_mm256_extractf128_ps(virialZ, 0)), *virialSumZ);
      *virialSumZ = wrapperFMA(energyFactorHigh, _mm256_cvtps_pd(_mm256_extractf128_ps(virialZ, 1)), *virialSumZ);
    }
#else
    // double precision

    fxacc = _mm256_add_pd(fxacc, fx);
    fyacc = _mm256_add_pd(fyacc, fy);
    fzacc = _mm256_add_pd(fzacc, fz);

    // if newton 3 is used subtract fD from particle j
    if constexpr (newton3) {
      const __m256d fx2 =
          remainderIsMasked ? _mm256_maskload_pd(&fx2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&fx2ptr[j]);
      const __m256d fy2 =
          remainderIsMasked ? _mm256_maskload_pd(&fy2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&fy2ptr[j]);
      const __m256d fz2 =
          remainderIsMasked ? _mm256_maskload_pd(&fz2ptr[j], _masks[rest - 1]) : _mm256_loadu_pd(&fz2ptr[j]);

      const __m256d fx2new = _mm256_sub_pd(fx2, fx);
      const __m256d fy2new = _mm256_sub_pd(fy2, fy);
      const __m256d fz2new = _mm256_sub_pd(fz2, fz);

      remainderIsMasked ? _mm256_maskstore_pd(&fx2ptr[j], _masks[rest - 1], fx2new)
                        : _mm256_storeu_pd(&fx2ptr[j], fx2new);
      remainderIsMasked ? _mm256_maskstore_pd(&fy2ptr[j], _masks[rest - 1], fy2new)
                        : _mm256_storeu_pd(&fy2ptr[j], fy2new);
      remainderIsMasked ? _mm256_maskstore_pd(&fz2ptr[j], _masks[rest - 1], fz2new)
                        : _mm256_storeu_pd(&fz2ptr[j], fz2new);
    }

    if constexpr (calculateGlobals) {
      // Global Virial
      const __m256d virialX = _mm256_mul_pd(fx, drx);
      const __m256d virialY = _mm256_mul_pd(fy, dry);
      const __m256d virialZ = _mm256_mul_pd(fz, drz);

      // Global Potential
      const __m256d potentialEnergy = wrapperFMA(epsilon24s, lj12m6, shift6s);

      const __m256d potentialEnergyMasked =
          remainderIsMasked
              ? _mm256_and_pd(potentialEnergy, _mm256_and_pd(cutoffDummyMask, _mm256_castsi256_pd(_masks[rest - 1])))
              : _mm256_and_pd(potentialEnergy, cutoffDummyMask);

      __m256d ownedMaskI =
          _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateI), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
      __m256d energyFactor = _mm256_blendv_pd(_zero, _one, ownedMaskI);
      if constexpr (newton3) {
        __m256d ownedMaskJ =
            _mm256_cmp_pd(_mm256_castsi256_pd(ownedStateJ), _mm256_castsi256_pd(_ownedStateOwnedMM256i), _CMP_EQ_UQ);
        energyFactor = _mm256_add_pd(energyFactor, _mm256_blendv_pd(_zero, _one, ownedMaskJ));
      }
      *potentialEnergySum = wrapperFMA(energyFactor, potentialEnergyMasked, *potentialEnergySum);
      *virialSumX = wrapperFMA(energyFactor, virialX, *virialSumX);
      *virialSumY = wrapperFMA(energyFactor, virialY, *virialSumY);
      *virialSumZ = wrapperFMA(energyFactor, virialZ, *virialSumZ);
    }
#endif
  }
#endif

 public:
  // clang-format off
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
   */
  // clang-format on
  inline void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                               const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                               bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

 private:
  template <bool newton3>
  inline void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                                   const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
#ifdef __AVX__
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();
    if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
      return;
    }

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

#if AUTOPAS_PRECISION_MODE == SPSP
    // accumulators
    __m256 virialSumX = _mm256_setzero_ps();
    __m256 virialSumY = _mm256_setzero_ps();
    __m256 virialSumZ = _mm256_setzero_ps();
    __m256 potentialEnergySum = _mm256_setzero_ps();
    __m256 fxacc = _mm256_setzero_ps();
    __m256 fyacc = _mm256_setzero_ps();
    __m256 fzacc = _mm256_setzero_ps();
#else
    // accumulators
    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d potentialEnergySum = _mm256_setzero_pd();
    __m256d fxacc = _mm256_setzero_pd();
    __m256d fyacc = _mm256_setzero_pd();
    __m256d fzacc = _mm256_setzero_pd();
#endif

#if AUTOPAS_PRECISION_MODE == SPSP || AUTOPAS_PRECISION_MODE == SPDP
    // broadcast particle 1
    const __m256 x1 = _mm256_broadcast_ss(&xptr[indexFirst]);
    const __m256 y1 = _mm256_broadcast_ss(&yptr[indexFirst]);
    const __m256 z1 = _mm256_broadcast_ss(&zptr[indexFirst]);

    // ownedStatePtr contains int32_t, so we broadcast these to make an __m256i.
    // _mm256_set1_epi32 broadcasts a 32-bit integer, we use this instruction to have 8 values!
    __m256i ownedStateI = _mm256_set1_epi32(static_cast<autopas::OwnershipType>(ownedStatePtr[indexFirst]));

    alignas(64) std::array<CalcType, vecLength> x2tmp{};
    alignas(64) std::array<CalcType, vecLength> y2tmp{};
    alignas(64) std::array<CalcType, vecLength> z2tmp{};
    alignas(64) std::array<AccuType, vecLength> fx2tmp{};
    alignas(64) std::array<AccuType, vecLength> fy2tmp{};
    alignas(64) std::array<AccuType, vecLength> fz2tmp{};
    alignas(64) std::array<size_t, vecLength> typeID2tmp{};
    alignas(64) std::array<autopas::OwnershipState, vecLength> ownedStates2tmp{};
#else
    // broadcast particle 1
    const __m256d x1 = _mm256_broadcast_sd(&xptr[indexFirst]);
    const __m256d y1 = _mm256_broadcast_sd(&yptr[indexFirst]);
    const __m256d z1 = _mm256_broadcast_sd(&zptr[indexFirst]);

    // ownedStatePtr contains int64_t, so we broadcast these to make an __m256i.
    // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
    __m256i ownedStateI = _mm256_set1_epi64x(static_cast<autopas::OwnershipType>(ownedStatePtr[indexFirst]));

    alignas(64) std::array<double, vecLength> x2tmp{};
    alignas(64) std::array<double, vecLength> y2tmp{};
    alignas(64) std::array<double, vecLength> z2tmp{};
    alignas(64) std::array<double, vecLength> fx2tmp{};
    alignas(64) std::array<double, vecLength> fy2tmp{};
    alignas(64) std::array<double, vecLength> fz2tmp{};
    alignas(64) std::array<size_t, vecLength> typeID2tmp{};
    alignas(64) std::array<autopas::OwnershipState, vecLength> ownedStates2tmp{};
#endif

    // load 4 or 8 neighbors
    size_t j = 0;
    // Loop over all neighbors as long as we can fill full vectors
    // (until `neighborList.size() - neighborList.size() % vecLength`)
    //
    // If b is a power of 2 the following holds:
    // a & ~(b - 1) == a - (a mod b)
    for (; j < (neighborList.size() & ~(vecLength - 1)); j += vecLength) {
      // AVX2 variant:
      // create buffer for 4 interaction particles
      // and fill buffers via gathering
      //      const __m256d x2tmp = _mm256_i64gather_pd(&xptr[j], _vindex, 1);
      //      const __m256d y2tmp = _mm256_i64gather_pd(&yptr[j], _vindex, 1);
      //      const __m256d z2tmp = _mm256_i64gather_pd(&zptr[j], _vindex, 1);
      //      const __m256d fx2tmp = _mm256_i64gather_pd(&fxptr[j], _vindex, 1);
      //      const __m256d fy2tmp = _mm256_i64gather_pd(&fyptr[j], _vindex, 1);
      //      const __m256d fz2tmp = _mm256_i64gather_pd(&fzptr[j], _vindex, 1);
      //      const __m256i typeID2tmp = _mm256_i64gather_epi64(&typeIDptr[j], _vindex, 1);

      for (size_t vecIndex = 0; vecIndex < vecLength; ++vecIndex) {
        x2tmp[vecIndex] = xptr[neighborList[j + vecIndex]];
        y2tmp[vecIndex] = yptr[neighborList[j + vecIndex]];
        z2tmp[vecIndex] = zptr[neighborList[j + vecIndex]];
        if constexpr (newton3) {
          fx2tmp[vecIndex] = fxptr[neighborList[j + vecIndex]];
          fy2tmp[vecIndex] = fyptr[neighborList[j + vecIndex]];
          fz2tmp[vecIndex] = fzptr[neighborList[j + vecIndex]];
        }
        typeID2tmp[vecIndex] = typeIDptr[neighborList[j + vecIndex]];
        ownedStates2tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
      }

      SoAKernel<newton3, false>(
          0, ownedStateI, reinterpret_cast<const autopas::OwnershipType *>(ownedStates2tmp.data()), x1, y1, z1,
          x2tmp.data(), y2tmp.data(), z2tmp.data(), fx2tmp.data(), fy2tmp.data(), fz2tmp.data(), &typeIDptr[indexFirst],
          typeID2tmp.data(), fxacc, fyacc, fzacc, &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, 0);

      if constexpr (newton3) {
        for (size_t vecIndex = 0; vecIndex < vecLength; ++vecIndex) {
          fxptr[neighborList[j + vecIndex]] = fx2tmp[vecIndex];
          fyptr[neighborList[j + vecIndex]] = fy2tmp[vecIndex];
          fzptr[neighborList[j + vecIndex]] = fz2tmp[vecIndex];
        }
      }
    }
    // Remainder loop
    // If b is a power of 2 the following holds:
    // a & (b - 1) == a mod b
    const auto rest = static_cast<int>(neighborList.size() & (vecLength - 1));
    if (rest > 0) {
      // AVX2 variant:
      // create buffer for 4 interaction particles
      // and fill buffers via gathering
      //      TODO: use masked load because there will not be enough data left for the whole gather
      //      const __m256d x2tmp = _mm256_i64gather_pd(&xptr[j], _vindex, 1);
      //      const __m256d y2tmp = _mm256_i64gather_pd(&yptr[j], _vindex, 1);
      //      const __m256d z2tmp = _mm256_i64gather_pd(&zptr[j], _vindex, 1);
      //      const __m256d fx2tmp = _mm256_i64gather_pd(&fxptr[j], _vindex, 1);
      //      const __m256d fy2tmp = _mm256_i64gather_pd(&fyptr[j], _vindex, 1);
      //      const __m256d fz2tmp = _mm256_i64gather_pd(&fzptr[j], _vindex, 1);
      //      const __m256d typeID2tmp = _mm256_i64gather_pd(&typeIDptr[j], _vindex, 1);

      for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
        x2tmp[vecIndex] = xptr[neighborList[j + vecIndex]];
        y2tmp[vecIndex] = yptr[neighborList[j + vecIndex]];
        z2tmp[vecIndex] = zptr[neighborList[j + vecIndex]];
        // if newton3 is used we need to load f of particle j so the kernel can update it too
        if constexpr (newton3) {
          fx2tmp[vecIndex] = fxptr[neighborList[j + vecIndex]];
          fy2tmp[vecIndex] = fyptr[neighborList[j + vecIndex]];
          fz2tmp[vecIndex] = fzptr[neighborList[j + vecIndex]];
        }
        typeID2tmp[vecIndex] = typeIDptr[neighborList[j + vecIndex]];
        ownedStates2tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
      }

      SoAKernel<newton3, true>(0, ownedStateI, reinterpret_cast<const autopas::OwnershipType *>(ownedStates2tmp.data()),
                               x1, y1, z1, x2tmp.data(), y2tmp.data(), z2tmp.data(), fx2tmp.data(), fy2tmp.data(),
                               fz2tmp.data(), &typeIDptr[indexFirst], typeID2tmp.data(), fxacc, fyacc, fzacc,
                               &virialSumX, &virialSumY, &virialSumZ, &potentialEnergySum, rest);

      if constexpr (newton3) {
        for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
          fxptr[neighborList[j + vecIndex]] = fx2tmp[vecIndex];
          fyptr[neighborList[j + vecIndex]] = fy2tmp[vecIndex];
          fzptr[neighborList[j + vecIndex]] = fz2tmp[vecIndex];
        }
      }
    }

#if AUTOPAS_PRECISION_MODE == SPSP
    // horizontally reduce fDacc to sumfD
    const __m256 hSumfxfy = _mm256_hadd_ps(fxacc, fyacc);
    const __m256 hSumfz = _mm256_hadd_ps(fzacc, fzacc);

    const __m128 hSumfxfyLow = _mm256_extractf128_ps(hSumfxfy, 0);
    const __m128 hSumfzLow = _mm256_extractf128_ps(hSumfz, 0);

    const __m128 hSumfxfyHigh = _mm256_extractf128_ps(hSumfxfy, 1);
    const __m128 hSumfzHigh = _mm256_extractf128_ps(hSumfz, 1);

    const __m128 sumfxfyVEC = _mm_add_ps(hSumfxfyLow, hSumfxfyHigh);
    const __m128 sumfzVEC = _mm_add_ps(hSumfzLow, hSumfzHigh);

    const __m128 hsumfxfyVEC = _mm_hadd_ps(sumfxfyVEC, sumfxfyVEC);
    const __m128 hsumfzVEC = _mm_hadd_ps(sumfzVEC, sumfzVEC);

    const float sumfx = hsumfxfyVEC[0];
    const float sumfy = hsumfxfyVEC[1];
    const float sumfz = _mm_cvtss_f32(hsumfzVEC);
#else
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
#endif

    fxptr[indexFirst] += sumfx;
    fyptr[indexFirst] += sumfy;
    fzptr[indexFirst] += sumfz;

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

#if AUTOPAS_PRECISION_MODE == SPSP
      // horizontally reduce virialSumX and virialSumY
      const __m256 hSumVirialxy = _mm256_hadd_ps(virialSumX, virialSumY);
      const __m128 hSumVirialxyLow = _mm256_extractf128_ps(hSumVirialxy, 0);
      const __m128 hSumVirialxyHigh = _mm256_extractf128_ps(hSumVirialxy, 1);
      const __m128 sumVirialxyVec = _mm_add_ps(hSumVirialxyHigh, hSumVirialxyLow);
      const __m128 hSumVirialxyVec = _mm_hadd_ps(sumVirialxyVec, sumVirialxyVec);

      // horizontally reduce virialSumZ and potentialEnergySum
      const __m256 hSumVirialzPotentialEnergy = _mm256_hadd_ps(virialSumZ, potentialEnergySum);
      const __m128 hSumVirialzPotentialEnergyLow = _mm256_extractf128_ps(hSumVirialzPotentialEnergy, 0);
      const __m128 hSumVirialzPotentialEnergyHigh = _mm256_extractf128_ps(hSumVirialzPotentialEnergy, 1);
      const __m128 sumVirialzPotentialEnergyVec =
          _mm_add_ps(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);
      const __m128 hSumVirialzPotentialEnergyVec =
          _mm_hadd_ps(sumVirialzPotentialEnergyVec, sumVirialzPotentialEnergyVec);

      // globals = {virialX, virialY, virialZ, potentialEnergy}
      float globals[4];
      _mm_store_ps(&globals[0], hSumVirialxyVec);
      _mm_store_ps(&globals[2], hSumVirialzPotentialEnergyVec);
#else
      // horizontally reduce virialSumX and virialSumY
      const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
      const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
      const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
      const __m128d sumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      // horizontally reduce virialSumZ and potentialEnergySum
      const __m256d hSumVirialzPotentialEnergy = _mm256_hadd_pd(virialSumZ, potentialEnergySum);
      const __m128d hSumVirialzPotentialEnergyLow = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 0);
      const __m128d hSumVirialzPotentialEnergyHigh = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 1);
      const __m128d sumVirialzPotentialEnergyVec =
          _mm_add_pd(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);

      // globals = {virialX, virialY, virialZ, potentialEnergy}
      double globals[4];
      _mm_store_pd(&globals[0], sumVirialxyVec);
      _mm_store_pd(&globals[2], sumVirialzPotentialEnergyVec);
#endif
      _aosThreadData[threadnum].virialSum[0] += globals[0];
      _aosThreadData[threadnum].virialSum[1] += globals[1];
      _aosThreadData[threadnum].virialSum[2] += globals[2];
      _aosThreadData[threadnum].potentialEnergySum += globals[3];
    }
    // interact with i with 4 or 8 neighbors
#endif  // __AVX__
  }

 public:
  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 9>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 6>{
        Particle::AttributeNames::id,   Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
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
   * Reset the global values.
   * Will set the global values to zero to prepare for the next iteration.
   */
  void initTraversal() final {
    _potentialEnergySum = 0.;
    _virialSum = {0., 0., 0.};
    _postProcessed = false;
    for (size_t i = 0; i < _aosThreadData.size(); ++i) {
      _aosThreadData[i].setZero();
    }
  }

  /**
   * Accumulates global values, e.g. potential energy and virial.
   * @param newton3
   */
  void endTraversal(bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (_postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }
    if (calculateGlobals) {
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        _potentialEnergySum += _aosThreadData[i].potentialEnergySum;
        _virialSum += _aosThreadData[i].virialSum;
      }
      // For each interaction, we added the full contribution for both particles. Divide by 2 here, so that each
      // contribution is only counted once per pair.
      _potentialEnergySum *= static_cast<AccuType>(0.5);
      _virialSum *= static_cast<AccuType>(0.5);

      // We have always calculated 6*potentialEnergy, so we divide by 6 here!
      _potentialEnergySum /= static_cast<AccuType>(6.);
      _postProcessed = true;

      AutoPasLog(INFO, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(INFO, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
    }
  }

  /**
   * Get the potential Energy
   * @return the potential Energy
   */
  AccuType getPotentialEnergy() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get potential energy even though calculateGlobals is false. If you want this functor to "
          "calculate "
          "global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Cannot get potential energy, because endTraversal was not called.");
    }
    return _potentialEnergySum;
  }

  /**
   * Get the virial
   * @return the virial
   */
  AccuType getVirial() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get virial even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Cannot get virial, because endTraversal was not called.");
    }
    return _virialSum[0] + _virialSum[1] + _virialSum[2];
  }

  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param epsilon24
   * @param sigmaSquared
   */
  void setParticleProperties(CalcType epsilon24, CalcType sigmaSquared) {
#ifdef __AVX__
#if AUTOPAS_PRECISION_MODE == SPSP || AUTOPAS_PRECISION_MODE == SPDP
    _epsilon24 = _mm256_set1_ps(epsilon24);
    _sigmaSquared = _mm256_set1_ps(sigmaSquared);
    if constexpr (applyShift) {
      _shift6 = _mm256_set1_ps(
          ParticlePropertiesLibrary<CalcType, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquared[0]));
    } else {
      _shift6 = _mm256_setzero_ps();
    }
#else
    _epsilon24 = _mm256_set1_pd(epsilon24);
    _sigmaSquared = _mm256_set1_pd(sigmaSquared);
    if constexpr (applyShift) {
      _shift6 = _mm256_set1_pd(
          ParticlePropertiesLibrary<CalcType, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquared[0]));
    } else {
      _shift6 = _mm256_setzero_pd();
    }
#endif
#endif

    _epsilon24AoS = epsilon24;
    _sigmaSquaredAoS = sigmaSquared;
    if constexpr (applyShift) {
      _shift6AoS = ParticlePropertiesLibrary<CalcType, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquaredAoS);
    } else {
      _shift6AoS = 0.;
    }
  }

 private:
#ifdef __AVX__
  /**
   * Wrapper function for FMA. If FMA is not supported it first executes the multiplication then the addition.
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
#endif

  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }

    // variables
    std::array<AccuType, 3> virialSum;
    AccuType potentialEnergySum;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 4 * sizeof(AccuType)) / sizeof(double)];
  };
  // make sure of the size of AoSThreadData
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

#ifdef __AVX__
#if AUTOPAS_PRECISION_MODE == SPSP || AUTOPAS_PRECISION_MODE == SPDP
  const __m256 _zero{_mm256_set1_ps(0.)};
  const __m256 _one{_mm256_set1_ps(1.)};
  const __m256i _vindex = _mm256_set_epi64x(0, 1, 3, 4);
  const __m256i _masks[7]{
      _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, -1),      _mm256_set_epi32(0, 0, 0, 0, 0, 0, -1, -1),
      _mm256_set_epi32(0, 0, 0, 0, 0, -1, -1, -1),    _mm256_set_epi32(0, 0, 0, 0, -1, -1, -1, -1),
      _mm256_set_epi32(0, 0, 0, -1, -1, -1, -1, -1),  _mm256_set_epi32(0, 0, -1, -1, -1, -1, -1, -1),
      _mm256_set_epi32(0, -1, -1, -1, -1, -1, -1, -1)};
  const __m256i _ownedStateDummyMM256i{0x0};
  const __m256i _ownedStateOwnedMM256i{_mm256_set1_epi32(static_cast<int32_t>(autopas::OwnershipState::owned))};
  const __m256 _cutoffSquared{};
  __m256 _shift6 = _mm256_setzero_ps();
  __m256 _epsilon24{};
  __m256 _sigmaSquared{};
#else
  const __m256d _zero{_mm256_set1_pd(0.)};
  const __m256d _one{_mm256_set1_pd(1.)};
  // currently unused
  const __m256i _vindex = _mm256_set_epi64x(0, 1, 3, 4);
  const __m256i _masks[3]{
      _mm256_set_epi64x(0, 0, 0, -1),
      _mm256_set_epi64x(0, 0, -1, -1),
      _mm256_set_epi64x(0, -1, -1, -1),
  };
  const __m256i _ownedStateDummyMM256i{0x0};
  const __m256i _ownedStateOwnedMM256i{_mm256_set1_epi64x(static_cast<int64_t>(autopas::OwnershipState::owned))};
  const __m256d _cutoffSquared{};
  __m256d _shift6 = _mm256_setzero_pd();
  __m256d _epsilon24{};
  __m256d _sigmaSquared{};
#endif
#endif

  const CalcType _cutoffSquaredAoS = 0;
  CalcType _epsilon24AoS, _sigmaSquaredAoS, _shift6AoS = 0;

  ParticlePropertiesLibrary<CalcType, size_t> *_PPLibrary = nullptr;

  // sum of the potential energy, only calculated if calculateGlobals is true
  AccuType _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<AccuType, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

// number of double values that fit into a vector register.
// MUST be power of 2 because some optimizations make this assumption
#if AUTOPAS_PRECISION_MODE == SPSP || AUTOPAS_PRECISION_MODE == SPDP
  constexpr static size_t vecLength = 8;
#else
  constexpr static size_t vecLength = 4;
#endif
};
}  // namespace mdLib
