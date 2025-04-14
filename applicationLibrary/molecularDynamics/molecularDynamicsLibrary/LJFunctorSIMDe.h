
#pragma once

#include <x86/avx.h>

#include <array>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
// for fused-multiply add needs separate import
#include <x86/fma.h>

namespace mdLib {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * This Version is implemented using SIMDe intrinsics.
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
class LJFunctorSIMDe
    : public autopas::PairwiseFunctor<Particle, LJFunctorSIMDe<Particle, applyShift, useMixing, useNewton3,
                                                               calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorSIMDe() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy unused, only there to make the signature different from the public constructor.
   */
  explicit LJFunctorSIMDe(double cutoff, void * /*dummy*/)
      : autopas::PairwiseFunctor<Particle, LJFunctorSIMDe<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
                                                          countFLOPs, relevantForTuning>>(cutoff),
        _cutoffsquare{simde_mm256_set1_pd(cutoff * cutoff)},
        _cutoffsquareAoS(cutoff * cutoff),
        _upotSum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
  }

 public:
  /**
   * Constructor for Functor with mixing disabled. When using this functor it is necessary to call
   * setParticleProperties() to set internal constants because it does not use a particle properties library.
   *
   * @note Only to be used with mixing == false.
   *
   * @param cutoff
   */
  explicit LJFunctorSIMDe(double cutoff) : LJFunctorSIMDe(cutoff, nullptr) {
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
  explicit LJFunctorSIMDe(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctorSIMDe(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  std::string getName() final { return "LJFunctorSIMDe"; }

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
    auto sigmasquare = _sigmaSquareAoS;
    auto epsilon24 = _epsilon24AoS;
    auto shift6 = _shift6AoS;
    if constexpr (useMixing) {
      sigmasquare = _PPLibrary->getMixingSigmaSquared(i.getTypeId(), j.getTypeId());
      epsilon24 = _PPLibrary->getMixing24Epsilon(i.getTypeId(), j.getTypeId());
      if constexpr (applyShift) {
        shift6 = _PPLibrary->getMixingShift6(i.getTypeId(), j.getTypeId());
      }
    }
    auto dr = i.getR() - j.getR();
    double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffsquareAoS) {
      return;
    }

    double invdr2 = 1. / dr2;
    double lj6 = sigmasquare * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    auto f = dr * fac;
    i.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(f);
    }
    if (calculateGlobals) {
      auto virial = dr * f;
      double upot = epsilon24 * lj12m6 + shift6;

      const int threadnum = autopas::autopas_get_thread_num();
      // for non-newton3 the division is in the post-processing step.
      if (newton3) {
        upot *= 0.5;
        virial *= (double)0.5;
      }
      if (i.isOwned()) {
        _aosThreadData[threadnum].upotSum += upot;
        _aosThreadData[threadnum].virialSum += virial;
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        _aosThreadData[threadnum].upotSum += upot;
        _aosThreadData[threadnum].virialSum += virial;
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
    if (soa.size() == 0) return;

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

    simde__m256d virialSumX = simde_mm256_setzero_pd();
    simde__m256d virialSumY = simde_mm256_setzero_pd();
    simde__m256d virialSumZ = simde_mm256_setzero_pd();
    simde__m256d upotSum = simde_mm256_setzero_pd();

    // reverse outer loop s.th. inner loop always beginns at aligned array start
    // typecast to detect underflow
    for (size_t i = soa.size() - 1; (long)i >= 0; --i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        // If the i-th particle is a dummy, skip this loop iteration.
        continue;
      }

      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");
      // ownedStatePtr contains int64_t, so we broadcast these to make an simde__m256i.
      // simde_mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
      simde__m256i ownedStateI = simde_mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr[i]));

      simde__m256d fxacc = simde_mm256_setzero_pd();
      simde__m256d fyacc = simde_mm256_setzero_pd();
      simde__m256d fzacc = simde_mm256_setzero_pd();

      const simde__m256d x1 = simde_mm256_broadcast_sd(&xptr[i]);
      const simde__m256d y1 = simde_mm256_broadcast_sd(&yptr[i]);
      const simde__m256d z1 = simde_mm256_broadcast_sd(&zptr[i]);

      size_t j = 0;
      // floor soa numParticles to multiple of vecLength
      // If b is a power of 2 the following holds:
      // a & ~(b -1) == a - (a mod b)
      for (; j < (i & ~(vecLength - 1)); j += 4) {
        SoAKernel<true, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr,
                               zptr, fxptr, fyptr, fzptr, &typeIDptr[i], typeIDptr, fxacc, fyacc, fzacc, &virialSumX,
                               &virialSumY, &virialSumZ, &upotSum, 0);
      }
      // If b is a power of 2 the following holds:
      // a & (b -1) == a mod b
      const int rest = (int)(i & (vecLength - 1));
      if (rest > 0) {
        SoAKernel<true, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr,
                              zptr, fxptr, fyptr, fzptr, &typeIDptr[i], typeIDptr, fxacc, fyacc, fzacc, &virialSumX,
                              &virialSumY, &virialSumZ, &upotSum, rest);
      }

      // horizontally reduce fDacc to sumfD
      const simde__m256d hSumfxfy = simde_mm256_hadd_pd(fxacc, fyacc);
      const simde__m256d hSumfz = simde_mm256_hadd_pd(fzacc, fzacc);

      const simde__m128d hSumfxfyLow = simde_mm256_extractf128_pd(hSumfxfy, 0);
      const simde__m128d hSumfzLow = simde_mm256_extractf128_pd(hSumfz, 0);

      const simde__m128d hSumfxfyHigh = simde_mm256_extractf128_pd(hSumfxfy, 1);
      const simde__m128d hSumfzHigh = simde_mm256_extractf128_pd(hSumfz, 1);

      const simde__m128d sumfxfyVEC = simde_mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
      const simde__m128d sumfzVEC = simde_mm_add_pd(hSumfzLow, hSumfzHigh);

      const double sumfx = sumfxfyVEC[0];
      const double sumfy = sumfxfyVEC[1];
      const double sumfz = simde_mm_cvtsd_f64(sumfzVEC);

      fxptr[i] += sumfx;
      fyptr[i] += sumfy;
      fzptr[i] += sumfz;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      // horizontally reduce virialSumX and virialSumY
      const simde__m256d hSumVirialxy = simde_mm256_hadd_pd(virialSumX, virialSumY);
      const simde__m128d hSumVirialxyLow = simde_mm256_extractf128_pd(hSumVirialxy, 0);
      const simde__m128d hSumVirialxyHigh = simde_mm256_extractf128_pd(hSumVirialxy, 1);
      const simde__m128d hSumVirialxyVec = simde_mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      // horizontally reduce virialSumZ and upotSum
      const simde__m256d hSumVirialzUpot = simde_mm256_hadd_pd(virialSumZ, upotSum);
      const simde__m128d hSumVirialzUpotLow = simde_mm256_extractf128_pd(hSumVirialzUpot, 0);
      const simde__m128d hSumVirialzUpotHigh = simde_mm256_extractf128_pd(hSumVirialzUpot, 1);
      const simde__m128d hSumVirialzUpotVec = simde_mm_add_pd(hSumVirialzUpotHigh, hSumVirialzUpotLow);

      // globals = {virialX, virialY, virialZ, uPot}
      double globals[4];
      simde_mm_store_pd(&globals[0], hSumVirialxyVec);
      simde_mm_store_pd(&globals[2], hSumVirialzUpotVec);

      double factor = 1.;
      // we assume newton3 to be enabled in this function call, thus we multiply by two if the value of newton3 is
      // false, since for newton3 disabled we divide by two later on.
      factor *= newton3 ? .5 : 1.;
      // In case we have a non-cell-wise owned state, we have multiplied everything by two, so we divide it by 2 again.
      _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
      _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
      _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
      _aosThreadData[threadnum].upotSum += globals[3] * factor;
    }
  }

  template <bool newton3>
  inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
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

    simde__m256d virialSumX = simde_mm256_setzero_pd();
    simde__m256d virialSumY = simde_mm256_setzero_pd();
    simde__m256d virialSumZ = simde_mm256_setzero_pd();
    simde__m256d upotSum = simde_mm256_setzero_pd();

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
        // If the i-th particle is a dummy, skip this loop iteration.
        continue;
      }

      simde__m256d fxacc = simde_mm256_setzero_pd();
      simde__m256d fyacc = simde_mm256_setzero_pd();
      simde__m256d fzacc = simde_mm256_setzero_pd();

      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");
      // ownedStatePtr1 contains int64_t, so we broadcast these to make an simde__m256i.
      // simde_mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
      simde__m256i ownedStateI = simde_mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr1[i]));

      const simde__m256d x1 = simde_mm256_broadcast_sd(&x1ptr[i]);
      const simde__m256d y1 = simde_mm256_broadcast_sd(&y1ptr[i]);
      const simde__m256d z1 = simde_mm256_broadcast_sd(&z1ptr[i]);

      // floor soa2 numParticles to multiple of vecLength
      unsigned int j = 0;
      for (; j < (soa2.size() & ~(vecLength - 1)); j += 4) {
        SoAKernel<newton3, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                  y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc,
                                  &virialSumX, &virialSumY, &virialSumZ, &upotSum, 0);
      }
      const int rest = (int)(soa2.size() & (vecLength - 1));
      if (rest > 0)
        SoAKernel<newton3, true>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                 y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc,
                                 &virialSumX, &virialSumY, &virialSumZ, &upotSum, rest);

      // horizontally reduce fDacc to sumfD
      const simde__m256d hSumfxfy = simde_mm256_hadd_pd(fxacc, fyacc);
      const simde__m256d hSumfz = simde_mm256_hadd_pd(fzacc, fzacc);

      const simde__m128d hSumfxfyLow = simde_mm256_extractf128_pd(hSumfxfy, 0);
      const simde__m128d hSumfzLow = simde_mm256_extractf128_pd(hSumfz, 0);

      const simde__m128d hSumfxfyHigh = simde_mm256_extractf128_pd(hSumfxfy, 1);
      const simde__m128d hSumfzHigh = simde_mm256_extractf128_pd(hSumfz, 1);

      const simde__m128d sumfxfyVEC = simde_mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
      const simde__m128d sumfzVEC = simde_mm_add_pd(hSumfzLow, hSumfzHigh);

      const double sumfx = sumfxfyVEC[0];
      const double sumfy = sumfxfyVEC[1];
      const double sumfz = simde_mm_cvtsd_f64(sumfzVEC);

      fx1ptr[i] += sumfx;
      fy1ptr[i] += sumfy;
      fz1ptr[i] += sumfz;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      // horizontally reduce virialSumX and virialSumY
      const simde__m256d hSumVirialxy = simde_mm256_hadd_pd(virialSumX, virialSumY);
      const simde__m128d hSumVirialxyLow = simde_mm256_extractf128_pd(hSumVirialxy, 0);
      const simde__m128d hSumVirialxyHigh = simde_mm256_extractf128_pd(hSumVirialxy, 1);
      const simde__m128d hSumVirialxyVec = simde_mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      // horizontally reduce virialSumZ and upotSum
      const simde__m256d hSumVirialzUpot = simde_mm256_hadd_pd(virialSumZ, upotSum);
      const simde__m128d hSumVirialzUpotLow = simde_mm256_extractf128_pd(hSumVirialzUpot, 0);
      const simde__m128d hSumVirialzUpotHigh = simde_mm256_extractf128_pd(hSumVirialzUpot, 1);
      const simde__m128d hSumVirialzUpotVec = simde_mm_add_pd(hSumVirialzUpotHigh, hSumVirialzUpotLow);

      // globals = {virialX, virialY, virialZ, uPot}
      double globals[4];
      simde_mm_store_pd(&globals[0], hSumVirialxyVec);
      simde_mm_store_pd(&globals[2], hSumVirialzUpotVec);

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
  }

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
   * @param upotSum
   * @param rest
   */
  template <bool newton3, bool remainderIsMasked>
  inline void SoAKernel(const size_t j, const simde__m256i ownedStateI, const int64_t *const __restrict ownedStatePtr2,
                        const simde__m256d &x1, const simde__m256d &y1, const simde__m256d &z1,
                        const double *const __restrict x2ptr, const double *const __restrict y2ptr,
                        const double *const __restrict z2ptr, double *const __restrict fx2ptr,
                        double *const __restrict fy2ptr, double *const __restrict fz2ptr,
                        const size_t *const typeID1ptr, const size_t *const typeID2ptr, simde__m256d &fxacc,
                        simde__m256d &fyacc, simde__m256d &fzacc, simde__m256d *virialSumX, simde__m256d *virialSumY,
                        simde__m256d *virialSumZ, simde__m256d *upotSum, const unsigned int rest = 0) {
    simde__m256d epsilon24s = _epsilon24;
    simde__m256d sigmaSquares = _sigmaSquare;
    simde__m256d shift6s = _shift6;
    if (useMixing) {
      // the first argument for set lands in the last bits of the register
      epsilon24s = simde_mm256_set_pd(
          not remainderIsMasked or rest > 3 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 3)) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 2)) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 1)) : 0,
          _PPLibrary->getMixing24Epsilon(*typeID1ptr, *(typeID2ptr + 0)));
      sigmaSquares = simde_mm256_set_pd(
          not remainderIsMasked or rest > 3 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 3)) : 0,
          not remainderIsMasked or rest > 2 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 2)) : 0,
          not remainderIsMasked or rest > 1 ? _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 1)) : 0,
          _PPLibrary->getMixingSigmaSquared(*typeID1ptr, *(typeID2ptr + 0)));
      if constexpr (applyShift) {
        shift6s = simde_mm256_set_pd(
            (not remainderIsMasked or rest > 3) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 3)) : 0,
            (not remainderIsMasked or rest > 2) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 2)) : 0,
            (not remainderIsMasked or rest > 1) ? _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 1)) : 0,
            _PPLibrary->getMixingShift6(*typeID1ptr, *(typeID2ptr + 0)));
      }
    }

    const simde__m256d x2 =
        remainderIsMasked ? simde_mm256_maskload_pd(&x2ptr[j], _masks[rest - 1]) : simde_mm256_loadu_pd(&x2ptr[j]);
    const simde__m256d y2 =
        remainderIsMasked ? simde_mm256_maskload_pd(&y2ptr[j], _masks[rest - 1]) : simde_mm256_loadu_pd(&y2ptr[j]);
    const simde__m256d z2 =
        remainderIsMasked ? simde_mm256_maskload_pd(&z2ptr[j], _masks[rest - 1]) : simde_mm256_loadu_pd(&z2ptr[j]);

    const simde__m256d drx = simde_mm256_sub_pd(x1, x2);
    const simde__m256d dry = simde_mm256_sub_pd(y1, y2);
    const simde__m256d drz = simde_mm256_sub_pd(z1, z2);

    const simde__m256d drx2 = simde_mm256_mul_pd(drx, drx);
    const simde__m256d dry2 = simde_mm256_mul_pd(dry, dry);
    const simde__m256d drz2 = simde_mm256_mul_pd(drz, drz);

    const simde__m256d dr2PART = simde_mm256_add_pd(drx2, dry2);
    const simde__m256d dr2 = simde_mm256_add_pd(dr2PART, drz2);

    // SIMDE_CMP_LE_OS == Less-Equal-then (ordered, signaling)
    // signaling = throw error if NaN is encountered
    // dr2 <= _cutoffsquare ? 0xFFFFFFFFFFFFFFFF : 0
    const simde__m256d cutoffMask = simde_mm256_cmp_pd(dr2, _cutoffsquare, SIMDE_CMP_LE_OS);

    // This requires that dummy is zero (otherwise when loading using a mask the owned state will not be zero)
    const simde__m256i ownedStateJ =
        remainderIsMasked ? simde_mm256_castpd_si256(simde_mm256_maskload_pd(
                                reinterpret_cast<double const *>(&ownedStatePtr2[j]), _masks[rest - 1]))
                          : simde_mm256_loadu_si256(reinterpret_cast<const simde__m256i *>(&ownedStatePtr2[j]));
    // This requires that dummy is the first entry in OwnershipState!
    const simde__m256d dummyMask = simde_mm256_cmp_pd(simde_mm256_castsi256_pd(ownedStateJ), _zero, SIMDE_CMP_NEQ_UQ);
    const simde__m256d cutoffDummyMask = simde_mm256_and_pd(cutoffMask, dummyMask);

    // if everything is masked away return from this function.
    if (simde_mm256_movemask_pd(cutoffDummyMask) == 0) {
      return;
    }

    const simde__m256d invdr2 = simde_mm256_div_pd(_one, dr2);
    const simde__m256d lj2 = simde_mm256_mul_pd(sigmaSquares, invdr2);
    const simde__m256d lj4 = simde_mm256_mul_pd(lj2, lj2);
    const simde__m256d lj6 = simde_mm256_mul_pd(lj2, lj4);
    const simde__m256d lj12 = simde_mm256_mul_pd(lj6, lj6);
    const simde__m256d lj12m6 = simde_mm256_sub_pd(lj12, lj6);
    const simde__m256d lj12m6alj12 = simde_mm256_add_pd(lj12m6, lj12);
    const simde__m256d lj12m6alj12e = simde_mm256_mul_pd(lj12m6alj12, epsilon24s);
    const simde__m256d fac = simde_mm256_mul_pd(lj12m6alj12e, invdr2);

    const simde__m256d facMasked =
        remainderIsMasked
            ? simde_mm256_and_pd(fac, simde_mm256_and_pd(cutoffDummyMask, simde_mm256_castsi256_pd(_masks[rest - 1])))
            : simde_mm256_and_pd(fac, cutoffDummyMask);

    const simde__m256d fx = simde_mm256_mul_pd(drx, facMasked);
    const simde__m256d fy = simde_mm256_mul_pd(dry, facMasked);
    const simde__m256d fz = simde_mm256_mul_pd(drz, facMasked);

    fxacc = simde_mm256_add_pd(fxacc, fx);
    fyacc = simde_mm256_add_pd(fyacc, fy);
    fzacc = simde_mm256_add_pd(fzacc, fz);

    // if newton 3 is used subtract fD from particle j
    if constexpr (newton3) {
      const simde__m256d fx2 =
          remainderIsMasked ? simde_mm256_maskload_pd(&fx2ptr[j], _masks[rest - 1]) : simde_mm256_loadu_pd(&fx2ptr[j]);
      const simde__m256d fy2 =
          remainderIsMasked ? simde_mm256_maskload_pd(&fy2ptr[j], _masks[rest - 1]) : simde_mm256_loadu_pd(&fy2ptr[j]);
      const simde__m256d fz2 =
          remainderIsMasked ? simde_mm256_maskload_pd(&fz2ptr[j], _masks[rest - 1]) : simde_mm256_loadu_pd(&fz2ptr[j]);

      const simde__m256d fx2new = simde_mm256_sub_pd(fx2, fx);
      const simde__m256d fy2new = simde_mm256_sub_pd(fy2, fy);
      const simde__m256d fz2new = simde_mm256_sub_pd(fz2, fz);

      remainderIsMasked ? simde_mm256_maskstore_pd(&fx2ptr[j], _masks[rest - 1], fx2new)
                        : simde_mm256_storeu_pd(&fx2ptr[j], fx2new);
      remainderIsMasked ? simde_mm256_maskstore_pd(&fy2ptr[j], _masks[rest - 1], fy2new)
                        : simde_mm256_storeu_pd(&fy2ptr[j], fy2new);
      remainderIsMasked ? simde_mm256_maskstore_pd(&fz2ptr[j], _masks[rest - 1], fz2new)
                        : simde_mm256_storeu_pd(&fz2ptr[j], fz2new);
    }

    if constexpr (calculateGlobals) {
      // Global Virial
      const simde__m256d virialX = simde_mm256_mul_pd(fx, drx);
      const simde__m256d virialY = simde_mm256_mul_pd(fy, dry);
      const simde__m256d virialZ = simde_mm256_mul_pd(fz, drz);

      // Global Potential
      const simde__m256d upot = wrapperFMA(epsilon24s, lj12m6, shift6s);

      const simde__m256d upotMasked =
          remainderIsMasked ? simde_mm256_and_pd(
                                  upot, simde_mm256_and_pd(cutoffDummyMask, simde_mm256_castsi256_pd(_masks[rest - 1])))
                            : simde_mm256_and_pd(upot, cutoffDummyMask);

      simde__m256d ownedMaskI = simde_mm256_cmp_pd(simde_mm256_castsi256_pd(ownedStateI),
                                                   simde_mm256_castsi256_pd(_ownedStateOwnedMM256i), SIMDE_CMP_EQ_UQ);
      simde__m256d energyFactor = simde_mm256_blendv_pd(_zero, _one, ownedMaskI);
      if constexpr (newton3) {
        simde__m256d ownedMaskJ = simde_mm256_cmp_pd(simde_mm256_castsi256_pd(ownedStateJ),
                                                     simde_mm256_castsi256_pd(_ownedStateOwnedMM256i), SIMDE_CMP_EQ_UQ);
        energyFactor = simde_mm256_add_pd(energyFactor, simde_mm256_blendv_pd(_zero, _one, ownedMaskJ));
      }
      *upotSum = wrapperFMA(energyFactor, upotMasked, *upotSum);
      *virialSumX = wrapperFMA(energyFactor, virialX, *virialSumX);
      *virialSumY = wrapperFMA(energyFactor, virialY, *virialSumY);
      *virialSumZ = wrapperFMA(energyFactor, virialZ, *virialSumZ);
    }
  }

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

    // accumulators
    simde__m256d virialSumX = simde_mm256_setzero_pd();
    simde__m256d virialSumY = simde_mm256_setzero_pd();
    simde__m256d virialSumZ = simde_mm256_setzero_pd();
    simde__m256d upotSum = simde_mm256_setzero_pd();
    simde__m256d fxacc = simde_mm256_setzero_pd();
    simde__m256d fyacc = simde_mm256_setzero_pd();
    simde__m256d fzacc = simde_mm256_setzero_pd();

    // broadcast particle 1
    const auto x1 = simde_mm256_broadcast_sd(&xptr[indexFirst]);
    const auto y1 = simde_mm256_broadcast_sd(&yptr[indexFirst]);
    const auto z1 = simde_mm256_broadcast_sd(&zptr[indexFirst]);
    // ownedStatePtr contains int64_t, so we broadcast these to make an simde__m256i.
    // simde_mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
    simde__m256i ownedStateI = simde_mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr[indexFirst]));

    alignas(64) std::array<double, vecLength> x2tmp{};
    alignas(64) std::array<double, vecLength> y2tmp{};
    alignas(64) std::array<double, vecLength> z2tmp{};
    alignas(64) std::array<double, vecLength> fx2tmp{};
    alignas(64) std::array<double, vecLength> fy2tmp{};
    alignas(64) std::array<double, vecLength> fz2tmp{};
    alignas(64) std::array<size_t, vecLength> typeID2tmp{};
    alignas(64) std::array<autopas::OwnershipState, vecLength> ownedStates2tmp{};

    // load 4 neighbors
    size_t j = 0;
    // Loop over all neighbors as long as we can fill full vectors
    // (until `neighborList.size() - neighborList.size() % vecLength`)
    //
    // If b is a power of 2 the following holds:
    // a & ~(b - 1) == a - (a mod b)
    for (; j < (neighborList.size() & ~(vecLength - 1)); j += vecLength) {
      // SIMDe2 variant:
      // create buffer for 4 interaction particles
      // and fill buffers via gathering
      //      const simde__m256d x2tmp = simde_mm256_i64gather_pd(&xptr[j], _vindex, 1);
      //      const simde__m256d y2tmp = simde_mm256_i64gather_pd(&yptr[j], _vindex, 1);
      //      const simde__m256d z2tmp = simde_mm256_i64gather_pd(&zptr[j], _vindex, 1);
      //      const simde__m256d fx2tmp = simde_mm256_i64gather_pd(&fxptr[j], _vindex, 1);
      //      const simde__m256d fy2tmp = simde_mm256_i64gather_pd(&fyptr[j], _vindex, 1);
      //      const simde__m256d fz2tmp = simde_mm256_i64gather_pd(&fzptr[j], _vindex, 1);
      //      const simde__m256i typeID2tmp = simde_mm256_i64gather_epi64(&typeIDptr[j], _vindex, 1);

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

      SoAKernel<newton3, false>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2tmp.data()), x1, y1, z1,
                                x2tmp.data(), y2tmp.data(), z2tmp.data(), fx2tmp.data(), fy2tmp.data(), fz2tmp.data(),
                                &typeIDptr[indexFirst], typeID2tmp.data(), fxacc, fyacc, fzacc, &virialSumX,
                                &virialSumY, &virialSumZ, &upotSum, 0);

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
      // SIMDe2 variant:
      // create buffer for 4 interaction particles
      // and fill buffers via gathering
      //      TODO: use masked load because there will not be enough data left for the whole gather
      //      const simde__m256d x2tmp = simde_mm256_i64gather_pd(&xptr[j], _vindex, 1);
      //      const simde__m256d y2tmp = simde_mm256_i64gather_pd(&yptr[j], _vindex, 1);
      //      const simde__m256d z2tmp = simde_mm256_i64gather_pd(&zptr[j], _vindex, 1);
      //      const simde__m256d fx2tmp = simde_mm256_i64gather_pd(&fxptr[j], _vindex, 1);
      //      const simde__m256d fy2tmp = simde_mm256_i64gather_pd(&fyptr[j], _vindex, 1);
      //      const simde__m256d fz2tmp = simde_mm256_i64gather_pd(&fzptr[j], _vindex, 1);
      //      const simde__m256d typeID2tmp = simde_mm256_i64gather_pd(&typeIDptr[j], _vindex, 1);

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

      SoAKernel<newton3, true>(0, ownedStateI, reinterpret_cast<const int64_t *>(ownedStates2tmp.data()), x1, y1, z1,
                               x2tmp.data(), y2tmp.data(), z2tmp.data(), fx2tmp.data(), fy2tmp.data(), fz2tmp.data(),
                               &typeIDptr[indexFirst], typeID2tmp.data(), fxacc, fyacc, fzacc, &virialSumX, &virialSumY,
                               &virialSumZ, &upotSum, rest);

      if constexpr (newton3) {
        for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
          fxptr[neighborList[j + vecIndex]] = fx2tmp[vecIndex];
          fyptr[neighborList[j + vecIndex]] = fy2tmp[vecIndex];
          fzptr[neighborList[j + vecIndex]] = fz2tmp[vecIndex];
        }
      }
    }

    // horizontally reduce fDacc to sumfD
    const simde__m256d hSumfxfy = simde_mm256_hadd_pd(fxacc, fyacc);
    const simde__m256d hSumfz = simde_mm256_hadd_pd(fzacc, fzacc);

    const simde__m128d hSumfxfyLow = simde_mm256_extractf128_pd(hSumfxfy, 0);
    const simde__m128d hSumfzLow = simde_mm256_extractf128_pd(hSumfz, 0);

    const simde__m128d hSumfxfyHigh = simde_mm256_extractf128_pd(hSumfxfy, 1);
    const simde__m128d hSumfzHigh = simde_mm256_extractf128_pd(hSumfz, 1);

    const simde__m128d sumfxfyVEC = simde_mm_add_pd(hSumfxfyLow, hSumfxfyHigh);
    const simde__m128d sumfzVEC = simde_mm_add_pd(hSumfzLow, hSumfzHigh);

    const double sumfx = sumfxfyVEC[0];
    const double sumfy = sumfxfyVEC[1];
    const double sumfz = simde_mm_cvtsd_f64(sumfzVEC);

    fxptr[indexFirst] += sumfx;
    fyptr[indexFirst] += sumfy;
    fzptr[indexFirst] += sumfz;

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      // horizontally reduce virialSumX and virialSumY
      const simde__m256d hSumVirialxy = simde_mm256_hadd_pd(virialSumX, virialSumY);
      const simde__m128d hSumVirialxyLow = simde_mm256_extractf128_pd(hSumVirialxy, 0);
      const simde__m128d hSumVirialxyHigh = simde_mm256_extractf128_pd(hSumVirialxy, 1);
      const simde__m128d hSumVirialxyVec = simde_mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      // horizontally reduce virialSumZ and upotSum
      const simde__m256d hSumVirialzUpot = simde_mm256_hadd_pd(virialSumZ, upotSum);
      const simde__m128d hSumVirialzUpotLow = simde_mm256_extractf128_pd(hSumVirialzUpot, 0);
      const simde__m128d hSumVirialzUpotHigh = simde_mm256_extractf128_pd(hSumVirialzUpot, 1);
      const simde__m128d hSumVirialzUpotVec = simde_mm_add_pd(hSumVirialzUpotHigh, hSumVirialzUpotLow);

      // globals = {virialX, virialY, virialZ, uPot}
      double globals[4];
      simde_mm_store_pd(&globals[0], hSumVirialxyVec);
      simde_mm_store_pd(&globals[2], hSumVirialzUpotVec);

      double factor = 1.;
      // we assume newton3 to be enabled in this function call, thus we multiply by two if the value of newton3 is
      // false, since for newton3 disabled we divide by two later on.
      factor *= newton3 ? .5 : 1.;
      // In case we have a non-cell-wise owned state, we have multiplied everything by two, so we divide it by 2 again.
      _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
      _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
      _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
      _aosThreadData[threadnum].upotSum += globals[3] * factor;
    }
    // interact with i with 4 neighbors
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
  void initTraversal() final {
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
  void endTraversal(bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (_postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }

    if (calculateGlobals) {
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        _upotSum += _aosThreadData[i].upotSum;
        _virialSum += _aosThreadData[i].virialSum;
      }
      if (not newton3) {
        // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
        // here.
        _upotSum *= 0.5;
        _virialSum *= 0.5;
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
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get upot even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException("Cannot get upot, because endTraversal was not called.");
    }
    return _upotSum;
  }

  /**
   * Get the virial
   * @return the virial
   */
  double getVirial() {
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
   * @param sigmaSquare
   */
  void setParticleProperties(double epsilon24, double sigmaSquare) {
    _epsilon24 = simde_mm256_set1_pd(epsilon24);
    _sigmaSquare = simde_mm256_set1_pd(sigmaSquare);
    if constexpr (applyShift) {
      _shift6 = simde_mm256_set1_pd(
          ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, _cutoffsquare[0]));
    } else {
      _shift6 = simde_mm256_setzero_pd();
    }

    _epsilon24AoS = epsilon24;
    _sigmaSquareAoS = sigmaSquare;
    if constexpr (applyShift) {
      _shift6AoS = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, _cutoffsquareAoS);
    } else {
      _shift6AoS = 0.;
    }
  }

 private:
  /**
   * Wrapper function for FMA. If FMA is not supported it executes first the multiplication then the addition.
   * @param factorA
   * @param factorB
   * @param summandC
   * @return A * B + C
   */
  inline simde__m256d wrapperFMA(const simde__m256d &factorA, const simde__m256d &factorB,
                                 const simde__m256d &summandC) {
    // const simde__m256d tmp = simde_mm256_mul_pd(factorA, factorB);
    // return simde_mm256_add_pd(summandC, tmp);
    return simde_mm256_fmadd_pd(factorA, factorB, summandC);
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

  const simde__m256d _zero{simde_mm256_set1_pd(0.)};
  const simde__m256d _one{simde_mm256_set1_pd(1.)};
  const simde__m256i _vindex = simde_mm256_set_epi64x(0, 1, 3, 4);
  const simde__m256i _masks[3]{
      simde_mm256_set_epi64x(0, 0, 0, -1),
      simde_mm256_set_epi64x(0, 0, -1, -1),
      simde_mm256_set_epi64x(0, -1, -1, -1),
  };
  const simde__m256i _ownedStateDummyMM256i{0x0};
  const simde__m256i _ownedStateOwnedMM256i{
      simde_mm256_set1_epi64x(static_cast<int64_t>(autopas::OwnershipState::owned))};
  const simde__m256d _cutoffsquare{};
  simde__m256d _shift6 = simde_mm256_setzero_pd();
  simde__m256d _epsilon24{};
  simde__m256d _sigmaSquare{};

  const double _cutoffsquareAoS = 0;
  double _epsilon24AoS, _sigmaSquareAoS, _shift6AoS = 0;

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
  // MUST be power of 2 because some optimizations make this assumption
  constexpr static size_t vecLength = 4;
};
}  // namespace mdLib
