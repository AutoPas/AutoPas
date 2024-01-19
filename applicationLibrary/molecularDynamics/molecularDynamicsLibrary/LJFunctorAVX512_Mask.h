/**
 * @file LJFunctorAVX512_Mask.h
 *
 * @date 17/01/2024
 * @author S. Newcome
 */
#pragma once
#ifndef __AVX512F__
#pragma message "Requested to compile LJFunctorAVX512 but AVX512 is not available!"
#else
#include <immintrin.h>
#endif

#include <array>

#include "ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This Version is implemented using AVX512 intrinsics using a mask to not add force contributions from pairs of molecules beyond the cutoff.
 * @tparam Particle The type of particle.
 * @tparam applyShift Switch for the lj potential to be truncated shifted.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true>
class LJFunctorAVX512_Mask
    : public autopas::Functor<
          Particle, LJFunctorAVX512_Mask<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorAVX512_Mask() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy unused, only there to make the signature different from the public constructor.
   */
  explicit LJFunctorAVX512_Mask(double cutoff, void * /*dummy*/)
#ifdef __AVX512F__
      : autopas::Functor<
            Particle, LJFunctorAVX512_Mask<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff),
        _cutoffSquared{_mm512_set1_pd(cutoff * cutoff)},
        _cutoffSquaredAoS(cutoff * cutoff),
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData(),
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
  }
#else
      : autopas::Functor<
            Particle, LJFunctorAVX512_Mask<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(
            cutoff) {
    autopas::utils::ExceptionHandler::exception("LJFunctorAVX512_Mask was compiled without AVX512 support!");
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
  explicit LJFunctorAVX512_Mask(double cutoff) : LJFunctorAVX512_Mask(cutoff, nullptr) {
//    static_assert(not useMixing,
//                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
//                  "mixing to false.");
    static_assert("LJFunctorAVX512_Mask requires a PPL");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like sigma, epsilon and shift.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit LJFunctorAVX512_Mask(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctorAVX512_Mask(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  inline void AoSFunctor(Particle &particleA, Particle &particleB, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (particleA.isDummy() or particleB.isDummy()) {
      return;
    }
    auto sigmaSquared = _sigmaSquaredAoS;
    auto epsilon24 = _epsilon24AoS;
    auto shift6 = _shift6AoS;
    if constexpr (useMixing) {
      sigmaSquared = _PPLibrary->getMixingSigmaSquared(particleA.getTypeId(), particleB.getTypeId());
      epsilon24 = _PPLibrary->getMixing24Epsilon(particleA.getTypeId(), particleB.getTypeId());
      if constexpr (applyShift) {
        shift6 = _PPLibrary->getMixingShift6(particleA.getTypeId(), particleB.getTypeId());
      }
    }
    auto displacement = particleA.getR() - particleB.getR();
    double distanceSquared = autopas::utils::ArrayMath::dot(displacement, displacement);

    if (distanceSquared > _cutoffSquaredAoS) {
      return;
    }

    double invdr2 = 1. / distanceSquared;
    double lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    auto f = displacement * fac;
    particleA.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      particleB.subF(f);
    }
    if (calculateGlobals) {
      auto virial = displacement * f;
      // Here we calculate either the potential energy * 6 or the potential energy * 12.
      // For newton3, this potential energy contribution is distributed evenly to the two molecules.
      // For non-newton3, the full potential energy is added to the one molecule.
      // The division by 6 is handled in endTraversal, as well as the division by two needed if newton3 is not used.
      double potentialEnergy6 = epsilon24 * lj12m6 + shift6;

      const int threadnum = autopas::autopas_get_thread_num();
      if (particleA.isOwned()) {
        if (newton3) {
          _aosThreadData[threadnum].potentialEnergySumN3 += potentialEnergy6 * 0.5;
          _aosThreadData[threadnum].virialSumN3 += virial * 0.5;
        } else {
          // for non-newton3 the division is in the post-processing step.
          _aosThreadData[threadnum].potentialEnergySumNoN3 += potentialEnergy6;
          _aosThreadData[threadnum].virialSumNoN3 += virial;
        }
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and particleB.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySumN3 += potentialEnergy6 * 0.5;
        _aosThreadData[threadnum].virialSumN3 += virial * 0.5;
      }
    }
  }

  /**
   * @copydoc autopas::Functor::SoAFunctorSingle()
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
   * @copydoc autopas::Functor::SoAFunctorPair()
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
   * Implementation of SoAFunctorSingle with newton3 as a template
   * @tparam newton3
   * @param soa
   */
  template <bool newton3>
  inline void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
#ifdef __AVX512F__
    if (soa.size() == 0) return;

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

    // Set up globals accumulators. Each double is then reduced to a single double at the end of this function.
    __m512d virialAccX = _mm512_setzero_pd();
    __m512d virialAccY = _mm512_setzero_pd();
    __m512d virialAccZ = _mm512_setzero_pd();
    __m512d potentialEnergyAcc = _mm512_setzero_pd();

    // reverse outer loop such that the inner loop always begins at aligned array start
    // typecast to detect underflow
    for (size_t i = soa.size() - 1; (long)i >= 0; --i) {
      if (ownedStatePtr[i] == autopas::OwnershipState::dummy) {
        // If the i-th particle is a dummy, skip this loop iteration.
        continue;
      }

      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");
      // ownedStatePtr contains int64_t, so we broadcast these to make an __m512i.
      // _mm512_set1_epi64 broadcasts a 64-bit integer, we use this instruction to have vecLength values!
      __m512i ownedStateI = _mm512_set1_epi64(static_cast<int64_t>(ownedStatePtr[i]));
      const __mmask8 isOwnedMaskI = _mm512_cmp_epi64_mask(ownedStateI, _ownedStateOwnedMM256i, _MM_CMPINT_EQ);

      __m512d fxacc = _mm512_setzero_pd();
      __m512d fyacc = _mm512_setzero_pd();
      __m512d fzacc = _mm512_setzero_pd();

      // Load positions of molecule i into full vector
      const __m512d x1 = _mm512_set1_pd(&xptr[i]);
      const __m512d y1 = _mm512_set1_pd(&yptr[i]);
      const __m512d z1 = _mm512_set1_pd(&zptr[i]);

      size_t j = 0;
      // floor soa numParticles to multiple of vecLength
      // If b is a power of 2 the following holds:
      // a & ~(b -1) == a - (a mod b)
      for (; j < (i & ~(vecLength - 1)); j += vecLength) {
        SoAKernel<true, false, false>(j, isOwnedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr,
                               zptr, fxptr, fyptr, fzptr, &typeIDptr[i], typeIDptr, fxacc, fyacc, fzacc, &virialAccX,
                               &virialAccY, &virialAccZ, &potentialEnergyAcc, 0);
      }
      // If b is a power of 2 the following holds:
      // a & (b -1) == a mod b
      const int remainder = (int)(i & (vecLength - 1));
      if (remainder > 0) {
        SoAKernel<true, true, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr), x1, y1, z1, xptr, yptr,
                              zptr, fxptr, fyptr, fzptr, &typeIDptr[i], typeIDptr, fxacc, fyacc, fzacc, &virialAccX,
                              &virialAccY, &virialAccZ, &potentialEnergyAcc, remainder);
      }

      // We accumulated the forces acting on molecule i across a full vector. We horizontally reduce this, fDacc, to sumfD
      const double sumfx = _mm512_reduce_add_pd(fxacc);
      const double sumfy = _mm512_reduce_add_pd(fyacc);
      const double sumfz = _mm512_reduce_add_pd(fzacc);

      fxptr[i] += sumfx;
      fyptr[i] += sumfy;
      fzptr[i] += sumfz;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      // Reduce global accumulators
      const double virialSumX = _mm512_reduce_add_pd(virialAccX);
      const double virialSumY = _mm512_reduce_add_pd(virialAccY);
      const double virialSumZ = _mm512_reduce_add_pd(virialAccZ);
      const double potentialEnergySum = _mm512_reduce_add_pd(potentialEnergyAcc);

      // we assume newton3 to be enabled in this function call, thus we multiply by two if the value of newton3 is
      // false, since for newton3 disabled we divide by two later on.
      if constexpr (newton3) {
        _aosThreadData[threadnum].virialSumN3[0] += virialSumX * 0.5;
        _aosThreadData[threadnum].virialSumN3[1] += virialSumY * 0.5;
        _aosThreadData[threadnum].virialSumN3[2] += virialSumZ * 0.5;
        _aosThreadData[threadnum].potentialEnergySumN3 += potentialEnergySum * 0.5;
      } else {
        _aosThreadData[threadnum].virialSumNoN3[0] += virialSumX;
        _aosThreadData[threadnum].virialSumNoN3[1] += virialSumY;
        _aosThreadData[threadnum].virialSumNoN3[2] += virialSumZ;
        _aosThreadData[threadnum].potentialEnergySumNoN3 += potentialEnergySum;
      }
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

    // Set up globals accumulators. Each double is then reduced to a single double at the end of this function.
    __m512d virialAccX = _mm512_setzero_pd();
    __m512d virialAccY = _mm512_setzero_pd();
    __m512d virialAccZ = _mm512_setzero_pd();
    __m512d potentialEnergyAcc = _mm512_setzero_pd();

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
        // If the i-th particle is a dummy, skip this loop iteration.
        continue;
      }

      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");
      // ownedStatePtr1 contains int64_t, so we broadcast these to make an __m512i.
      // _mm512_set1_epi64 broadcasts a 64-bit integer, we use this instruction to have vecLength values!
      __m512i ownedStateI = _mm512_set1_epi64(static_cast<int64_t>(ownedStatePtr1[i]));
      const __mmask8 isOwnedMaskI = _mm512_cmp_epi64_mask(ownedStateI, _ownedStateOwnedMM256i, _MM_CMPINT_EQ);

      __m512d fxacc = _mm512_setzero_pd();
      __m512d fyacc = _mm512_setzero_pd();
      __m512d fzacc = _mm512_setzero_pd();

      // Load positions of molecule i into full vector
      const __m512d x1 = _mm512_set1_pd(&x1ptr[i]);
      const __m512d y1 = _mm512_set1_pd(&y1ptr[i]);
      const __m512d z1 = _mm512_set1_pd(&z1ptr[i]);

      // floor soa2 numParticles to multiple of vecLength
      size_t j = 0;
      for (; j < (soa2.size() & ~(vecLength - 1)); j += vecLength) {
        SoAKernel<newton3, false, false>(j, ownedStateI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                  y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc,
                                  &virialAccX, &virialAccY, &virialAccZ, &potentialEnergyAcc, 0);
      }
      const int rest = (int)(soa2.size() & (vecLength - 1));
      if (rest > 0)
        SoAKernel<newton3, true, false>(j, isOwnedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, x2ptr,
                                 y2ptr, z2ptr, fx2ptr, fy2ptr, fz2ptr, typeID1ptr, typeID2ptr, fxacc, fyacc, fzacc,
                                 &virialAccX, &virialAccY, &virialAccZ, &potentialEnergyAcc, rest);

      // We accumulated the forces acting on molecule i across a full vector. We horizontally reduce this, fDacc, to sumfD
      const double sumfx = _mm512_reduce_add_pd(fxacc);
      const double sumfy = _mm512_reduce_add_pd(fyacc);
      const double sumfz = _mm512_reduce_add_pd(fzacc);

      fx1ptr[i] += sumfx;
      fy1ptr[i] += sumfy;
      fz1ptr[i] += sumfz;
    }

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      // Reduce global accumulators
      const double virialSumX = _mm512_reduce_add_pd(virialAccX);
      const double virialSumY = _mm512_reduce_add_pd(virialAccY);
      const double virialSumZ = _mm512_reduce_add_pd(virialAccZ);
      const double potentialEnergySum = _mm512_reduce_add_pd(potentialEnergyAcc);

      // we assume newton3 to be enabled in this function call, thus we multiply by two if the value of newton3 is
      // false, since for newton3 disabled we divide by two later on.
      if constexpr (newton3) {
        _aosThreadData[threadnum].virialSumN3[0] += virialSumX * 0.5;
        _aosThreadData[threadnum].virialSumN3[1] += virialSumY * 0.5;
        _aosThreadData[threadnum].virialSumN3[2] += virialSumZ * 0.5;
        _aosThreadData[threadnum].potentialEnergySumN3 += potentialEnergySum * 0.5;
      } else {
        _aosThreadData[threadnum].virialSumNoN3[0] += virialSumX;
        _aosThreadData[threadnum].virialSumNoN3[1] += virialSumY;
        _aosThreadData[threadnum].virialSumNoN3[2] += virialSumZ;
        _aosThreadData[threadnum].potentialEnergySumNoN3 += potentialEnergySum;
      }
    }
#endif
  }

#ifdef __AVX512F__
  /**
   * Actual inner kernel of the SoAFunctors.
   *
   * @tparam newton3
   * @tparam remainderIsMasked If false the full vector length is used. Otherwise the last entries are masked away
   * depending on the argument "rest".
   * @tparam isVerlet true if Kernel is called from SoAFunctorVerlet. If true, SoA index "j" is not used. Instead a __m512i
   * of neighborIndices is used in combination with gather & scatter instructions. The kernel does not check that a non-nullptr pointer is used.
   *
   * @param j
   * @param baseEnergyFactor Either 0 or 1 depending on if molecule 1 is owned. Used in globals calculation
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
   * @param remainder
   */
  template <bool newton3, bool remainderCase, bool isVerlet>
  inline void SoAKernel(const size_t j, const __mmask8 *isOwnedMask1, const int64_t *const __restrict ownedStatePtr2,
                        const __m512d &x1, const __m512d &y1, const __m512d &z1, const double *const __restrict x2ptr,
                        const double *const __restrict y2ptr, const double *const __restrict z2ptr,
                        double *const __restrict fx2ptr, double *const __restrict fy2ptr,
                        double *const __restrict fz2ptr, const double *const __restrict mixingPtr, const size_t *const typeID2ptr,
                        __m512d &fxacc, __m512d &fyacc, __m512d &fzacc, __m512d *virialSumX, __m512d *virialSumY,
                        __m512d *virialSumZ, __m512d *potentialEnergySum, const unsigned int remainder = 0, __m512i *neighborIndices = nullptr) {

    const __mmask8 remainderMask = remainderCase ? _remainderMasks[8-remainder] : __mmask8(255);

    // todo: can we convert this ugly logic into a lambda function with zero overhead

    const __m512d x2 = remainderCase ? (isVerlet ?
                                                 _mm512_mask_i64gather_pd(_zero, remainderMask, *neighborIndices, &x2ptr[0], 8) :
                                                _mm512_maskz_loadu_pd(remainderMask, &x2ptr[j])) :
                                     (isVerlet ?
                                               _mm512_i64gather_pd(*neighborIndices, &x2ptr[0], 8) :
                                               _mm512_loadu_pd(&x2ptr[j]));

    const __m512d y2 = remainderCase ? (isVerlet ?
                                                 _mm512_mask_i64gather_pd(_zero, remainderMask, *neighborIndices, &y2ptr[0], 8) :
                                                 _mm512_maskz_loadu_pd(remainderMask, &y2ptr[j])) :
                                     (isVerlet ?
                                               _mm512_i64gather_pd(*neighborIndices, &y2ptr[0], 8) :
                                               _mm512_loadu_pd(&y2ptr[j]));

    const __m512d z2 = remainderCase ? (isVerlet ?
                                                 _mm512_mask_i64gather_pd(_zero, remainderMask, *neighborIndices, &z2ptr[0], 8) :
                                                 _mm512_maskz_loadu_pd(remainderMask, &z2ptr[j])) :
                                     (isVerlet ?
                                               _mm512_i64gather_pd(*neighborIndices, &z2ptr[0], 8) :
                                               _mm512_loadu_pd(&z2ptr[j]));

    const __m512d displacementX = _mm512_sub_pd(x1, x2);
    const __m512d displacementY = _mm512_sub_pd(y1, y2);
    const __m512d displacementZ = _mm512_sub_pd(z1, z2);

    const __m512d distanceSquaredX = _mm512_mul_pd(displacementX, displacementX);
    const __m512d distanceSquaredY = _mm512_mul_pd(displacementY, displacementY);
    const __m512d distanceSquaredZ = _mm512_mul_pd(displacementZ, displacementZ);

    const __m512d distanceSquared =
        _mm512_add_pd(distanceSquaredX, _mm512_add_pd(distanceSquaredY, distanceSquaredZ));

    // mask = distanceSquared < _cutoffSquared

    // _CMP_LE_OS == Less-Equal-then (ordered, signaling)
    // ordered = gives "false" if NaN is encountered
    // signaling = throw error if NaN is encountered
    const __mmask8 cutoffMask = _mm512_cmp_pd_mask(distanceSquared, _cutoffSquared, _CMP_LE_OS);

    const __m512i ownershipState2 = remainderCase ? (isVerlet ?
                                                 _mm512_mask_i64gather_epi64(_zero, remainderMask, *neighborIndices, &ownedStatePtr2[0], 8) :
                                                 _mm512_maskz_loadu_epi64(remainderMask, &ownedStatePtr2[j])) :
                                     (isVerlet ?
                                               _mm512_i64gather_epi64(*neighborIndices, &ownedStatePtr2[0], 8) :
                                               _mm512_loadu_epi64(&ownedStatePtr2[j]));

    // isNotDummy = ownershipStateB != dummy
    const __mmask8 isNotDummy = _mm512_cmp_epi64_mask(ownershipState2, _ownedStateDummyMM512i, _MM_CMPINT_NE);

    // combinedMask = mask AND isNotDummy
    const __mmask8 combinedMask = _mm512_kand(cutoffMask, isNotDummy);

    // if everything is masked away return from this function.
    if (combinedMask == 0) {
      return;
    }

    // gather mixing parameters.
    // Uses base address of mixingPtr (+0/1/2), gathered data is offset by siteTypeIndicesScaled x 64 bits x 8.
    // We need scaled site-types due to the way the mixing data is stored - todo change this
    const __m512i typeIndices = remainderCase ? (isVerlet ?
                                                              _mm512_mask_i64gather_epi64(_zero, remainderMask, *neighborIndices, &typeID2ptr[0], 8) :
                                                              _mm512_maskz_loadu_epi64(remainderMask, &typeID2ptr[j])) :
                                                  (isVerlet ?
                                                            _mm512_i64gather_epi64(*neighborIndices, &typeID2ptr[0], 8) :
                                                            _mm512_loadu_epi64(&typeID2ptr[j]));

    const __m512i typeIndicesScaled = _mm512_mullox_epi64(typeIndices, _three);

    const __m512d epsilon24 = remainderCase ? _mm512_mask_i64gather_pd(_zero, remainderMask, typeIndicesScaled, mixingPtr, 8) : _mm512_i64gather_pd(typeIndicesScaled, mixingPtr, 8);
    const __m512d sigmaSquared = remainderCase ? _mm512_mask_i64gather_pd(_zero, remainderMask, typeIndicesScaled, mixingPtr+1, 8) : _mm512_i64gather_pd(typeIndicesScaled, mixingPtr+1, 8);
    const __m512d shift6 = remainderCase ? _mm512_mask_i64gather_pd(_zero, remainderMask, typeIndicesScaled, mixingPtr+2, 8) : _mm512_i64gather_pd(typeIndicesScaled, mixingPtr+2, 8);

    const __m512d invDistSquared = _mm512_div_pd(_one, distanceSquared);
    const __m512d lj2 = _mm512_mul_pd(sigmaSquared, invDistSquared); // = (sigma/dist)^2
    const __m512d lj6 = _mm512_mul_pd(_mm512_mul_pd(lj2, lj2), lj2); // = (sigma/dist)^6
    const __m512d lj12 = _mm512_mul_pd(lj6, lj6); // = (sigma/dist)^12
    const __m512d lj12m6 = _mm512_sub_pd(lj12, lj6); // = (sigma/dist)^12 - (sigma/dist)^6

    const __m512d twolj12m6 = _mm512_add_pd(lj12m6, lj12); // = 2 * (sigma/dist)^12 - (sigma/dist)^6
    const __m512d lj12m6alj12mEps24 = _mm512_mul_pd(twolj12m6, epsilon24); // = 24 * epsilon * ( twolj12m6 )
    const __m512d scalar = _mm512_mul_pd(lj12m6alj12mEps24, invDistSquared); // = LJ scalar to be multiplied to displacement

    // Determine forces and mask out molecules beyond cutoff or are dummies
    const __m512d forceX = _mm512_maskz_mul_pd((combinedMask), scalar, displacementX);
    const __m512d forceY = _mm512_maskz_mul_pd((combinedMask), scalar, displacementY);
    const __m512d forceZ = _mm512_maskz_mul_pd((combinedMask), scalar, displacementZ);

    fxacc = remainderCase ? _mm512_mask_add_pd(fxacc, remainderMask, fxacc, forceX) : _mm512_add_pd( fxacc, forceX);
    fyacc = remainderCase ? _mm512_mask_add_pd(fyacc, remainderMask, fyacc, forceY) : _mm512_add_pd( fyacc, forceY);
    fzacc = remainderCase ? _mm512_mask_add_pd(fzacc, remainderMask, fzacc, forceZ) : _mm512_add_pd( fzacc, forceZ);


    // Newton 3 optimization upon particle js
    if constexpr (newton3) {
      if constexpr (remainderCase) {
        if constexpr (isVerlet) {
          const __m512d forceSumJX = _mm512_mask_i64gather_pd(_zero, remainderMask, *neighborIndices, &fx2ptr[0], 8);
          const __m512d newForceSumJX = _mm512_sub_pd(forceSumJX, forceX);
          _mm512_mask_i64scatter_pd(&fx2ptr[j], remainderMask, *neighborIndices, newForceSumJX, 8);

          const __m512d forceSumJY = _mm512_mask_i64gather_pd(_zero, remainderMask, *neighborIndices, &fy2ptr[0], 8);
          const __m512d newForceSumJY = _mm512_sub_pd(forceSumJY, forceY);
          _mm512_mask_i64scatter_pd(&fy2ptr[j], remainderMask, *neighborIndices, newForceSumJY, 8);

          const __m512d forceSumJZ = _mm512_mask_i64gather_pd(_zero, remainderMask, *neighborIndices, &fz2ptr[0], 8);
          const __m512d newForceSumJZ = _mm512_sub_pd(forceSumJZ, forceZ);
          _mm512_mask_i64scatter_pd(&fz2ptr[j], remainderMask, *neighborIndices, newForceSumJZ, 8);
        } else {
          const __m512d forceSumJX = _mm512_maskz_loadu_pd(remainderMask, &fx2ptr[j]);
          const __m512d newForceSumJX = _mm512_sub_pd(forceSumJX, forceX);
          _mm512_mask_storeu_pd(&fx2ptr[j], remainderMask, newForceSumJX);

          const __m512d forceSumJY = _mm512_maskz_loadu_pd(remainderMask, &fy2ptr[j]);
          const __m512d newForceSumJY = _mm512_sub_pd(forceSumJY, forceY);
          _mm512_mask_storeu_pd(&fy2ptr[j], remainderMask, newForceSumJY);

          const __m512d forceSumJZ = _mm512_maskz_loadu_pd(remainderMask, &fz2ptr[j]);
          const __m512d newForceSumJZ = _mm512_sub_pd(forceSumJZ, forceZ);
          _mm512_mask_storeu_pd(&fz2ptr[j], remainderMask, newForceSumJZ);
        }
      } else {
        if constexpr (isVerlet) {
          const __m512d forceSumJX = _mm512_i64gather_pd(*neighborIndices, &fx2ptr[0], 8);
          const __m512d newForceSumJX = _mm512_sub_pd(forceSumJX, forceX);
          _mm512_i64scatter_pd(&fx2ptr[j], *neighborIndices, newForceSumJX, 8);

          const __m512d forceSumJY = _mm512_i64gather_pd(*neighborIndices, &fy2ptr[0], 8);
          const __m512d newForceSumJY = _mm512_sub_pd(forceSumJY, forceY);
          _mm512_i64scatter_pd(&fy2ptr[j], *neighborIndices, newForceSumJY, 8);

          const __m512d forceSumJZ = _mm512_i64gather_pd(*neighborIndices, &fz2ptr[0], 8);
          const __m512d newForceSumJZ = _mm512_sub_pd(forceSumJZ, forceZ);
          _mm512_i64scatter_pd(&fz2ptr[j], *neighborIndices, newForceSumJZ, 8);
        } else {
          const __m512d forceSumJX = _mm512_loadu_pd(&fx2ptr[j]);
          const __m512d newForceSumJX = _mm512_sub_pd(forceSumJX, forceX);
          _mm512_storeu_pd(&fx2ptr[j], newForceSumJX);

          const __m512d forceSumJY = _mm512_loadu_pd(&fy2ptr[j]);
          const __m512d newForceSumJY = _mm512_sub_pd(forceSumJY, forceY);
          _mm512_storeu_pd(&fy2ptr[j], newForceSumJY);

          const __m512d forceSumJZ = _mm512_loadu_pd(&fz2ptr[j]);
          const __m512d newForceSumJZ = _mm512_sub_pd(forceSumJZ, forceZ);
          _mm512_storeu_pd(&fz2ptr[j], newForceSumJZ);
        }
      }

    }

    if constexpr (calculateGlobals) {
      // Global Virial
      const __m512d virialX = remainderCase ? _mm512_maskz_mul_pd(remainderMask, forceX, displacementX) : _mm512_mul_pd(forceX, displacementX);
      const __m512d virialY = remainderCase ? _mm512_maskz_mul_pd(remainderMask, forceY, displacementY) : _mm512_mul_pd(forceY, displacementY);
      const __m512d virialZ = remainderCase ? _mm512_maskz_mul_pd(remainderMask, forceZ, displacementZ) : _mm512_mul_pd(forceZ, displacementZ);

      // Global Potential
      const __m512d potentialEnergy = remainderCase ? _mm512_maskz_fmadd_pd(remainderMask, epsilon24, lj12m6, shift6) : _mm512_fmadd_pd(epsilon24, lj12m6, shift6);

      // Add above to accumulators iff molecule 1 is owned, using isOwnedMask1
      *virialSumX = _mm512_mask_add_pd(*virialSumX, *isOwnedMask1, *virialSumX, virialX);
      *virialSumY = _mm512_mask_add_pd(*virialSumY, *isOwnedMask1, *virialSumY, virialY);
      *virialSumZ = _mm512_mask_add_pd(*virialSumZ, *isOwnedMask1, *virialSumZ, virialZ);
      *potentialEnergySum = _mm512_mask_add_pd(*potentialEnergySum, *isOwnedMask1, *potentialEnergySum, potentialEnergy);

      // If newton3, do the same again iff molecules in 2 are owned
      if constexpr (newton3) {
        const __mmask8 isOwnedMask2 = _mm512_cmp_epi64_mask(ownershipState2, _ownedStateOwnedMM256i, _MM_CMPINT_EQ);

        *virialSumX = _mm512_mask_add_pd(*virialSumX, isOwnedMask2, *virialSumX, virialX);
        *virialSumY = _mm512_mask_add_pd(*virialSumY, isOwnedMask2, *virialSumY, virialY);
        *virialSumZ = _mm512_mask_add_pd(*virialSumZ, isOwnedMask2, *virialSumZ, virialZ);
        *potentialEnergySum = _mm512_mask_add_pd(*potentialEnergySum, isOwnedMask2, *potentialEnergySum, potentialEnergy);
      }

    }
  }
#endif

 public:
  // clang-format off
  /**
   * @copydoc autopas::Functor::SoAFunctorVerlet()
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
#ifdef __AVX512F__
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
    __m256d virialSumX = _mm256_setzero_pd();
    __m256d virialSumY = _mm256_setzero_pd();
    __m256d virialSumZ = _mm256_setzero_pd();
    __m256d potentialEnergySum = _mm256_setzero_pd();
    __m256d fxacc = _mm256_setzero_pd();
    __m256d fyacc = _mm256_setzero_pd();
    __m256d fzacc = _mm256_setzero_pd();

    // broadcast particle 1
    const auto x1 = _mm256_broadcast_sd(&xptr[indexFirst]);
    const auto y1 = _mm256_broadcast_sd(&yptr[indexFirst]);
    const auto z1 = _mm256_broadcast_sd(&zptr[indexFirst]);
    // ownedStatePtr contains int64_t, so we broadcast these to make an __m256i.
    // _mm256_set1_epi64x broadcasts a 64-bit integer, we use this instruction to have 4 values!
    __m256i ownedStateI = _mm256_set1_epi64x(static_cast<int64_t>(ownedStatePtr[indexFirst]));


    // load vecLength neighbors
    size_t j = 0;
    // Loop over all neighbors as long as we can fill full vectors
    // (until `neighborList.size() - neighborList.size() % vecLength`)
    //
    // If b is a power of 2 the following holds:
    // a & ~(b - 1) == a - (a mod b)
    for (; j < (neighborList.size() & ~(vecLength - 1)); j += vecLength) {
      const __m512i neighborIndices = _mm512_loadu_epi64(&neighborList[j]);
      SoAKernel<newton3, false, true>(0, ownedStateI, ownedStatePtr, x1, y1, z1,
                                xptr, yptr, zptr, fxptr, fyptr, fzptr,
                                &typeIDptr[indexFirst], typeIDptr, fxacc, fyacc, fzacc, &virialSumX,
                                &virialSumY, &virialSumZ, &potentialEnergySum, 0, neighborIndices);
    }
    // Remainder loop
    // If b is a power of 2 the following holds:
    // a & (b - 1) == a mod b
    const auto remainder = static_cast<int>(neighborList.size() & (vecLength - 1));
    if (remainder > 0) {
      const __mmask8 remainderMask = _remainderMasks[8-remainder];
      const __m512i neighborIndices = _mm512_maskz_loadu_epi64(remainderMask, &neighborList[j]);

      SoAKernel<newton3, true>(0, ownedStateI, ownedStatePtr, x1, y1, z1,
                               xptr, yptr, zptr, fxptr, fyptr, fzptr,
                               &typeIDptr[indexFirst], typeIDptr, fxacc, fyacc, fzacc, &virialSumX, &virialSumY,
                               &virialSumZ, &potentialEnergySum, remainder, neighborIndices);
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

    fxptr[indexFirst] += sumfx;
    fyptr[indexFirst] += sumfy;
    fzptr[indexFirst] += sumfz;

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_thread_num();

      // horizontally reduce virialSumX and virialSumY
      const __m256d hSumVirialxy = _mm256_hadd_pd(virialSumX, virialSumY);
      const __m128d hSumVirialxyLow = _mm256_extractf128_pd(hSumVirialxy, 0);
      const __m128d hSumVirialxyHigh = _mm256_extractf128_pd(hSumVirialxy, 1);
      const __m128d hSumVirialxyVec = _mm_add_pd(hSumVirialxyHigh, hSumVirialxyLow);

      // horizontally reduce virialSumZ and potentialEnergySum
      const __m256d hSumVirialzPotentialEnergy = _mm256_hadd_pd(virialSumZ, potentialEnergySum);
      const __m128d hSumVirialzPotentialEnergyLow = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 0);
      const __m128d hSumVirialzPotentialEnergyHigh = _mm256_extractf128_pd(hSumVirialzPotentialEnergy, 1);
      const __m128d hSumVirialzPotentialEnergyVec =
          _mm_add_pd(hSumVirialzPotentialEnergyHigh, hSumVirialzPotentialEnergyLow);

      // globals = {virialX, virialY, virialZ, potentialEnergy}
      double globals[4];
      _mm_store_pd(&globals[0], hSumVirialxyVec);
      _mm_store_pd(&globals[2], hSumVirialzPotentialEnergyVec);

      // we assume newton3 to be enabled in this function call, thus we multiply by two if the value of newton3 is
      // false, since for newton3 disabled we divide by two later on.
      if (newton3) {
        _aosThreadData[threadnum].virialSumN3[0] += globals[0] * 0.5;
        _aosThreadData[threadnum].virialSumN3[1] += globals[1] * 0.5;
        _aosThreadData[threadnum].virialSumN3[2] += globals[2] * 0.5;
        _aosThreadData[threadnum].potentialEnergySumN3 += globals[3] * 0.5;
      } else {
        _aosThreadData[threadnum].virialSumNoN3[0] += globals[0];
        _aosThreadData[threadnum].virialSumNoN3[1] += globals[1];
        _aosThreadData[threadnum].virialSumNoN3[2] += globals[2];
        _aosThreadData[threadnum].potentialEnergySumNoN3 += globals[3];
      }
    }
    // interact with i with 4 neighbors
#endif  // __AVX512F__
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
   * Get the number of flops used per kernel call for a given particle pair. This should count the
   * floating point operations needed for two particles that lie within a cutoff radius, having already calculated the
   * distance.
   * @param molAType molecule A's type id
   * @param molBType molecule B's type id
   * @param newton3 is newton3 applied.
   * @note molAType and molBType make no difference for LJFunctor, but are kept to have a consistent interface for other
   * functors where they may.
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall(size_t molAType, size_t molBType, bool newton3) {
    // Kernel: 12 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply scale) sum
    // Adding to particle forces: 6 or 3 depending newton3
    // Total = 12 + (6 or 3) = 18 or 15
    return newton3 ? 18ul : 15ul;
  }

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
      // We distinguish between non-newton3 and newton3 functor calls. Newton3 calls are accumulated directly.
      // Non-newton3 calls are accumulated temporarily and later divided by 2.
      double potentialEnergySumNoN3Acc = 0;
      std::array<double, 3> virialSumNoN3Acc = {0, 0, 0};
      for (size_t i = 0; i < _aosThreadData.size(); ++i) {
        potentialEnergySumNoN3Acc += _aosThreadData[i].potentialEnergySumNoN3;
        _potentialEnergySum += _aosThreadData[i].potentialEnergySumN3;

        virialSumNoN3Acc += _aosThreadData[i].virialSumNoN3;
        _virialSum += _aosThreadData[i].virialSumN3;
      }
      // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
      // here.
      potentialEnergySumNoN3Acc *= 0.5;
      virialSumNoN3Acc *= 0.5;

      _potentialEnergySum += potentialEnergySumNoN3Acc;
      _virialSum += virialSumNoN3Acc;

      // we have always calculated 6*potentialEnergy, so we divide by 6 here!
      _potentialEnergySum /= 6.;
      _postProcessed = true;

      AutoPasLog(TRACE, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(TRACE, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
    }
  }

  /**
   * Get the potential Energy
   * @return the potential Energy
   */
  double getPotentialEnergy() {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get potential energy even though calculateGlobals is false. If you want this functor to calculate "
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
   * @param sigmaSquared
   */
  void setParticleProperties(double epsilon24, double sigmaSquared) {
#ifdef __AVX512F__
//    _epsilon24 = _mm256_set1_pd(epsilon24);
//    _sigmaSquared = _mm256_set1_pd(sigmaSquared);
//    if constexpr (applyShift) {
//      _shift6 = _mm256_set1_pd(
//          ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquared[0]));
//    } else {
//      _shift6 = _mm256_setzero_pd();
//    }
#endif

    _epsilon24AoS = epsilon24;
    _sigmaSquaredAoS = sigmaSquared;
    if constexpr (applyShift) {
      _shift6AoS = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquaredAoS);
    } else {
      _shift6AoS = 0.;
    }
  }

 private:
  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData()
        : virialSumNoN3{0., 0., 0.},
          virialSumN3{0., 0., 0.},
          potentialEnergySumNoN3{0.},
          potentialEnergySumN3{0.},
          __remainingTo64{} {}
    void setZero() {
      virialSumNoN3 = {0., 0., 0.};
      virialSumN3 = {0., 0., 0.};
      potentialEnergySumNoN3 = 0.;
      potentialEnergySumN3 = 0.;
    }

    // variables
    std::array<double, 3> virialSumNoN3;
    std::array<double, 3> virialSumN3;
    double potentialEnergySumNoN3;
    double potentialEnergySumN3;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 8 * sizeof(double)) / sizeof(double)];
  };
  // make sure of the size of AoSThreadData
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

#ifdef __AVX512F__
  const __m512d _zero{_mm512_set1_pd(0.)};
  const __m512d _one{_mm512_set1_pd(1.)};
  const __m512d _three{_mm512_set1_pd(3.)};
  /**
   * Masks for the remainder cases to avoid loading data beyond the end of the vectors.
   * Masks are generated with decimal numbers, whose binary equivalent consists of 8 0s or 1s, representing which registers
   * are masked when that mask is used.
   */
  std::array<__mmask8,8> _remainderMasks{
      __mmask8(255), // = 11111111
      __mmask8(127), // = 01111111
      __mmask8(63),  // = 00111111
      __mmask8(31),  // = 00011111
      __mmask8(15),  // = 00001111
      __mmask8(7),   // = 00000111
      __mmask8(3),   // = 00000011
      __mmask8(1),   // = 00000001
  };
  const __m512i _ownedStateDummyMM512i{_mm512_set1_epi64(static_cast<int64_t>(autopas::OwnershipState::dummy))};
  const __m512i _ownedStateOwnedMM256i{_mm512_set1_epi64(static_cast<int64_t>(autopas::OwnershipState::owned))};
  const __m512d _cutoffSquared{};
  __m512d _shift6 = _mm512_setzero_pd();
  __m512d _epsilon24{};
  __m512d _sigmaSquared{};
#endif

  const double _cutoffSquaredAoS = 0;
  double _epsilon24AoS, _sigmaSquaredAoS, _shift6AoS = 0;

  ParticlePropertiesLibrary<double, size_t> *_PPLibrary = nullptr;

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadData> _aosThreadData;

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;

  // number of double values that fit into a vector register, i.e. for AVX512 this is 8.
  constexpr static size_t vecLength = 8;
};
}  // namespace mdLib
