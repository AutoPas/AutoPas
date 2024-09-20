// Created by Luis Gall on 27.03.2024
//modified by ivander alson tanjaya 08.09.2024

#pragma once
#include <chrono>
#include "ParticlePropertiesLibrary.h"
#include "VectorizationPatterns.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"
#include "hwy/print-inl.h"

// #undef HWY_TARGET_INCLUDE
// #define HWY_TARGET_INCLUDE "molecularDynamicsLibrary/LJFunctorHWY.h"
// #include <hwy/foreach_target.h>
#include <hwy/highway.h>

// HWY_BEFORE_NAMESPACE();
namespace mdLib {
//namespace HWY_NAMESPACE {

namespace highway = hwy::HWY_NAMESPACE;

// architecture specific information
using VectorDouble = decltype(highway::Zero(tag_double));
using VectorLong = decltype(highway::Zero(tag_long));
using MaskDouble = decltype(highway::FirstN(tag_double, 1));
using MaskLong = decltype(highway::FirstN(tag_long, 2));
//time counter
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool relevantForTuning = true, VectorizationPattern vecPattern = VectorizationPattern::p1xVec>

class LJFunctorSmoothHWYGS
    : public autopas::Functor<Particle,
                              LJFunctorSmoothHWYGS<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  LJFunctorSmoothHWYGS() = delete;

 private:
  explicit LJFunctorSmoothHWYGS(double cutoff, double innerCutoff, void*)
      : autopas::Functor<Particle,
                         LJFunctorSmoothHWYGS<Particle, applyShift, useMixing, useNewton3, calculateGlobals, relevantForTuning>>(cutoff),
        _cutoffSquareAoS{cutoff * cutoff},
        _cutoffAoS{cutoff},
        _innerCutoffAoS{innerCutoff},
        _innerCutoffSquareAoS{innerCutoff*innerCutoff},
        _cutoffDiffAoS{cutoff - innerCutoff},
        _cutoffDiffCubedInvAoS{ 1 / ((cutoff - innerCutoff)*(cutoff - innerCutoff)*(cutoff - innerCutoff))},
        _tripleOuterMinusInnerAoS{3* cutoff - innerCutoff},

        _cutoffSquared{highway::Set(tag_double,cutoff*cutoff)},
        _cutoff{highway::Set(tag_double, cutoff)},
        _innerCutoff{highway::Set(tag_double, innerCutoff)},
        _innerCutoffSquared{highway::Set(tag_double, innerCutoff*innerCutoff)},
        _cutoffDiff{highway::Set(tag_double,cutoff - innerCutoff)},
        _cutoffDiffCubedInv{highway::Set(tag_double, 1 / ((cutoff - innerCutoff)*(cutoff - innerCutoff)*(cutoff - innerCutoff)))},
        _tripleOuterMinusInner{highway::Set(tag_double,3 * cutoff - innerCutoff)},

        _uPotSum{0.},
        _virialSum{0.,0.,0.},
        _aosThreadData{},
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }

    initialize();

    //std::cout << hwy::TargetName(HWY_TARGET) << std::endl; // AutoPasLog(INFO, "Highway Wrapper initialized with a register size of ({}) with architecture {}.", _vecLengthDouble, hwy::TargetName(HWY_TARGET));
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
  explicit LJFunctorSmoothHWYGS(double cutoff, double innerCutoff) : LJFunctorSmoothHWYGS(cutoff,innerCutoff, nullptr) {
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
  explicit LJFunctorSmoothHWYGS(double cutoff, double innerCutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctorSmoothHWYGS(cutoff, innerCutoff,nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final { return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both; }

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

    double smoothing = 1;
    if (dr2 > _cutoffSquareAoS) {
      return;
    }
    else if(dr2 > _innerCutoffSquareAoS){
      double drtrue = std::sqrt(dr2);
      double temp = drtrue -_innerCutoffAoS;
      double _tripleOuterMinusInnerAoSLocal = _tripleOuterMinusInnerAoS;
      double _cutoffDiffCubedInvAoSLocal = _cutoffDiffCubedInvAoS;
      smoothing = 1 - (temp * temp) * (_tripleOuterMinusInnerAoSLocal - 2 * drtrue) * _cutoffDiffCubedInvAoSLocal; //justify change
    }

    double invdr2 = 1. / dr2;
    double lj6 = sigmasquare * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    auto f = dr * fac * smoothing;
    i.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(f);
    }
    if (calculateGlobals) {
      auto virial = dr * f;
      double upot = smoothing * epsilon24 * lj12m6;

      const int threadnum = autopas::autopas_get_thread_num();
      // for non-newton3 the division is in the post-processing step.
      if (newton3) {
        upot *= 0.5;
        virial *= (double)0.5;
      }
      if (i.isOwned()) {
        _aosThreadData[threadnum].uPotSum += upot;
        _aosThreadData[threadnum].virialSum += virial;
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        _aosThreadData[threadnum].uPotSum += upot;
        _aosThreadData[threadnum].virialSum += virial;
      }
    }
  }

  /**
                 * @copydoc Functor::SoAFunctorSingle(SoAView<SoAArraysType> soa, bool newton3)
                 * This functor will always do a newton3 like traversal of the soa.
                 * However, it still needs to know about newton3 to correctly add up the global values.
   */
  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (newton3) {
      SoAFunctorSingleImplGS<true>(soa);
    } else {
      SoAFunctorSingleImplGS<false>(soa);
    }
  }

  // clang-format off
                /**
                 * @copydoc Functor::SoAFunctorPair(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3)
                 */
  // clang-format on
  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImplGS<true>(soa1, soa2);
    } else {
      SoAFunctorPairImplGS<false>(soa1, soa2);
    }
  }

 private:

  inline void initialize() {

    for (size_t n = 0; n <_vecLengthDouble-1;++n) {
      restMasksDouble[n] = highway::FirstN(tag_double, n+1);
      restMasksLong[n] = highway::FirstN(tag_long, n+1);
      _remainderMask[n] =  highway::Not(highway::FirstN(tag_double,n));
      asc[n] = n;
    }
    asc[_vecLengthDouble-1] = _vecLengthDouble-1;
    _ascendingIndices = highway::LoadU(tag_double,asc);
    _remainderMask[_vecLengthDouble-1] =  highway::Not(highway::FirstN(tag_double,_vecLengthDouble-1));
  }

  template <bool reversed>
  inline bool checkFirstLoopCondition(const size_t i, const size_t stop) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      return reversed ? (long)i >= 0 : (i < stop);
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      return reversed ? (long)i >= 1 : (i < stop - 1);
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      // TODO : Implement
      return false;
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // TODO : Implement
      return false;
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
      // TODO : Implement
      return false;
    }
    else {
      return false;
    }
  }

  inline void decrementFirstLoop(size_t &i) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      --i;
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      i -= 2;
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
      // TODO : Implement
    }
    else {

    }
  }

  inline void incrementFirstLoop(unsigned int &i) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      ++i;
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      i += 2;
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
      // TODO : Implement
    }
    else {

    }
  }

  template <bool reversed>
  inline int obtainFirstLoopRest(const int i, const int stop) {
    return reversed ? (i < 0 ? 0 : i+1) : stop-i;
  }

  inline bool checkSecondLoopCondition(const size_t i, const size_t j) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      return j < (i & ~(_vecLengthDouble - 1));
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      return j < (i & ~(_vecLengthDouble/2 - 1));
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      // TODO : Implement
      return false;
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // TODO : Implement
      return false;
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
      // TODO : Implement
      return false;
    }
    else {
      return false;
    }
  }

  inline void incrementSecondLoop(unsigned int &j) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      j += _vecLengthDouble;
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      j += _vecLengthDouble/2;
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
      // TODO : Implement
    }
    else {

    }
  }

  inline int obtainSecondLoopRest(const size_t i) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      return (int)(i & (_vecLengthDouble - 1));
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      return (int)(i & (_vecLengthDouble/2 - 1));
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      // TODO : implement
      return -1;
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // TODO : implement
      return -1;
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
      // TODO : implement
      return -1;
    }
    else {
      return -1;
    }
  }

  template <bool remainder, bool reversed>
  inline void fillIRegisters(const size_t i, const double *const __restrict xPtr, const double *const __restrict yPtr,
                             const double *const __restrict zPtr, const autopas::OwnershipState *const __restrict ownedStatePtr,
                             VectorDouble& x1, VectorDouble& y1, VectorDouble& z1, MaskDouble& ownedMaskI, const size_t restI) {

    VectorDouble ownedStateIDouble = _zeroDouble;

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      int64_t owned = static_cast<int64_t>(ownedStatePtr[i]);
      ownedStateIDouble = highway::Set(tag_double, static_cast<double>(owned));

      x1 = highway::Set(tag_double, xPtr[i]);
      y1 = highway::Set(tag_double, yPtr[i]);
      z1 = highway::Set(tag_double, zPtr[i]);
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      int64_t ownedFirst = static_cast<int64_t>(ownedStatePtr[i]);
      ownedStateIDouble = highway::Set(tag_double, static_cast<double>(ownedFirst));

      x1 = highway::Set(tag_double, xPtr[i]);
      y1 = highway::Set(tag_double, yPtr[i]);
      z1 = highway::Set(tag_double, zPtr[i]);

      VectorDouble tmpOwnedI = _zeroDouble;
      VectorDouble tmpX1 = _zeroDouble;
      VectorDouble tmpY1 = _zeroDouble;
      VectorDouble tmpZ1 = _zeroDouble;

      if constexpr (!remainder) {
        int index = reversed ? i-1 : i+1;
        int64_t ownedSecond = static_cast<int64_t>(ownedStatePtr[index]);
        tmpOwnedI = highway::Set(tag_double, static_cast<double>(ownedSecond));
        tmpX1 = highway::Set(tag_double, xPtr[index]);
        tmpY1 = highway::Set(tag_double, yPtr[index]);
        tmpZ1 = highway::Set(tag_double, zPtr[index]);
      }

      ownedStateIDouble = highway::ConcatLowerLower(tag_double, tmpOwnedI, ownedStateIDouble);
      x1 = highway::ConcatLowerLower(tag_double, tmpX1, x1);
      y1 = highway::ConcatLowerLower(tag_double, tmpY1, y1);
      z1 = highway::ConcatLowerLower(tag_double, tmpZ1, z1);
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      // TODO : implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // TODO : implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
      // TODO : implement
    }

    ownedMaskI = highway::Ne(ownedStateIDouble, _zeroDouble);
  }

  template <bool remainder>
  inline void handleNewton3Reduction(const VectorDouble& fx, const VectorDouble& fy, const VectorDouble& fz,
                                     double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                                     double *const __restrict fz2Ptr, const size_t j, const size_t rest) {

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      const VectorDouble fx2 = remainder ? highway::MaskedLoad(restMasksDouble[rest-1], tag_double, &fx2Ptr[j])
                                         : highway::LoadU(tag_double, &fx2Ptr[j]);
      const VectorDouble fy2 = remainder ? highway::MaskedLoad(restMasksDouble[rest-1], tag_double, &fy2Ptr[j])
                                         : highway::LoadU(tag_double, &fy2Ptr[j]);
      const VectorDouble fz2 = remainder ? highway::MaskedLoad(restMasksDouble[rest-1], tag_double, &fz2Ptr[j])
                                         : highway::LoadU(tag_double, &fz2Ptr[j]);

      const VectorDouble fx2New = fx2 - fx;
      const VectorDouble fy2New = fy2 - fy;
      const VectorDouble fz2New = fz2 - fz;

      remainder ? highway::BlendedStore(fx2New, restMasksDouble[rest-1], tag_double, &fx2Ptr[j])
                : highway::StoreU(fx2New, tag_double, &fx2Ptr[j]);
      remainder ? highway::BlendedStore(fy2New, restMasksDouble[rest-1], tag_double, &fy2Ptr[j])
                : highway::StoreU(fy2New, tag_double, &fy2Ptr[j]);
      remainder ? highway::BlendedStore(fz2New, restMasksDouble[rest-1], tag_double, &fz2Ptr[j])
                : highway::StoreU(fz2New, tag_double, &fz2Ptr[j]);
    }

    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
        const auto lowerFx = highway::LowerHalf(fx);
        const auto lowerFy = highway::LowerHalf(fy);
        const auto lowerFz = highway::LowerHalf(fz);

        const auto upperFx = highway::UpperHalf(tag_half_double, fx);
        const auto upperFy = highway::UpperHalf(tag_half_double, fy);
        const auto upperFz = highway::UpperHalf(tag_half_double, fz);

        const auto fxCombined = lowerFx + upperFx;
        const auto fyCombined = lowerFy + upperFy;
        const auto fzCombined = lowerFz + upperFz;

        const auto fxCombinedExt = highway::ZeroExtendVector(tag_double, fxCombined);
        const auto fyCombinedExt = highway::ZeroExtendVector(tag_double, fyCombined);
        const auto fzCombinedExt = highway::ZeroExtendVector(tag_double, fzCombined);

        const auto mask = restMasksDouble[remainder ? rest-1 : _vecLengthDouble/2-1];
        VectorDouble fx2 = highway::MaskedLoad(mask, tag_double, &fx2Ptr[j]);
        VectorDouble fy2 = highway::MaskedLoad(mask, tag_double, &fy2Ptr[j]);
        VectorDouble fz2 = highway::MaskedLoad(mask, tag_double, &fz2Ptr[j]);

        auto newFx = fx2 - fxCombinedExt;
        auto newFy = fy2 - fyCombinedExt;
        auto newFz = fz2 - fzCombinedExt;

        highway::BlendedStore(newFx, mask, tag_double, &fx2Ptr[j]);
        highway::BlendedStore(newFy, mask, tag_double, &fy2Ptr[j]);
        highway::BlendedStore(newFz, mask, tag_double, &fz2Ptr[j]);
    }


    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
      // TODO : Implement
    }
  }


  template <bool remainder>
  inline void handleNewton3ReductionGS(VectorDouble indices, const VectorDouble& fx, const VectorDouble& fy, const VectorDouble& fz,
                                       double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                                       double *const __restrict fz2Ptr, const size_t j, const size_t rest) {

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      const VectorDouble fx2 = remainder ? highway::MaskedGatherIndex(restMasksDouble[rest-1],tag_double,fx2Ptr,highway::ConvertTo(tag_long,indices))
                                         : highway::GatherIndex(tag_double,fx2Ptr,highway::ConvertTo(tag_long,indices));
      const VectorDouble fy2 = remainder ? highway::MaskedGatherIndex(restMasksDouble[rest-1],tag_double,fy2Ptr,highway::ConvertTo(tag_long,indices))
                                         : highway::GatherIndex(tag_double,fy2Ptr,highway::ConvertTo(tag_long,indices));
      const VectorDouble fz2 = remainder ? highway::MaskedGatherIndex(restMasksDouble[rest-1],tag_double,fz2Ptr,highway::ConvertTo(tag_long,indices))
                                         : highway::GatherIndex(tag_double,fz2Ptr,highway::ConvertTo(tag_long,indices));

      const VectorDouble fx2New = fx2 - fx;
      const VectorDouble fy2New = fy2 - fy;
      const VectorDouble fz2New = fz2 - fz;

      remainder ? highway::MaskedScatterIndex(fx2New,restMasksDouble[rest-1],tag_double,fx2Ptr,highway::ConvertTo(tag_long,indices))
                : highway::ScatterIndex(fx2New,tag_double,fx2Ptr,highway::ConvertTo(tag_long,indices));
      remainder ? highway::MaskedScatterIndex(fy2New,restMasksDouble[rest-1],tag_double,fy2Ptr,highway::ConvertTo(tag_long,indices))
                : highway::ScatterIndex(fy2New,tag_double,fy2Ptr,highway::ConvertTo(tag_long,indices));
      remainder ? highway::MaskedScatterIndex(fz2New,restMasksDouble[rest-1],tag_double,fz2Ptr,highway::ConvertTo(tag_long,indices))
                : highway::ScatterIndex(fz2New,tag_double,fz2Ptr,highway::ConvertTo(tag_long,indices));
    }

    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
        const auto lowerFx = highway::LowerHalf(fx);
        const auto lowerFy = highway::LowerHalf(fy);
        const auto lowerFz = highway::LowerHalf(fz);

        const auto upperFx = highway::UpperHalf(tag_half_double, fx);
        const auto upperFy = highway::UpperHalf(tag_half_double, fy);
        const auto upperFz = highway::UpperHalf(tag_half_double, fz);

        const auto fxCombined = lowerFx + upperFx;
        const auto fyCombined = lowerFy + upperFy;
        const auto fzCombined = lowerFz + upperFz;

        const auto fxCombinedExt = highway::ZeroExtendVector(tag_double, fxCombined);
        const auto fyCombinedExt = highway::ZeroExtendVector(tag_double, fyCombined);
        const auto fzCombinedExt = highway::ZeroExtendVector(tag_double, fzCombined);

        const auto mask = restMasksDouble[remainder ? rest-1 : _vecLengthDouble/2-1];
        VectorDouble fx2 = highway::IfThenElseZero(mask, highway::GatherIndex(tag_double,fx2Ptr,highway::ConvertTo(tag_long,indices)));
        VectorDouble fy2 = highway::IfThenElseZero(mask, highway::GatherIndex(tag_double,fy2Ptr,highway::ConvertTo(tag_long,indices)));
        VectorDouble fz2 = highway::IfThenElseZero(mask, highway::GatherIndex(tag_double,fz2Ptr,highway::ConvertTo(tag_long,indices)));

        auto newFx = fx2 - fxCombinedExt;
        auto newFy = fy2 - fyCombinedExt;
        auto newFz = fz2 - fzCombinedExt;

        highway::ScatterIndex(highway::IfThenElseZero(mask, newFx),tag_double,fx2Ptr,highway::ConvertTo(tag_long,indices));
        highway::ScatterIndex(highway::IfThenElseZero(mask, newFy),tag_double,fx2Ptr,highway::ConvertTo(tag_long,indices));
        highway::ScatterIndex(highway::IfThenElseZero(mask, newFz),tag_double,fx2Ptr,highway::ConvertTo(tag_long,indices));
    }

    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
      // TODO : Implement
    }
  }

  template <bool reversed, bool remainder>
  inline void reduceAccumulatedForce(const size_t i, double *const __restrict fxPtr, double *const __restrict fyPtr,
                                     double *const __restrict fzPtr, const VectorDouble& fxAcc, const VectorDouble& fyAcc, const VectorDouble& fzAcc, const int restI) {

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
      fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
      fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
    }

    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
        const auto lowerFxAcc = highway::LowerHalf(fxAcc);
        const auto lowerFyAcc = highway::LowerHalf(fyAcc);
        const auto lowerFzAcc = highway::LowerHalf(fzAcc);

        const auto lowerFxAccExt = highway::ZeroExtendVector(tag_double, lowerFxAcc);
        const auto lowerFyAccExt = highway::ZeroExtendVector(tag_double, lowerFyAcc);
        const auto lowerFzAccExt = highway::ZeroExtendVector(tag_double, lowerFzAcc);

        fxPtr[i] += highway::ReduceSum(tag_double, lowerFxAccExt);
        fyPtr[i] += highway::ReduceSum(tag_double, lowerFyAccExt);
        fzPtr[i] += highway::ReduceSum(tag_double, lowerFzAccExt);

        if constexpr (!remainder) {
            const auto upperFxAcc = highway::UpperHalf(tag_half_double, fxAcc);
            const auto upperFyAcc = highway::UpperHalf(tag_half_double, fyAcc);
            const auto upperFzAcc = highway::UpperHalf(tag_half_double, fzAcc);

            const auto upperFxAccExt = highway::ZeroExtendVector(tag_double, upperFxAcc);
            const auto upperFyAccExt = highway::ZeroExtendVector(tag_double, upperFyAcc);
            const auto upperFzAccExt = highway::ZeroExtendVector(tag_double, upperFzAcc);

            int index = reversed ? i-1 : i+1;
            fxPtr[index] += highway::ReduceSum(tag_double, upperFxAccExt);
            fyPtr[index] += highway::ReduceSum(tag_double, upperFyAccExt);
            fzPtr[index] += highway::ReduceSum(tag_double, upperFzAccExt);
        }
    }

    else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // TODO : Implement
    }
    else if constexpr (vecPattern == VectorizationPattern::pVecxVec) {
      // TODO : Implement
    }
  }

  template <bool newton3>
  inline void computeGlobals(const VectorDouble& virialSumX, const VectorDouble& virialSumY,
                             const VectorDouble& virialSumZ, const VectorDouble& uPotSum) {
    const int threadnum = autopas::autopas_get_thread_num();

    double globals[4] {
        highway::ReduceSum(tag_double, virialSumX),
        highway::ReduceSum(tag_double, virialSumY),
        highway::ReduceSum(tag_double, virialSumZ),
        highway::ReduceSum(tag_double, uPotSum)
    };

    double factor = 1.;
    if constexpr (newton3) {
      factor = 0.5;
    }

    _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
    _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
    _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
    _aosThreadData[threadnum].uPotSum += globals[3] * factor;
  }

  template <bool reversed, bool newton3, bool remainderI>
  inline void handleILoopBody(const size_t i, const double *const __restrict xPtr1, const double *const __restrict yPtr1,
                              const double *const __restrict zPtr1, const autopas::OwnershipState *const __restrict ownedStatePtr1,
                              const double *const __restrict xPtr2, const double *const __restrict yPtr2,
                              const double *const __restrict zPtr2, const autopas::OwnershipState *const __restrict ownedStatePtr2,
                              double *const __restrict fxPtr1, double *const __restrict fyPtr1, double *const __restrict fzPtr1,
                              double *const __restrict fxPtr2, double *const __restrict fyPtr2, double *const __restrict fzPtr2,
                              const size_t * const __restrict typeIDptr1, const size_t * const __restrict typeIDptr2,
                              VectorDouble& virialSumX, VectorDouble& virialSumY, VectorDouble& virialSumZ, VectorDouble& uPotSum,
                              const size_t restI, const size_t jStop) {

    VectorDouble fxAcc = _zeroDouble;
    VectorDouble fyAcc = _zeroDouble;
    VectorDouble fzAcc = _zeroDouble;

    MaskDouble ownedMaskI;

    VectorDouble x1 = _zeroDouble;
    VectorDouble y1 = _zeroDouble;
    VectorDouble z1 = _zeroDouble;

    fillIRegisters<remainderI, reversed>(i, xPtr1, yPtr1, zPtr1, ownedStatePtr1, x1, y1, z1, ownedMaskI, restI);

    unsigned int j = 0;
    for (; checkSecondLoopCondition(jStop, j); incrementSecondLoop(j)) {

      SoAKernel<newton3, remainderI, false>(j, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                            x1, y1, z1, xPtr2, yPtr2, zPtr2, fxPtr2, fyPtr2, fzPtr2, &typeIDptr1[i], typeIDptr2,
                                            fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, restI, 0);
    }

    const int restJ = obtainSecondLoopRest(jStop);
    if (restJ > 0) {

      SoAKernel<newton3, remainderI, true>(j, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                           x1, y1, z1, xPtr2, yPtr2, zPtr2, fxPtr2, fyPtr2, fzPtr2, &typeIDptr1[i], typeIDptr2,
                                           fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, restI, restJ);
    }

    reduceAccumulatedForce<reversed, remainderI>(i, fxPtr1, fyPtr1, fzPtr1, fxAcc, fyAcc, fzAcc, restI);
  }

  inline VectorDouble alignr(VectorDouble a,  VectorDouble b, int shift){
    alignas(64) double concatenated[2 * _vecLengthDouble];

    // Store vectors 'a' and 'b' into concatenated array
    Store(a, tag_double, concatenated);
    Store(b, tag_double, concatenated + _vecLengthDouble);

    // Create a vector to hold the result
    VectorDouble result;

    // Load the result starting from the appropriate position
    result = highway::LoadU(tag_double,concatenated+shift);
    return result;
  }


  template <bool reversed, bool newton3, bool remainderI>
  inline void handleILoopBodyGS(const size_t i, const double *const __restrict xPtr1, const double *const __restrict yPtr1,
                                const double *const __restrict zPtr1, const autopas::OwnershipState *const __restrict ownedStatePtr1,
                                const double *const __restrict xPtr2, const double *const __restrict yPtr2,
                                const double *const __restrict zPtr2, const autopas::OwnershipState *const __restrict ownedStatePtr2,
                                double *const __restrict fxPtr1, double *const __restrict fyPtr1, double *const __restrict fzPtr1,
                                double *const __restrict fxPtr2, double *const __restrict fyPtr2, double *const __restrict fzPtr2,
                                const size_t * const __restrict typeIDptr1, const size_t * const __restrict typeIDptr2,
                                VectorDouble& virialSumX, VectorDouble& virialSumY, VectorDouble& virialSumZ, VectorDouble& uPotSum,
                                const size_t restI, const size_t jStop) {

    VectorDouble fxAcc = _zeroDouble;
    VectorDouble fyAcc = _zeroDouble;
    VectorDouble fzAcc = _zeroDouble;

    MaskDouble ownedMaskI;

    VectorDouble x1 = _zeroDouble;
    VectorDouble y1 = _zeroDouble;
    VectorDouble z1 = _zeroDouble;

    fillIRegisters<remainderI, reversed>(i, xPtr1, yPtr1, zPtr1, ownedStatePtr1, x1, y1, z1, ownedMaskI, restI);

    VectorDouble interactionIndices =_zeroDouble;
    VectorDouble interactionInnerIndices =_zeroDouble;
    long numAssignedRegisters = 0; // = number of non-empty registers in interactionIndices
    long numAssignedInnerRegisters = 0; // = number of non-empty registers in interactionIndices that are inside the inner cutoff

    unsigned int j = 0;
    for (; checkSecondLoopCondition(jStop, j); incrementSecondLoop(j)) {
      VectorDouble loopIndices = highway::Add(highway::Set(tag_double,j),_ascendingIndices);
      VectorDouble x2;
      VectorDouble y2;
      VectorDouble z2;
      VectorDouble ownedStateJDouble;

      fillJRegisters<false>(j, xPtr2, yPtr2, zPtr2, reinterpret_cast<const int64_t *>(ownedStatePtr2), x2, y2, z2, ownedStateJDouble, 0);

      // distance calculations
      const auto drX = x1 - x2;
      const auto drY = y1 - y2;
      const auto drZ = z1 - z2;

      const auto drX2 = drX * drX;
      const auto drY2 = drY * drY;
      const auto drZ2 = drZ * drZ;

      const auto dr2 = drX2 + drY2 + drZ2;

      const auto cutoffMask = highway::Le(dr2, _cutoffSquared);
      const auto dummyMask = highway::And(ownedMaskI, highway::Ne(ownedStateJDouble, _ownedStateDummy));
      const auto zeroMask = highway::Ne(dr2, _zeroDouble);
      auto cutoffDummyMask = highway::And(zeroMask, highway::And(cutoffMask, dummyMask));
      //huge difference!!!! 4-6 times speedup
      if (highway::AllFalse(tag_double, cutoffDummyMask)) {
            return;
      }
      const auto innerCutoffMask = highway::Ge(dr2,_innerCutoffSquared);
      auto innerCutoffDummyMask = highway::And(zeroMask, highway::And(innerCutoffMask, dummyMask));

      const auto temp = highway::And(innerCutoffDummyMask, cutoffDummyMask);
      const auto temp2 = highway::Xor(innerCutoffDummyMask, cutoffDummyMask);
      cutoffDummyMask = highway::And(cutoffDummyMask,temp2);
      innerCutoffDummyMask = temp;

      const long popCountMask = highway::CountTrue(tag_double, cutoffDummyMask);
      const long popCountInnerMask = highway::CountTrue(tag_double,innerCutoffDummyMask);

      //gathering indices for the normal functor
      if (numAssignedRegisters + popCountMask < _vecLengthDouble) {
        VectorDouble newInteractionIndices = highway::Compress(loopIndices,cutoffDummyMask);
        newInteractionIndices = highway::IfThenElseZero(highway::FirstN(tag_double,popCountMask),newInteractionIndices);
        #if HWY_TARGET <= HWY_AVX3
        interactionIndices = alignr(interactionIndices,newInteractionIndices,popCountMask); // there is a native intrinsic for this in avx512 but im not sure how to call it using highway
        #else
        interactionIndices = alignr(interactionIndices,newInteractionIndices,popCountMask);
        #endif
        numAssignedRegisters += popCountMask;
      }
      else {
        VectorDouble newInteractionIndices = highway::Compress(loopIndices,cutoffDummyMask);
        newInteractionIndices = highway::IfThenElseZero(highway::FirstN(tag_double,popCountMask),newInteractionIndices);
        interactionIndices = alignr(interactionIndices,newInteractionIndices,_vecLengthDouble-numAssignedRegisters);

        SoAKernelGS<newton3, remainderI, false>(interactionIndices,
                                                j, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2,
                                                zPtr2, fxPtr2, fyPtr2, fzPtr2, &typeIDptr1[i], typeIDptr2, fxAcc, fyAcc, fzAcc, virialSumX,
                                                virialSumY, virialSumZ, uPotSum, restI, 0);

        const auto _alreadyProcessedMask =  _remainderMask[_vecLengthDouble-numAssignedRegisters];
        interactionIndices = highway::IfThenElseZero(_alreadyProcessedMask,newInteractionIndices);
        interactionIndices = alignr(_zeroDouble,interactionIndices,popCountMask);
        numAssignedRegisters = popCountMask - (_vecLengthDouble - numAssignedRegisters);
      }

      //gathering indices for the smoothed functor
      if (numAssignedInnerRegisters + popCountInnerMask < _vecLengthDouble) {
        VectorDouble newInteractionInnerIndices = highway::Compress(loopIndices, innerCutoffDummyMask);
        newInteractionInnerIndices = highway::IfThenElseZero(highway::FirstN(tag_double,popCountInnerMask),newInteractionInnerIndices);
        interactionInnerIndices = alignr(interactionInnerIndices,newInteractionInnerIndices,popCountInnerMask);
        numAssignedInnerRegisters += popCountInnerMask;
      }
      else {
        VectorDouble newInteractionInnerIndices = highway::Compress(loopIndices, innerCutoffDummyMask);
        newInteractionInnerIndices = highway::IfThenElseZero(highway::FirstN(tag_double,popCountInnerMask),newInteractionInnerIndices);
        interactionInnerIndices = alignr(interactionInnerIndices,newInteractionInnerIndices,_vecLengthDouble-numAssignedInnerRegisters);

        SoAKernelSmoothGS<newton3, remainderI, false>(interactionInnerIndices,
                                                      j, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2,
                                                      zPtr2, fxPtr2, fyPtr2, fzPtr2, &typeIDptr1[i], typeIDptr2, fxAcc, fyAcc, fzAcc, virialSumX,
                                                      virialSumY, virialSumZ, uPotSum, restI, 0);

        const auto _alreadyProcessedMask =  _remainderMask[_vecLengthDouble-numAssignedInnerRegisters];
        interactionInnerIndices = highway::IfThenElseZero(_alreadyProcessedMask,newInteractionInnerIndices);
        interactionInnerIndices = alignr(_zeroDouble,interactionInnerIndices,popCountInnerMask);
        numAssignedInnerRegisters = popCountInnerMask - (_vecLengthDouble - numAssignedInnerRegisters);
      }

    }
    //process the remaining indices for both
    if(numAssignedRegisters > 0){
      interactionIndices = alignr(interactionIndices,_zeroDouble,_vecLengthDouble-numAssignedRegisters);
      SoAKernelGS<newton3, remainderI, true>(interactionIndices,j, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                             x1, y1, z1, xPtr2, yPtr2, zPtr2, fxPtr2, fyPtr2, fzPtr2, &typeIDptr1[i], typeIDptr2,
                                             fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, restI, numAssignedRegisters);
    }
    if(numAssignedInnerRegisters > 0){
      interactionInnerIndices = alignr(interactionInnerIndices,_zeroDouble,_vecLengthDouble-numAssignedInnerRegisters);
      SoAKernelSmoothGS<newton3, remainderI, true>(interactionInnerIndices,j, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                                   x1, y1, z1, xPtr2, yPtr2, zPtr2, fxPtr2, fyPtr2, fzPtr2, &typeIDptr1[i], typeIDptr2,
                                                   fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, restI, numAssignedInnerRegisters);
    }
    const int restJ = obtainSecondLoopRest(jStop);
    if (restJ > 0) {
      SoAKernelSmooth<newton3, remainderI, true>(j, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2),
                                                 x1, y1, z1, xPtr2, yPtr2, zPtr2, fxPtr2, fyPtr2, fzPtr2, &typeIDptr1[i], typeIDptr2,
                                                 fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, restI, restJ);
    }
    reduceAccumulatedForce<reversed, remainderI>(i, fxPtr1, fyPtr1, fzPtr1, fxAcc, fyAcc, fzAcc, restI);
  }

  template <bool newton3>
  inline void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {

    if (soa.size() == 0) return;

    // obtain iterators for the various values
    const auto *const __restrict xPtr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxPtr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

    // initialize and declare vector variables
    auto virialSumX = _zeroDouble;
    auto virialSumY = _zeroDouble;
    auto virialSumZ = _zeroDouble;
    auto uPotSum = _zeroDouble;

    size_t i = soa.size() - 1;
    for (; checkFirstLoopCondition<true>(i, 0); decrementFirstLoop(i)) {
      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");

      handleILoopBody<true, true, false>(i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr, ownedStatePtr,
                                         fxPtr, fyPtr, fzPtr, fxPtr, fyPtr, fzPtr, typeIDptr, typeIDptr, virialSumX, virialSumY, virialSumZ, uPotSum, 0, i);
    }

    if constexpr (vecPattern != VectorizationPattern::p1xVec) {
      // Rest I can't occur in 1xVec case
      const int restI = obtainFirstLoopRest<true>(i, -1);

      if (restI > 0) {
        handleILoopBody<true, true, true>(i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr, ownedStatePtr,
                                          fxPtr, fyPtr, fzPtr, fxPtr, fyPtr, fzPtr, typeIDptr, typeIDptr, virialSumX, virialSumY, virialSumZ, uPotSum, restI, i);
      }
    }

    if constexpr (calculateGlobals) {
      computeGlobals<newton3>(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  template <bool newton3>
  inline void SoAFunctorSingleImplGS(autopas::SoAView<SoAArraysType> soa) {

    if (soa.size() == 0) return;

    // obtain iterators for the various values
    const auto *const __restrict xPtr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fxPtr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle::AttributeNames::typeId>();

    // initialize and declare vector variables
    auto virialSumX = _zeroDouble;
    auto virialSumY = _zeroDouble;
    auto virialSumZ = _zeroDouble;
    auto uPotSum = _zeroDouble;

    size_t i = soa.size() - 1;
    for (; checkFirstLoopCondition<true>(i, 0); decrementFirstLoop(i)) {

      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");

      handleILoopBodyGS<true, true, false>(i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr, ownedStatePtr,
                                           fxPtr, fyPtr, fzPtr, fxPtr, fyPtr, fzPtr, typeIDptr, typeIDptr, virialSumX, virialSumY, virialSumZ, uPotSum, 0, i);
    }

    if constexpr (vecPattern != VectorizationPattern::p1xVec) {
      // Rest I can't occur in 1xVec case
      const int restI = obtainFirstLoopRest<true>(i, -1);

      if (restI > 0) {
        handleILoopBodyGS<true, true, true>(i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr, ownedStatePtr,
                                            fxPtr, fyPtr, fzPtr, fxPtr, fyPtr, fzPtr, typeIDptr, typeIDptr, virialSumX, virialSumY, virialSumZ, uPotSum, restI, i);
      }
    }

    if constexpr (calculateGlobals) {
      computeGlobals<newton3>(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }



  template <bool newton3>
  inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {

    if (soa1.size() == 0 || soa2.size() == 0) {
      return;
    }

    const auto *const __restrict x1Ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1Ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1Ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2Ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2Ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2Ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1Ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1Ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1Ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2Ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2Ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2Ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeID1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    const auto *const __restrict typeID2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

    VectorDouble virialSumX = _zeroDouble;
    VectorDouble virialSumY = _zeroDouble;
    VectorDouble virialSumZ = _zeroDouble;
    VectorDouble uPotSum = _zeroDouble;

    unsigned int i = 0;
    for (; checkFirstLoopCondition<false>(i, soa1.size()); incrementFirstLoop(i)) {
      handleILoopBody<false, newton3, false>(i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2,
                                             fx1Ptr, fy1Ptr, fz1Ptr, fx2Ptr, fy2Ptr, fz2Ptr, typeID1ptr, typeID2ptr, virialSumX, virialSumY, virialSumZ, uPotSum, 0, soa2.size());
    }

    if constexpr (vecPattern != VectorizationPattern::p1xVec) {
      // Rest I can't occur in 1xVec case
      const int restI = obtainFirstLoopRest<false>(i, soa1.size());
      if (restI > 0) {
        handleILoopBody<false, newton3, true>(i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2,
                                              fx1Ptr, fy1Ptr, fz1Ptr, fx2Ptr, fy2Ptr, fz2Ptr, typeID1ptr, typeID2ptr, virialSumX, virialSumY, virialSumZ, uPotSum, restI, soa2.size());
      }
    }

    if constexpr (calculateGlobals) {
      computeGlobals<newton3>(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  template <bool newton3>
  inline void SoAFunctorPairImplGS(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {

    if (soa1.size() == 0 || soa2.size() == 0) {
      return;
    }

    const auto *const __restrict x1Ptr = soa1.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y1Ptr = soa1.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z1Ptr = soa1.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict x2Ptr = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict y2Ptr = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict z2Ptr = soa2.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

    auto *const __restrict fx1Ptr = soa1.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy1Ptr = soa1.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz1Ptr = soa1.template begin<Particle::AttributeNames::forceZ>();
    auto *const __restrict fx2Ptr = soa2.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fy2Ptr = soa2.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fz2Ptr = soa2.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeID1ptr = soa1.template begin<Particle::AttributeNames::typeId>();
    const auto *const __restrict typeID2ptr = soa2.template begin<Particle::AttributeNames::typeId>();

    VectorDouble virialSumX = _zeroDouble;
    VectorDouble virialSumY = _zeroDouble;
    VectorDouble virialSumZ = _zeroDouble;
    VectorDouble uPotSum = _zeroDouble;

    unsigned int i = 0;
    for (; checkFirstLoopCondition<false>(i, soa1.size()); incrementFirstLoop(i)) {
      handleILoopBodyGS<false, newton3, false>(i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2,
                                               fx1Ptr, fy1Ptr, fz1Ptr, fx2Ptr, fy2Ptr, fz2Ptr, typeID1ptr, typeID2ptr, virialSumX, virialSumY, virialSumZ, uPotSum, 0, soa2.size());
    }

    if constexpr (vecPattern != VectorizationPattern::p1xVec) {
      // Rest I can't occur in 1xVec case
      const int restI = obtainFirstLoopRest<false>(i, soa1.size());
      if (restI > 0) {
        handleILoopBodyGS<false, newton3, true>(i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2,
                                                fx1Ptr, fy1Ptr, fz1Ptr, fx2Ptr, fy2Ptr, fz2Ptr, typeID1ptr, typeID2ptr, virialSumX, virialSumY, virialSumZ, uPotSum, restI, soa2.size());
      }
    }

    if constexpr (calculateGlobals) {
      computeGlobals<newton3>(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  template <bool remainder>
  inline void fillJRegisters(const size_t j, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
                             const double *const __restrict z2Ptr, const int64_t *const __restrict ownedStatePtr2,
                             VectorDouble& x2, VectorDouble& y2, VectorDouble& z2, VectorDouble& ownedStateJDouble, const unsigned int rest) {

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {

      if constexpr (remainder) {
        const auto restMaskDouble = restMasksDouble[rest-1];
        const auto restMaskLong = restMasksLong[rest-1];

        x2 = highway::MaskedLoad(restMaskDouble, tag_double, &x2Ptr[j]);
        y2 = highway::MaskedLoad(restMaskDouble, tag_double, &y2Ptr[j]);
        z2 = highway::MaskedLoad(restMaskDouble, tag_double, &z2Ptr[j]);

        const VectorLong ownedStateJ = highway::MaskedLoad(restMaskLong, tag_long, &ownedStatePtr2[j]);
        ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
      }
      else {
        x2 = highway::LoadU(tag_double, &x2Ptr[j]);
        y2 = highway::LoadU(tag_double, &y2Ptr[j]);
        z2 = highway::LoadU(tag_double, &z2Ptr[j]);

        const VectorLong ownedStateJ = highway::LoadU(tag_long, &ownedStatePtr2[j]);
        ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
      }
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {

      const auto restMaskDouble = restMasksDouble[remainder ? rest-1 : _vecLengthDouble/2-1];
      const auto restMaskLong = restMasksLong[remainder ? rest-1 : _vecLengthDouble/2-1];

      VectorLong ownedStateJ = highway::MaskedLoad(restMaskLong, tag_long, &ownedStatePtr2[j]);
      x2 = highway::MaskedLoad(restMaskDouble, tag_double, &x2Ptr[j]);
      y2 = highway::MaskedLoad(restMaskDouble, tag_double, &y2Ptr[j]);
      z2 = highway::MaskedLoad(restMaskDouble, tag_double, &z2Ptr[j]);

      // "broadcast" lower half to upper half
      ownedStateJ = highway::ConcatLowerLower(tag_long, ownedStateJ, ownedStateJ);
      ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
      x2 = highway::ConcatLowerLower(tag_double, x2, x2);
      y2 = highway::ConcatLowerLower(tag_double, y2, y2);
      z2 = highway::ConcatLowerLower(tag_double, z2, z2);
    }
  }

  template <bool remainder>
  inline void fillJRegistersGS(const VectorDouble indices, const size_t j, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
                               const double *const __restrict z2Ptr, const int64_t *const __restrict ownedStatePtr2,
                               VectorDouble& x2, VectorDouble& y2, VectorDouble& z2, VectorDouble& ownedStateJDouble, const unsigned int rest) {

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {

      if constexpr (remainder) {
        const auto restMaskDouble = restMasksDouble[rest-1];
        const auto restMaskLong = restMasksLong[rest-1];
        //check this
        x2 = highway::MaskedGatherIndex(restMaskDouble,tag_double,x2Ptr,highway::ConvertTo(tag_long,indices));
        y2 = highway::MaskedGatherIndex(restMaskDouble,tag_double,y2Ptr,highway::ConvertTo(tag_long,indices));
        z2 = highway::MaskedGatherIndex(restMaskDouble,tag_double,z2Ptr,highway::ConvertTo(tag_long,indices));

        const VectorLong ownedStateJ = highway::MaskedGatherIndex(restMaskLong,tag_long,ownedStatePtr2,highway::ConvertTo(tag_long,indices));
        ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
      }
      else {
        x2 = highway::GatherIndex(tag_double,x2Ptr,highway::ConvertTo(tag_long,indices));
        y2 = highway::GatherIndex(tag_double,y2Ptr,highway::ConvertTo(tag_long,indices));
        z2 = highway::GatherIndex(tag_double,z2Ptr,highway::ConvertTo(tag_long,indices));

        const VectorLong ownedStateJ = highway::GatherIndex(tag_long,ownedStatePtr2,highway::ConvertTo(tag_long,indices));
        ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
      }
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {

      const auto restMaskDouble = restMasksDouble[remainder ? rest-1 : _vecLengthDouble/2-1];
      const auto restMaskLong = restMasksLong[remainder ? rest-1 : _vecLengthDouble/2-1];

      VectorLong ownedStateJ = highway::IfThenElseZero(restMaskLong, highway::GatherIndex(tag_long,ownedStatePtr2,highway::ConvertTo(tag_long,indices)));
      x2 = highway::IfThenElseZero(restMaskDouble, highway::GatherIndex(tag_double,x2Ptr,highway::ConvertTo(tag_long,indices)));
      y2 = highway::IfThenElseZero(restMaskDouble, highway::GatherIndex(tag_double,y2Ptr,highway::ConvertTo(tag_long,indices)));
      z2 = highway::IfThenElseZero(restMaskDouble, highway::GatherIndex(tag_double,z2Ptr,highway::ConvertTo(tag_long,indices)));

      // "broadcast" lower half to upper half
      ownedStateJ = highway::ConcatLowerLower(tag_long, ownedStateJ, ownedStateJ);
      ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
      x2 = highway::ConcatLowerLower(tag_double, x2, x2);
      y2 = highway::ConcatLowerLower(tag_double, y2, y2);
      z2 = highway::ConcatLowerLower(tag_double, z2, z2);
    }
  }

  template <bool remainderI, bool remainderJ>
  inline void fillPhysicsRegisters(const size_t *const typeID1Ptr, const size_t *const typeID2Ptr,
                                   VectorDouble& epsilon24s, VectorDouble& sigmaSquareds, VectorDouble& shift6s,
                                   const unsigned int restI, const unsigned int restJ) {
    double epsilons[_vecLengthDouble] = {0.};
    double sigmas[_vecLengthDouble] = {0.};
    double shifts[_vecLengthDouble] = {0.};

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      for (int n = 0; n < (remainderJ ? restJ : _vecLengthDouble); ++n) {
        epsilons[n] = _PPLibrary->getMixing24Epsilon(*typeID1Ptr, *(typeID2Ptr + n));
        sigmas[n] = _PPLibrary->getMixingSigmaSquared(*typeID1Ptr, *(typeID2Ptr + n));
        if constexpr (applyShift) {
          shifts[n] = _PPLibrary->getMixingShift6(*typeID1Ptr, *(typeID2Ptr + n));
        }
      }
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      for (int i = 0; i < (remainderI ? 1 : 2); ++i) {
        for (int j = 0; j < (remainderJ ? restJ : _vecLengthDouble/2); ++j) {
          epsilons[i* (_vecLengthDouble/2) + j] = _PPLibrary->getMixing24Epsilon(*(typeID1Ptr+i), *(typeID2Ptr + j));
          sigmas[i* (_vecLengthDouble/2) + j] = _PPLibrary->getMixingSigmaSquared(*(typeID1Ptr+i), *(typeID2Ptr + j));

          if constexpr (applyShift) {
            shifts[i* (_vecLengthDouble/2) + j] = _PPLibrary->getMixingShift6(*(typeID1Ptr+i), *(typeID2Ptr + j));
          }
        }
      }
    }

    epsilon24s = highway::LoadU(tag_double, epsilons);
    sigmaSquareds = highway::LoadU(tag_double, sigmas);
    if constexpr (applyShift) {
      shift6s = highway::LoadU(tag_double, shifts);
    }
  }

  // not sure what the changes here need to be for the gather scatter approach
  template <bool remainderI, bool remainderJ>
  inline void fillPhysicsRegistersGS(const VectorDouble indices,const size_t *const typeID1Ptr, const size_t *const typeID2Ptr,
                                     VectorDouble& epsilon24s, VectorDouble& sigmaSquareds, VectorDouble& shift6s,
                                     const unsigned int restI, const unsigned int restJ) {
    double epsilons[_vecLengthDouble] = {0.};
    double sigmas[_vecLengthDouble] = {0.};
    double shifts[_vecLengthDouble] = {0.};

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      for (int n = 0; n < (remainderJ ? restJ : _vecLengthDouble); ++n) {
        epsilons[n] = _PPLibrary->getMixing24Epsilon(*typeID1Ptr, *(typeID2Ptr + n));
        sigmas[n] = _PPLibrary->getMixingSigmaSquared(*typeID1Ptr, *(typeID2Ptr + n));
        if constexpr (applyShift) {
          shifts[n] = _PPLibrary->getMixingShift6(*typeID1Ptr, *(typeID2Ptr + n));
        }
      }
    }
    else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      for (int i = 0; i < (remainderI ? 1 : 2); ++i) {
        for (int j = 0; j < (remainderJ ? restJ : _vecLengthDouble/2); ++j) {
          epsilons[i* (_vecLengthDouble/2) + j] = _PPLibrary->getMixing24Epsilon(*(typeID1Ptr+i), *(typeID2Ptr + j));
          sigmas[i* (_vecLengthDouble/2) + j] = _PPLibrary->getMixingSigmaSquared(*(typeID1Ptr+i), *(typeID2Ptr + j));

          if constexpr (applyShift) {
            shifts[i* (_vecLengthDouble/2) + j] = _PPLibrary->getMixingShift6(*(typeID1Ptr+i), *(typeID2Ptr + j));
          }
        }
      }
    }

    epsilon24s = highway::LoadU(tag_double, epsilons);
    sigmaSquareds = highway::LoadU(tag_double, sigmas);
    if constexpr (applyShift) {
      shift6s = highway::LoadU(tag_double, shifts);
    }
  }
  template <bool newton3, bool remainderI, bool remainderJ>
  inline void SoAKernelSmoothGS(const VectorDouble gatheredIndices, const size_t j, const MaskDouble& ownedMaskI, const int64_t *const __restrict ownedStatePtr2,
                                const VectorDouble &x1, const VectorDouble &y1, const VectorDouble &z1, const double *const __restrict x2Ptr,
                                const double *const __restrict y2Ptr, const double *const __restrict z2Ptr,
                                double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                                double *const __restrict fz2Ptr, const size_t *const typeID1Ptr, const size_t *const typeID2Ptr,
                                VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc, VectorDouble &virialSumX, VectorDouble& virialSumY,
                                VectorDouble &virialSumZ, VectorDouble &uPotSum, const unsigned int restI, const unsigned int restJ) {

    VectorDouble epsilon24s = _epsilon24;
    VectorDouble sigmaSquareds = _sigmaSquared;
    VectorDouble shift6s = _shift6;

    if constexpr (useMixing) {
      fillPhysicsRegistersGS<remainderI, remainderJ>(gatheredIndices,typeID1Ptr, typeID2Ptr, epsilon24s, sigmaSquareds, shift6s, restI, restJ);
    }

    VectorDouble x2;
    VectorDouble y2;
    VectorDouble z2;
    VectorDouble ownedStateJDouble;

    fillJRegistersGS<remainderJ>(gatheredIndices,j, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, x2, y2, z2, ownedStateJDouble, restJ);


    // distance calculations
    const auto drX = x1 - x2;
    const auto drY = y1 - y2;
    const auto drZ = z1 - z2;

    const auto drX2 = drX * drX;
    const auto drY2 = drY * drY;
    const auto drZ2 = drZ * drZ;

    const auto dr2 = drX2 + drY2 + drZ2;

    const auto dummyMask = highway::And(ownedMaskI, highway::Ne(ownedStateJDouble, _ownedStateDummy));

    /*
    const auto cutoffMask = highway::Le(dr2, _cutoffSquared);
    const auto innerCutoffMask = highway::Ge(dr2, _innerCutoff);
    const auto zeroMask = highway::Ne(dr2, _zeroDouble);
    const auto cutoffDummyMask = highway::And(zeroMask, highway::And(cutoffMask, dummyMask));
    const auto innerCutoffDummyMask = highway::And(zeroMask, highway::And(innerCutoffMask, dummyMask));
     */

    //compute smoothing factor
    const auto drtrue = highway::Sqrt(dr2);
    const auto temp = drtrue - _innerCutoff;
    const auto temp2 = temp*temp;
    const auto twodrtrue = _twoDouble * drtrue;
    const auto tripleOuterMinusInnerAoSLocalm2drt = _tripleOuterMinusInner - twodrtrue;
    const auto numerator = temp2 * tripleOuterMinusInnerAoSLocalm2drt;
    const auto fraction = numerator * _cutoffDiffCubedInv;
    const auto smoothingTerm = _oneDouble - fraction;
    const auto smoothingMasked = highway::IfThenElse(dummyMask,smoothingTerm,_oneDouble);

    // compute LJ Potential
    const auto invDr2 = _oneDouble / dr2;
    const auto lj2 = sigmaSquareds * invDr2;
    const auto lj4 = lj2 * lj2;
    const auto lj6 = lj2 * lj4;
    const auto lj12 = lj6 * lj6;
    const auto lj12m6 = lj12 - lj6;
    const auto lj12m6alj12 = lj12m6 + lj12;
    const auto lj12m6alj12e = lj12m6alj12 * epsilon24s;
    const auto fac = lj12m6alj12e * invDr2 * smoothingMasked;

    const auto facMasked = highway::IfThenElseZero(dummyMask, fac);

    const VectorDouble fx = drX * facMasked;
    const VectorDouble fy = drY * facMasked;
    const VectorDouble fz = drZ * facMasked;

    fxAcc = fxAcc + fx;
    fyAcc = fyAcc + fy;
    fzAcc = fzAcc + fz;

    if constexpr (newton3) {
      handleNewton3ReductionGS<remainderJ>(gatheredIndices,fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, restJ);
    }

    if constexpr (calculateGlobals) {
      auto virialX = fx * drX;
      auto virialY = fy * drY;
      auto virialZ = fz * drZ;

      auto uPot = highway::MulAdd(epsilon24s, lj12m6,shift6s);
      uPot = highway::Mul(uPot, smoothingTerm);
      auto uPotMasked = highway::IfThenElseZero(dummyMask, uPot);

      auto energyFactor = highway::IfThenElse(dummyMask, _oneDouble, _zeroDouble);

      if constexpr (newton3) {
        energyFactor = energyFactor + highway::IfThenElse(dummyMask, _oneDouble, _zeroDouble);
      }

      uPotSum = highway::MulAdd(energyFactor, uPotMasked, uPotSum);
      virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
      virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
      virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
    }
  }

  template <bool newton3, bool remainderI, bool remainderJ>
  inline void SoAKernelGS(const VectorDouble gatheredIndices, const size_t j, const MaskDouble& ownedMaskI, const int64_t *const __restrict ownedStatePtr2,
                          const VectorDouble &x1, const VectorDouble &y1, const VectorDouble &z1, const double *const __restrict x2Ptr,
                          const double *const __restrict y2Ptr, const double *const __restrict z2Ptr,
                          double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                          double *const __restrict fz2Ptr, const size_t *const typeID1Ptr, const size_t *const typeID2Ptr,
                          VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc, VectorDouble &virialSumX, VectorDouble& virialSumY,
                          VectorDouble &virialSumZ, VectorDouble &uPotSum, const unsigned int restI, const unsigned int restJ) {

    VectorDouble epsilon24s = _epsilon24;
    VectorDouble sigmaSquareds = _sigmaSquared;
    VectorDouble shift6s = _shift6;

    if constexpr (useMixing) {
      fillPhysicsRegistersGS<remainderI, remainderJ>(gatheredIndices,typeID1Ptr, typeID2Ptr, epsilon24s, sigmaSquareds, shift6s, restI, restJ);
    }

    VectorDouble x2;
    VectorDouble y2;
    VectorDouble z2;
    VectorDouble ownedStateJDouble;

    fillJRegistersGS<remainderJ>(gatheredIndices,j, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, x2, y2, z2, ownedStateJDouble, restJ);

    // distance calculations
    const auto drX = x1 - x2;
    const auto drY = y1 - y2;
    const auto drZ = z1 - z2;

    const auto drX2 = drX * drX;
    const auto drY2 = drY * drY;
    const auto drZ2 = drZ * drZ;

    const auto dr2 = drX2 + drY2 + drZ2;

    const auto dummyMask = highway::And(ownedMaskI, highway::Ne(ownedStateJDouble, _ownedStateDummy));

    /*
    const auto cutoffMask = highway::Le(dr2, _cutoffSquared);
    const auto zeroMask = highway::Ne(dr2, _zeroDouble);
    const auto cutoffDummyMask = highway::And(zeroMask, highway::And(cutoffMask, dummyMask));
     */

    // compute LJ Potential
    const auto invDr2 = _oneDouble / dr2;
    const auto lj2 = sigmaSquareds * invDr2;
    const auto lj4 = lj2 * lj2;
    const auto lj6 = lj2 * lj4;
    const auto lj12 = lj6 * lj6;
    const auto lj12m6 = lj12 - lj6;
    const auto lj12m6alj12 = lj12m6 + lj12;
    const auto lj12m6alj12e = lj12m6alj12 * epsilon24s;
    const auto fac = lj12m6alj12e * invDr2;

    const auto facMasked = highway::IfThenElseZero(dummyMask, fac);

    const VectorDouble fx = drX * facMasked;
    const VectorDouble fy = drY * facMasked;
    const VectorDouble fz = drZ * facMasked;

    fxAcc = fxAcc + fx;
    fyAcc = fyAcc + fy;
    fzAcc = fzAcc + fz;

    if constexpr (newton3) {
      handleNewton3ReductionGS<remainderJ>(gatheredIndices,fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, restJ);
    }

    if constexpr (calculateGlobals) {
      auto virialX = fx * drX;
      auto virialY = fy * drY;
      auto virialZ = fz * drZ;

      auto uPot = highway::MulAdd(epsilon24s, lj12m6,shift6s);
      auto uPotMasked = highway::IfThenElseZero(dummyMask, uPot);

      auto energyFactor = highway::IfThenElse(dummyMask, _oneDouble, _zeroDouble);

      if constexpr (newton3) {
        energyFactor = energyFactor + highway::IfThenElse(dummyMask, _oneDouble, _zeroDouble);
      }

      uPotSum = highway::MulAdd(energyFactor, uPotMasked, uPotSum);
      virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
      virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
      virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
    }
  }

  template <bool newton3, bool remainderI, bool remainderJ>
  inline void SoAKernel(const size_t j, const MaskDouble& ownedMaskI, const int64_t *const __restrict ownedStatePtr2,
                        const VectorDouble &x1, const VectorDouble &y1, const VectorDouble &z1, const double *const __restrict x2Ptr,
                        const double *const __restrict y2Ptr, const double *const __restrict z2Ptr,
                        double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                        double *const __restrict fz2Ptr, const size_t *const typeID1Ptr, const size_t *const typeID2Ptr,
                        VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc, VectorDouble &virialSumX, VectorDouble& virialSumY,
                        VectorDouble &virialSumZ, VectorDouble &uPotSum, const unsigned int restI, const unsigned int restJ) {

    VectorDouble epsilon24s = _epsilon24;
    VectorDouble sigmaSquareds = _sigmaSquared;
    VectorDouble shift6s = _shift6;

    if constexpr (useMixing) {
      fillPhysicsRegisters<remainderI, remainderJ>(typeID1Ptr, typeID2Ptr, epsilon24s, sigmaSquareds, shift6s, restI, restJ);
    }

    VectorDouble x2;
    VectorDouble y2;
    VectorDouble z2;
    VectorDouble ownedStateJDouble;

    fillJRegisters<remainderJ>(j, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, x2, y2, z2, ownedStateJDouble, restJ);

    // distance calculations
    const auto drX = x1 - x2;
    const auto drY = y1 - y2;
    const auto drZ = z1 - z2;

    const auto drX2 = drX * drX;
    const auto drY2 = drY * drY;
    const auto drZ2 = drZ * drZ;

    const auto dr2 = drX2 + drY2 + drZ2;

    const auto cutoffMask = highway::Le(dr2, _cutoffSquared);
    const auto zeroMask = highway::Ne(dr2, _zeroDouble);
    const auto dummyMask = highway::And(ownedMaskI, highway::Ne(ownedStateJDouble, _ownedStateDummy));
    const auto cutoffDummyMask = highway::And(zeroMask, highway::And(cutoffMask, dummyMask));

    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
      return;
    }

    // compute LJ Potential
    const auto invDr2 = _oneDouble / dr2;
    const auto lj2 = sigmaSquareds * invDr2;
    const auto lj4 = lj2 * lj2;
    const auto lj6 = lj2 * lj4;
    const auto lj12 = lj6 * lj6;
    const auto lj12m6 = lj12 - lj6;
    const auto lj12m6alj12 = lj12m6 + lj12;
    const auto lj12m6alj12e = lj12m6alj12 * epsilon24s;
    const auto fac = lj12m6alj12e * invDr2;

    const auto facMasked = highway::IfThenElseZero(cutoffDummyMask, fac);

    const VectorDouble fx = drX * facMasked;
    const VectorDouble fy = drY * facMasked;
    const VectorDouble fz = drZ * facMasked;

    fxAcc = fxAcc + fx;
    fyAcc = fyAcc + fy;
    fzAcc = fzAcc + fz;

    if constexpr (newton3) {
      handleNewton3Reduction<remainderJ>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, restJ);
    }

    if constexpr (calculateGlobals) {
      auto virialX = fx * drX;
      auto virialY = fy * drY;
      auto virialZ = fz * drZ;

      auto uPot = highway::MulAdd(epsilon24s, lj12m6, shift6s);
      auto uPotMasked = highway::IfThenElseZero(cutoffDummyMask, uPot);

      auto energyFactor = highway::IfThenElse(dummyMask, _oneDouble, _zeroDouble);

      if constexpr (newton3) {
        energyFactor = energyFactor + highway::IfThenElse(dummyMask, _oneDouble, _zeroDouble);
      }

      uPotSum = highway::MulAdd(energyFactor, uPotMasked, uPotSum);
      virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
      virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
      virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
    }
  }
  template <bool newton3, bool remainderI, bool remainderJ>
  inline void SoAKernelSmooth(const size_t j, const MaskDouble& ownedMaskI, const int64_t *const __restrict ownedStatePtr2,
                              const VectorDouble &x1, const VectorDouble &y1, const VectorDouble &z1, const double *const __restrict x2Ptr,
                              const double *const __restrict y2Ptr, const double *const __restrict z2Ptr,
                              double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                              double *const __restrict fz2Ptr, const size_t *const typeID1Ptr, const size_t *const typeID2Ptr,
                              VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc, VectorDouble &virialSumX, VectorDouble& virialSumY,
                              VectorDouble &virialSumZ, VectorDouble &uPotSum, const unsigned int restI, const unsigned int restJ) {

    VectorDouble epsilon24s = _epsilon24;
    VectorDouble sigmaSquareds = _sigmaSquared;
    VectorDouble shift6s = _shift6;

    if constexpr (useMixing) {
      fillPhysicsRegisters<remainderI, remainderJ>(typeID1Ptr, typeID2Ptr, epsilon24s, sigmaSquareds, shift6s, restI, restJ);
    }

    VectorDouble x2;
    VectorDouble y2;
    VectorDouble z2;
    VectorDouble ownedStateJDouble;

    fillJRegisters<remainderJ>(j, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, x2, y2, z2, ownedStateJDouble, restJ);

    // distance calculations
    const auto drX = x1 - x2;
    const auto drY = y1 - y2;
    const auto drZ = z1 - z2;

    const auto drX2 = drX * drX;
    const auto drY2 = drY * drY;
    const auto drZ2 = drZ * drZ;

    const auto dr2 = drX2 + drY2 + drZ2;

    const auto cutoffMask = highway::Le(dr2, _cutoffSquared);
    const auto innerCutoffMask = highway::Ge(dr2, _innerCutoffSquared);
    const auto zeroMask = highway::Ne(dr2, _zeroDouble);
    const auto dummyMask = highway::And(ownedMaskI, highway::Ne(ownedStateJDouble, _ownedStateDummy));
    const auto cutoffDummyMask = highway::And(zeroMask, highway::And(cutoffMask, dummyMask));
    const auto innerCutoffDummyMask = highway::And(zeroMask, highway::And(innerCutoffMask, dummyMask));

    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
      return;
    }
    //compute smoothing factor
    const auto drtrue = highway::Sqrt(dr2);
    const auto temp = drtrue - _innerCutoff;
    const auto temp2 = temp*temp;
    const auto twodrtrue = _twoDouble * drtrue;
    const auto tripleOuterMinusInnerAoSLocalm2drt = _tripleOuterMinusInner - twodrtrue;
    const auto numerator = temp2 * tripleOuterMinusInnerAoSLocalm2drt;
    const auto fraction = numerator * _cutoffDiffCubedInv;
    const auto smoothingTerm = _oneDouble - fraction;
    const auto smoothingMasked = highway::IfThenElse(innerCutoffDummyMask,smoothingTerm,_oneDouble);


    // compute LJ Potential
    const auto invDr2 = _oneDouble / dr2;
    const auto lj2 = sigmaSquareds * invDr2;
    const auto lj4 = lj2 * lj2;
    const auto lj6 = lj2 * lj4;
    const auto lj12 = lj6 * lj6;
    const auto lj12m6 = lj12 - lj6;
    const auto lj12m6alj12 = lj12m6 + lj12;
    const auto lj12m6alj12e = lj12m6alj12 * epsilon24s;
    const auto fac = lj12m6alj12e * invDr2 * smoothingMasked;


    const auto facMasked = highway::IfThenElseZero(cutoffDummyMask, fac);

    const VectorDouble fx = drX * facMasked;
    const VectorDouble fy = drY * facMasked;
    const VectorDouble fz = drZ * facMasked;

    fxAcc = fxAcc + fx;
    fyAcc = fyAcc + fy;
    fzAcc = fzAcc + fz;

    if constexpr (newton3) {
      handleNewton3Reduction<remainderJ>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, restJ);
    }

    if constexpr (calculateGlobals) {
      auto virialX = fx * drX;
      auto virialY = fy * drY;
      auto virialZ = fz * drZ;


      auto uPot = highway::MulAdd(epsilon24s, lj12m6,shift6s);
      uPot = highway::Mul(uPot, smoothingMasked);
      auto uPotMasked = highway::IfThenElseZero(cutoffDummyMask, uPot);

      auto energyFactor = highway::IfThenElse(dummyMask, _oneDouble, _zeroDouble);

      if constexpr (newton3) {
        energyFactor = energyFactor + highway::IfThenElse(dummyMask, _oneDouble, _zeroDouble);
      }

      uPotSum = highway::MulAdd(energyFactor, uPotMasked, uPotSum);
      virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
      virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
      virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
    }
  }
 public:
  // clang-format off
            /**
             * @copydoc Functor::SoAFunctorVerlet(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
             * @note If you want to parallelize this by openmp, please ensure that there
             * are no dependencies, i.e. introduce colors and specify iFrom and iTo accordingly.
             */
  // clang-format on
  inline void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                               const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                               bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      //SoAFunctorSingleImplGS<true>(soa);
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      //SoAFunctorSingleImplGS<false>(soa);
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

    const auto *const __restrict xPtr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict fxPtr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeIDPtr = soa.template begin<Particle::AttributeNames::typeId>();

    VectorDouble virialSumX = highway::Zero(tag_double);
    VectorDouble virialSumY = highway::Zero(tag_double);
    VectorDouble virialSumZ = highway::Zero(tag_double);
    VectorDouble uPotSum = highway::Zero(tag_double);
    VectorDouble fxAcc = highway::Zero(tag_double);
    VectorDouble fyAcc = highway::Zero(tag_double);
    VectorDouble fzAcc = highway::Zero(tag_double);

    const VectorDouble x1 = highway::Set(tag_double, xPtr[indexFirst]);
    const VectorDouble y1 = highway::Set(tag_double, yPtr[indexFirst]);
    const VectorDouble z1 = highway::Set(tag_double, zPtr[indexFirst]);
    const int64_t ownedI = static_cast<int64_t>(ownedStatePtr[indexFirst]);
    const VectorDouble ownedStateI = highway::Set(tag_double, static_cast<double>(ownedI));
    const MaskDouble ownedMaskI = highway::Ne(ownedStateI, _zeroDouble);

    alignas(64) std::array<double, _vecLengthDouble> x2Tmp{};
    alignas(64) std::array<double, _vecLengthDouble> y2Tmp{};
    alignas(64) std::array<double, _vecLengthDouble> z2Tmp{};
    alignas(64) std::array<double, _vecLengthDouble> fx2Tmp{};
    alignas(64) std::array<double, _vecLengthDouble> fy2Tmp{};
    alignas(64) std::array<double, _vecLengthDouble> fz2Tmp{};
    alignas(64) std::array<size_t, _vecLengthDouble> typeID2Tmp{};
    alignas(64) std::array<autopas::OwnershipState, _vecLengthDouble> ownedStates2Tmp{};

    size_t j = 0;

    for (; j < (neighborList.size() & ~(_vecLengthDouble - 1)); j += _vecLengthDouble) {

      // load neighbor particles in consecutive array
      for (long vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
        x2Tmp[vecIndex] = xPtr[neighborList[j + vecIndex]];
        y2Tmp[vecIndex] = yPtr[neighborList[j + vecIndex]];
        z2Tmp[vecIndex] = zPtr[neighborList[j + vecIndex]];
        if constexpr (newton3) {
          fx2Tmp[vecIndex] = fxPtr[neighborList[j + vecIndex]];
          fy2Tmp[vecIndex] = fyPtr[neighborList[j + vecIndex]];
          fz2Tmp[vecIndex] = fzPtr[neighborList[j + vecIndex]];
        }
        typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
        ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
      }

      SoAKernelSmooth<newton3, false, false>(0, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()), x1, y1, z1,
                                       x2Tmp.data(), y2Tmp.data(), z2Tmp.data(), fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(),
                                       &typeIDPtr[indexFirst], typeID2Tmp.data(), fxAcc, fyAcc, fzAcc, virialSumX,
                                       virialSumY, virialSumZ, uPotSum, 0, 0);

      if constexpr (newton3) {
        for (size_t vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
          fxPtr[neighborList[j + vecIndex]] = fx2Tmp[vecIndex];
          fyPtr[neighborList[j + vecIndex]] = fy2Tmp[vecIndex];
          fzPtr[neighborList[j + vecIndex]] = fz2Tmp[vecIndex];
        }
      }
    }

    const int rest = static_cast<int>(neighborList.size() & (_vecLengthDouble - 1));

    if (rest > 0) {
      for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
        x2Tmp[vecIndex] = xPtr[neighborList[j + vecIndex]];
        y2Tmp[vecIndex] = yPtr[neighborList[j + vecIndex]];
        z2Tmp[vecIndex] = zPtr[neighborList[j + vecIndex]];
        if constexpr (newton3) {
          fx2Tmp[vecIndex] = fxPtr[neighborList[j + vecIndex]];
          fy2Tmp[vecIndex] = fyPtr[neighborList[j + vecIndex]];
          fz2Tmp[vecIndex] = fzPtr[neighborList[j + vecIndex]];
        }
        typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
        ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
      }

      SoAKernelSmooth<newton3, false, true>(0, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()), x1, y1, z1,
                                      x2Tmp.data(), y2Tmp.data(), z2Tmp.data(), fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(),
                                      &typeIDPtr[indexFirst], typeID2Tmp.data(), fxAcc, fyAcc, fzAcc, virialSumX,
                                     virialSumY, virialSumZ, uPotSum, 0, rest);

      if constexpr (newton3) {

        for (long vecIndex = 0; vecIndex < _vecLengthDouble && vecIndex < rest; ++vecIndex) {
          fxPtr[neighborList[j + vecIndex]] = fx2Tmp[vecIndex];
          fyPtr[neighborList[j + vecIndex]] = fy2Tmp[vecIndex];
          fzPtr[neighborList[j + vecIndex]] = fz2Tmp[vecIndex];
        }
      }
    }

    fxPtr[indexFirst] += highway::ReduceSum(tag_double, fxAcc);
    fyPtr[indexFirst] += highway::ReduceSum(tag_double, fyAcc);
    fzPtr[indexFirst] += highway::ReduceSum(tag_double, fzAcc);

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_num_threads();

      double globals[] = {
          highway::ReduceSum(tag_double, virialSumX),
          highway::ReduceSum(tag_double, virialSumY),
          highway::ReduceSum(tag_double, virialSumZ),
          highway::ReduceSum(tag_double, uPotSum)
      };

      double factor = 1.;
      factor *= newton3 ? .5 : 1.;
      _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
      _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
      _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
      _aosThreadData[threadnum].uPotSum += globals[3] * factor;
    }
  }

  template <bool newton3>
  inline void SoAFunctorVerletImplGS(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                                   const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();
    if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
      return;
    }

    const auto *const __restrict xPtr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict fxPtr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle::AttributeNames::forceZ>();

    const auto *const __restrict typeIDPtr = soa.template begin<Particle::AttributeNames::typeId>();

    VectorDouble virialSumX = highway::Zero(tag_double);
    VectorDouble virialSumY = highway::Zero(tag_double);
    VectorDouble virialSumZ = highway::Zero(tag_double);
    VectorDouble uPotSum = highway::Zero(tag_double);
    VectorDouble fxAcc = highway::Zero(tag_double);
    VectorDouble fyAcc = highway::Zero(tag_double);
    VectorDouble fzAcc = highway::Zero(tag_double);

    const VectorDouble x1 = highway::Set(tag_double, xPtr[indexFirst]);
    const VectorDouble y1 = highway::Set(tag_double, yPtr[indexFirst]);
    const VectorDouble z1 = highway::Set(tag_double, zPtr[indexFirst]);
    const int64_t ownedI = static_cast<int64_t>(ownedStatePtr[indexFirst]);
    const VectorDouble ownedStateI = highway::Set(tag_double, static_cast<double>(ownedI));
    const MaskDouble ownedMaskI = highway::Ne(ownedStateI, _zeroDouble);

    alignas(64) std::array<double, _vecLengthDouble> x2Tmp{};
    alignas(64) std::array<double, _vecLengthDouble> y2Tmp{};
    alignas(64) std::array<double, _vecLengthDouble> z2Tmp{};
    alignas(64) std::array<double, _vecLengthDouble> fx2Tmp{};
    alignas(64) std::array<double, _vecLengthDouble> fy2Tmp{};
    alignas(64) std::array<double, _vecLengthDouble> fz2Tmp{};
    alignas(64) std::array<size_t, _vecLengthDouble> typeID2Tmp{};
    alignas(64) std::array<autopas::OwnershipState, _vecLengthDouble> ownedStates2Tmp{};

    size_t j = 0;

    for (; j < (neighborList.size() & ~(_vecLengthDouble - 1)); j += _vecLengthDouble) {

      // load neighbor particles in consecutive array
      for (long vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
        x2Tmp[vecIndex] = xPtr[neighborList[j + vecIndex]];
        y2Tmp[vecIndex] = yPtr[neighborList[j + vecIndex]];
        z2Tmp[vecIndex] = zPtr[neighborList[j + vecIndex]];
        if constexpr (newton3) {
          fx2Tmp[vecIndex] = fxPtr[neighborList[j + vecIndex]];
          fy2Tmp[vecIndex] = fyPtr[neighborList[j + vecIndex]];
          fz2Tmp[vecIndex] = fzPtr[neighborList[j + vecIndex]];
        }
        typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
        ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
      }

      SoAKernelSmooth<newton3, false, false>(0, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()), x1, y1, z1,
                                       x2Tmp.data(), y2Tmp.data(), z2Tmp.data(), fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(),
                                       &typeIDPtr[indexFirst], typeID2Tmp.data(), fxAcc, fyAcc, fzAcc, virialSumX,
                                       virialSumY, virialSumZ, uPotSum, 0, 0);

      if constexpr (newton3) {

        for (size_t vecIndex = 0; vecIndex < _vecLengthDouble; ++vecIndex) {
          fxPtr[neighborList[j + vecIndex]] = fx2Tmp[vecIndex];
          fyPtr[neighborList[j + vecIndex]] = fy2Tmp[vecIndex];
          fzPtr[neighborList[j + vecIndex]] = fz2Tmp[vecIndex];
        }
      }
    }

    const int rest = static_cast<int>(neighborList.size() & (_vecLengthDouble - 1));

    if (rest > 0) {
      for (size_t vecIndex = 0; vecIndex < rest; ++vecIndex) {
        x2Tmp[vecIndex] = xPtr[neighborList[j + vecIndex]];
        y2Tmp[vecIndex] = yPtr[neighborList[j + vecIndex]];
        z2Tmp[vecIndex] = zPtr[neighborList[j + vecIndex]];
        if constexpr (newton3) {
          fx2Tmp[vecIndex] = fxPtr[neighborList[j + vecIndex]];
          fy2Tmp[vecIndex] = fyPtr[neighborList[j + vecIndex]];
          fz2Tmp[vecIndex] = fzPtr[neighborList[j + vecIndex]];
        }
        typeID2Tmp[vecIndex] = typeIDPtr[neighborList[j + vecIndex]];
        ownedStates2Tmp[vecIndex] = ownedStatePtr[neighborList[j + vecIndex]];
      }

      SoAKernelSmooth<newton3, false, true>(0, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStates2Tmp.data()), x1, y1, z1,
                                      x2Tmp.data(), y2Tmp.data(), z2Tmp.data(), fx2Tmp.data(), fy2Tmp.data(), fz2Tmp.data(),
                                      &typeIDPtr[indexFirst], typeID2Tmp.data(), fxAcc, fyAcc, fzAcc, virialSumX,
                                      virialSumY, virialSumZ, uPotSum, 0, rest);

      if constexpr (newton3) {

        for (long vecIndex = 0; vecIndex < _vecLengthDouble && vecIndex < rest; ++vecIndex) {
          fxPtr[neighborList[j + vecIndex]] = fx2Tmp[vecIndex];
          fyPtr[neighborList[j + vecIndex]] = fy2Tmp[vecIndex];
          fzPtr[neighborList[j + vecIndex]] = fz2Tmp[vecIndex];
        }
      }
    }

    fxPtr[indexFirst] += highway::ReduceSum(tag_double, fxAcc);
    fyPtr[indexFirst] += highway::ReduceSum(tag_double, fyAcc);
    fzPtr[indexFirst] += highway::ReduceSum(tag_double, fzAcc);

    if constexpr (calculateGlobals) {
      const int threadnum = autopas::autopas_get_num_threads();

      double globals[] = {
          highway::ReduceSum(tag_double, virialSumX),
          highway::ReduceSum(tag_double, virialSumY),
          highway::ReduceSum(tag_double, virialSumZ),
          highway::ReduceSum(tag_double, uPotSum)
      };

      double factor = 1.;
      factor *= newton3 ? .5 : 1.;
      _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
      _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
      _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
      _aosThreadData[threadnum].uPotSum += globals[3] * factor;
    }
  }

 public:
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
  void initTraversal() final {
    _uPotSum = 0.;
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
        _uPotSum += _aosThreadData[i].uPotSum;
        _virialSum += _aosThreadData[i].virialSum;
      }
      if (not newton3) {
        // if the newton3 optimization is disabled we have added every energy contribution twice, so we divide by 2
        // here.
        _uPotSum *= 0.5;
        _virialSum *= 0.5;
      }
      // we have always calculated 6*upot, so we divide by 6 here!
      _uPotSum /= 6.;
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
    return _uPotSum;
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
      throw autopas::utils::ExceptionHandler::AutoPasException("Cannot get virial, because endTraversal was not called.");
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
    _epsilon24 = highway::Set(tag_double, epsilon24);
    _sigmaSquared = highway::Set(tag_double, sigmaSquare);
    if constexpr (applyShift) {
      _shift6 = highway::Set(tag_double, ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, highway::GetLane(_cutoffSquared)));
    } else {
      _shift6 = _zeroDouble;
    }
    _epsilon24AoS = epsilon24;
    _sigmaSquareAoS = sigmaSquare;
    if constexpr (applyShift) {
      _shift6AoS = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, _cutoffSquareAoS);
    } else {
      _shift6AoS = 0.;
    }
  }

 private:

  class AoSThreadData
  {
   public:
    AoSThreadData() : virialSum{0.,0.,0.}, uPotSum{0.} {}

    void setZero() {
      virialSum = {0.,0.,0.};
      uPotSum = 0.;
    }

    std::array<double, 3> virialSum {};
    double uPotSum {0};

   private:
    double _remainingTo64[4];
  };
  // static_assert(sizeof(AoSThreadData) & 64 == 0, "AoSThreadData has wrong size");

  // helper variables for the LJ-calculation
  const VectorDouble _zeroDouble {highway::Zero(tag_double)};
  const VectorLong _zeroLong {highway::Zero(tag_long)};
  const VectorDouble _oneDouble {highway::Set(tag_double, 1.)};
  const VectorDouble _twoDouble {highway::Set(tag_double, 2.)};
  const VectorLong _oneLong {highway::Set(tag_long, 1)};
  double asc[_vecLengthDouble];
  VectorDouble _ascendingIndices{};
  const VectorDouble _ownedStateDummy{highway::Zero(tag_double)};
  const VectorDouble _cutoffSquared {};
  const VectorDouble _cutoff {};
  const VectorDouble _innerCutoffSquared {};
  const VectorDouble _innerCutoff {};
  const VectorDouble _cutoffDiffCubedInv {};
  const VectorDouble _cutoffDiff {};
  const VectorDouble _tripleOuterMinusInner {};
  VectorDouble _shift6 {highway::Zero(tag_double)};
  VectorDouble _epsilon24 {highway::Zero(tag_double)};
  VectorDouble _sigmaSquared {highway::Zero(tag_double)};

  MaskDouble _remainderMask[_vecLengthDouble];
  MaskDouble restMasksDouble[_vecLengthDouble-1];
  MaskLong restMasksLong[_vecLengthDouble-1];

  const double _cutoffSquareAoS {0.};
  const double _cutoffAoS {0.};
  const double _innerCutoffSquareAoS {0.};
  const double _innerCutoffAoS {0.};
  const double _cutoffDiffCubedInvAoS {0.};
  const double _cutoffDiffAoS {0.};
  const double _tripleOuterMinusInnerAoS {0.};
  double _epsilon24AoS, _sigmaSquareAoS, _shift6AoS = 0.;
  ParticlePropertiesLibrary<double, size_t>* _PPLibrary = nullptr;
  double _uPotSum {0.};
  std::array<double, 3> _virialSum;
  std::vector<AoSThreadData> _aosThreadData;
  bool _postProcessed;
  public:
  duration<double, std::milli> timeSmoothAoS{0};
  duration<double, std::milli> timeSmoothSoA{0};
  duration<double, std::milli> timeNoSmoothSoA{0};
};
// } // Highway
} // mdLib
// HWY_AFTER_NAMESPACE();