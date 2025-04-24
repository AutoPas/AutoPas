/**
 * @file LJFunctorHWY.h
 *
 * @date 27.03.2024
 * @author Luis Gall
 */

#pragma once

#include <hwy/highway.h>

#include <algorithm>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"

namespace mdLib {

namespace highway = hwy::HWY_NAMESPACE;
/** Higwhay tag for full double register */
const highway::ScalableTag<double> tag_double;
/** Highway tag for full long register */
const highway::ScalableTag<int64_t> tag_long;
/** Nuber of double values in a full register */
const size_t _vecLengthDouble{highway::Lanes(tag_double)};
/** Type for a Double vector register */
using VectorDouble = decltype(highway::Zero(tag_double));
/** Type for a Long vector register */
using VectorLong = decltype(highway::Zero(tag_long));
/** Highway tag for a half-filled double register */
const highway::Half<highway::DFromV<VectorDouble>> tag_double_half;
/** Type for a Double Mask */
using MaskDouble = decltype(highway::FirstN(tag_double, 1));
/** Type for a Long Mask */
using MaskLong = decltype(highway::FirstN(tag_long, 2));
/** Vectorization Pattern Type */
using VectorizationPattern = autopas::VectorizationPatternOption::Value;

/**
 * A functor to handle lennard-jones interactions between two particles (molecules)
 * This functor uses the SIMD abstraction library Google Highway to provide architecture independent vectorization
 * @tparam Particle The type of particle
 * @tparam applyShift Switch for the lj potential to be truncated shifted
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 * @tparam countFLOPs counts FLOPs and hitrate. Not implemented for this functor. Please use the AutoVec functor.*/
template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>

class LJFunctorHWY
    : public autopas::PairwiseFunctor<Particle, LJFunctorHWY<Particle, applyShift, useMixing, useNewton3,
                                                             calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorHWY() = delete;

 private:
  /**
   * Internal, actual constructor
   * @param cutoff
   * @note param dummy unused, only there to make the signature different from the public constructor
   */
  explicit LJFunctorHWY(double cutoff, void * /*dummy*/)
      : autopas::PairwiseFunctor<Particle, LJFunctorHWY<Particle, applyShift, useMixing, useNewton3, calculateGlobals,
                                                        countFLOPs, relevantForTuning>>(cutoff),
        _cutoffSquared{highway::Set(tag_double, cutoff * cutoff)},
        _cutoffSquareAoS{cutoff * cutoff},
        _uPotSum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData{},
        _postProcessed{false} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      AutoPasLog(DEBUG, "Using LJFunctorHWY with countFLOPs but FLOP counting is not implemented.");
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
  explicit LJFunctorHWY(double cutoff) : LJFunctorHWY(cutoff, nullptr) {
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
  explicit LJFunctorHWY(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctorHWY(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  std::string getName() final { return "LJFunctorHWY"; }

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

    if (dr2 > _cutoffSquareAoS) {
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
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always do a newton3 like traversal of the soa.
   */
  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, const bool newton3) final {
    if (newton3) {
      SoAFunctorSingleImpl<true, VectorizationPattern::p1xVec>(soa);
    } else {
      SoAFunctorSingleImpl<false, VectorizationPattern::p1xVec>(soa);
    }
  }

  // clang-format off
  /**
  * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
  */
  // clang-format on
  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                             bool newton3) final {
    switch (_vecPattern) {
      case VectorizationPattern::p1xVec: {
        if (newton3) {
          SoAFunctorPairImpl<true, VectorizationPattern::p1xVec>(soa1, soa2);
        } else {
          SoAFunctorPairImpl<false, VectorizationPattern::p1xVec>(soa1, soa2);
        }
        break;
      }
      case VectorizationPattern::p2xVecDiv2: {
        if (newton3) {
          SoAFunctorPairImpl<true, VectorizationPattern::p2xVecDiv2>(soa1, soa2);
        } else {
          SoAFunctorPairImpl<false, VectorizationPattern::p2xVecDiv2>(soa1, soa2);
        }
        break;
      }
      case VectorizationPattern::pVecDiv2x2: {
        if (newton3) {
          SoAFunctorPairImpl<true, VectorizationPattern::pVecDiv2x2>(soa1, soa2);
        } else {
          SoAFunctorPairImpl<false, VectorizationPattern::pVecDiv2x2>(soa1, soa2);
        }
        break;
      }
      case VectorizationPattern::pVecx1: {
        if (newton3) {
          SoAFunctorPairImpl<true, VectorizationPattern::pVecx1>(soa1, soa2);
        } else {
          SoAFunctorPairImpl<false, VectorizationPattern::pVecx1>(soa1, soa2);
        }
        break;
      }
      default:
        break;
    }
  }

 private:
  template <bool reversed, VectorizationPattern vecPattern>
  inline bool checkFirstLoopCondition(const long i, const long stop) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      return reversed ? (long)i >= 0 : (i < stop);
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      return reversed ? (long)i >= 1 : (i < (stop - 1));
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      if constexpr (reversed) {
        return (long)i >= _vecLengthDouble / 2;
      } else {
        const long criterion = stop - _vecLengthDouble / 2 + 1;
        return i < criterion;
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      if constexpr (reversed) {
        return (long)i >= _vecLengthDouble;
      } else {
        const long criterion = stop - _vecLengthDouble + 1;
        return i < criterion;
      }
    } else {
      return false;
    }
  }

  template <VectorizationPattern vecPattern>
  inline void decrementFirstLoop(size_t &i) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      --i;
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      i -= 2;
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      i -= _vecLengthDouble / 2;
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      i -= _vecLengthDouble;
    }
  }

  template <VectorizationPattern vecPattern>
  inline void incrementFirstLoop(unsigned int &i) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      ++i;
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      i += 2;
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      i += _vecLengthDouble / 2;
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      i += _vecLengthDouble;
    }
  }

  template <bool reversed>
  inline int obtainFirstLoopRest(const int i, const int stop) {
    return reversed ? (i < 0 ? 0 : i + 1) : stop - i;
  }

  template <VectorizationPattern vecPattern>
  inline bool checkSecondLoopCondition(const size_t i, const size_t j) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      return j < (i & ~(_vecLengthDouble - 1));
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      return j < (i & ~(_vecLengthDouble / 2 - 1));
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      return j < (i & ~(1));
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      return j < i;
    } else {
      return false;
    }
  }

  template <VectorizationPattern vecPattern>
  inline void incrementSecondLoop(unsigned int &j) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      j += _vecLengthDouble;
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      j += _vecLengthDouble / 2;
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      j += 2;
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      ++j;
    }
  }

  template <VectorizationPattern vecPattern>
  inline int obtainSecondLoopRest(const size_t i) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      return (int)(i & (_vecLengthDouble - 1));
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      return (int)(i & (_vecLengthDouble / 2 - 1));
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      return (int)(i & (1));
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      return 0;
    } else {
      return -1;
    }
  }

  template <bool remainder, bool reversed, VectorizationPattern vecPattern>
  inline void fillIRegisters(const size_t i, const double *const __restrict xPtr, const double *const __restrict yPtr,
                             const double *const __restrict zPtr,
                             const autopas::OwnershipState *const __restrict ownedStatePtr, VectorDouble &x1,
                             VectorDouble &y1, VectorDouble &z1, MaskDouble &ownedMaskI, const size_t restI) {
    VectorDouble ownedStateIDouble = _zeroDouble;

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      int64_t owned = static_cast<int64_t>(ownedStatePtr[i]);
      ownedStateIDouble = highway::Set(tag_double, static_cast<double>(owned));

      x1 = highway::Set(tag_double, xPtr[i]);
      y1 = highway::Set(tag_double, yPtr[i]);
      z1 = highway::Set(tag_double, zPtr[i]);
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
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
        int index = reversed ? i - 1 : i + 1;
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
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      const int index = reversed ? (remainder ? 0 : i - _vecLengthDouble / 2 + 1) : i;
      const int lanes = remainder ? restI : _vecLengthDouble / 2;

      const VectorLong ownedJ =
          highway::LoadN(tag_long, reinterpret_cast<const int64_t *>(&ownedStatePtr[index]), lanes);
      ownedStateIDouble = highway::ConvertTo(tag_double, ownedJ);

      x1 = highway::LoadN(tag_double, &xPtr[index], lanes);
      y1 = highway::LoadN(tag_double, &yPtr[index], lanes);
      z1 = highway::LoadN(tag_double, &zPtr[index], lanes);

      ownedStateIDouble = highway::ConcatLowerLower(tag_double, ownedStateIDouble, ownedStateIDouble);
      x1 = highway::ConcatLowerLower(tag_double, x1, x1);
      y1 = highway::ConcatLowerLower(tag_double, y1, y1);
      z1 = highway::ConcatLowerLower(tag_double, z1, z1);
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      int index = reversed ? (remainder ? 0 : i - _vecLengthDouble + 1) : i;

      if constexpr (remainder) {
        x1 = highway::LoadN(tag_double, &xPtr[index], restI);
        y1 = highway::LoadN(tag_double, &yPtr[index], restI);
        z1 = highway::LoadN(tag_double, &zPtr[index], restI);

        const VectorLong ownedStateI =
            highway::LoadN(tag_long, reinterpret_cast<const int64_t *>(&ownedStatePtr[index]), restI);
        ownedStateIDouble = highway::ConvertTo(tag_double, ownedStateI);
      } else {
        x1 = highway::LoadU(tag_double, &xPtr[index]);
        y1 = highway::LoadU(tag_double, &yPtr[index]);
        z1 = highway::LoadU(tag_double, &zPtr[index]);

        const VectorLong ownedStateI =
            highway::LoadU(tag_long, reinterpret_cast<const int64_t *>(&ownedStatePtr[index]));
        ownedStateIDouble = highway::ConvertTo(tag_double, ownedStateI);
      }
    }

    ownedMaskI = highway::Ne(ownedStateIDouble, _zeroDouble);
  }

  template <bool remainder, bool reversed, VectorizationPattern vecPattern>
  inline void handleNewton3Reduction(const VectorDouble &fx, const VectorDouble &fy, const VectorDouble &fz,
                                     double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                                     double *const __restrict fz2Ptr, const size_t i, const size_t j,
                                     const size_t rest) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      const VectorDouble fx2 =
          remainder ? highway::LoadN(tag_double, &fx2Ptr[j], rest) : highway::LoadU(tag_double, &fx2Ptr[j]);
      const VectorDouble fy2 =
          remainder ? highway::LoadN(tag_double, &fy2Ptr[j], rest) : highway::LoadU(tag_double, &fy2Ptr[j]);
      const VectorDouble fz2 =
          remainder ? highway::LoadN(tag_double, &fz2Ptr[j], rest) : highway::LoadU(tag_double, &fz2Ptr[j]);

      const VectorDouble fx2New = fx2 - fx;
      const VectorDouble fy2New = fy2 - fy;
      const VectorDouble fz2New = fz2 - fz;

      remainder ? highway::StoreN(fx2New, tag_double, &fx2Ptr[j], rest)
                : highway::StoreU(fx2New, tag_double, &fx2Ptr[j]);
      remainder ? highway::StoreN(fy2New, tag_double, &fy2Ptr[j], rest)
                : highway::StoreU(fy2New, tag_double, &fy2Ptr[j]);
      remainder ? highway::StoreN(fz2New, tag_double, &fz2Ptr[j], rest)
                : highway::StoreU(fz2New, tag_double, &fz2Ptr[j]);
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      const auto lowerFx = highway::LowerHalf(fx);
      const auto lowerFy = highway::LowerHalf(fy);
      const auto lowerFz = highway::LowerHalf(fz);

      const auto upperFx = highway::UpperHalf(tag_double_half, fx);
      const auto upperFy = highway::UpperHalf(tag_double_half, fy);
      const auto upperFz = highway::UpperHalf(tag_double_half, fz);

      const auto fxCombined = lowerFx + upperFx;
      const auto fyCombined = lowerFy + upperFy;
      const auto fzCombined = lowerFz + upperFz;

      const auto fxCombinedExt = highway::ZeroExtendVector(tag_double, fxCombined);
      const auto fyCombinedExt = highway::ZeroExtendVector(tag_double, fyCombined);
      const auto fzCombinedExt = highway::ZeroExtendVector(tag_double, fzCombined);

      const int lanes = remainder ? rest : _vecLengthDouble / 2;

      VectorDouble fx2 = highway::LoadN(tag_double, &fx2Ptr[j], lanes);
      VectorDouble fy2 = highway::LoadN(tag_double, &fy2Ptr[j], lanes);
      VectorDouble fz2 = highway::LoadN(tag_double, &fz2Ptr[j], lanes);

      auto newFx = fx2 - fxCombinedExt;
      auto newFy = fy2 - fyCombinedExt;
      auto newFz = fz2 - fzCombinedExt;

      highway::StoreN(newFx, tag_double, &fx2Ptr[j], lanes);
      highway::StoreN(newFy, tag_double, &fy2Ptr[j], lanes);
      highway::StoreN(newFz, tag_double, &fz2Ptr[j], lanes);
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      const auto lowerFx = highway::LowerHalf(fx);
      const auto lowerFy = highway::LowerHalf(fy);
      const auto lowerFz = highway::LowerHalf(fz);

      auto lowerFxExt = highway::ZeroExtendVector(tag_double, lowerFx);
      auto lowerFyExt = highway::ZeroExtendVector(tag_double, lowerFy);
      auto lowerFzExt = highway::ZeroExtendVector(tag_double, lowerFz);

      if constexpr (reversed) {
        if (j > i - _vecLengthDouble / 2) {
          const auto mask = overlapMasks[remainder ? rest : _vecLengthDouble / 2];
          lowerFxExt = highway::IfThenElseZero(mask, lowerFxExt);
          lowerFyExt = highway::IfThenElseZero(mask, lowerFyExt);
          lowerFzExt = highway::IfThenElseZero(mask, lowerFzExt);
        }
      }

      fx2Ptr[j] -= highway::ReduceSum(tag_double, lowerFxExt);
      fy2Ptr[j] -= highway::ReduceSum(tag_double, lowerFyExt);
      fz2Ptr[j] -= highway::ReduceSum(tag_double, lowerFzExt);

      if constexpr (!remainder) {
        const auto upperFx = highway::UpperHalf(tag_double_half, fx);
        const auto upperFy = highway::UpperHalf(tag_double_half, fy);
        const auto upperFz = highway::UpperHalf(tag_double_half, fz);

        auto upperFxExt = highway::ZeroExtendVector(tag_double, upperFx);
        auto upperFyExt = highway::ZeroExtendVector(tag_double, upperFy);
        auto upperFzExt = highway::ZeroExtendVector(tag_double, upperFz);

        if constexpr (reversed) {
          if (j >= i - _vecLengthDouble / 2) {
            const auto mask = overlapMasks[remainder ? rest : _vecLengthDouble / 2];
            upperFxExt = highway::IfThenElseZero(mask, upperFxExt);
            upperFyExt = highway::IfThenElseZero(mask, upperFyExt);
            upperFzExt = highway::IfThenElseZero(mask, upperFzExt);
          }
        }

        fx2Ptr[j + 1] -= highway::ReduceSum(tag_double, upperFxExt);
        fy2Ptr[j + 1] -= highway::ReduceSum(tag_double, upperFyExt);
        fz2Ptr[j + 1] -= highway::ReduceSum(tag_double, upperFzExt);
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      fx2Ptr[j] -= highway::ReduceSum(tag_double, fx);
      fy2Ptr[j] -= highway::ReduceSum(tag_double, fy);
      fz2Ptr[j] -= highway::ReduceSum(tag_double, fz);
    }
  }

  template <bool reversed, bool remainder, VectorizationPattern vecPattern>
  inline void reduceAccumulatedForce(const size_t i, double *const __restrict fxPtr, double *const __restrict fyPtr,
                                     double *const __restrict fzPtr, const VectorDouble &fxAcc,
                                     const VectorDouble &fyAcc, const VectorDouble &fzAcc, const int restI) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      fxPtr[i] += highway::ReduceSum(tag_double, fxAcc);
      fyPtr[i] += highway::ReduceSum(tag_double, fyAcc);
      fzPtr[i] += highway::ReduceSum(tag_double, fzAcc);
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
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
        const auto upperFxAcc = highway::UpperHalf(tag_double_half, fxAcc);
        const auto upperFyAcc = highway::UpperHalf(tag_double_half, fyAcc);
        const auto upperFzAcc = highway::UpperHalf(tag_double_half, fzAcc);

        const auto upperFxAccExt = highway::ZeroExtendVector(tag_double, upperFxAcc);
        const auto upperFyAccExt = highway::ZeroExtendVector(tag_double, upperFyAcc);
        const auto upperFzAccExt = highway::ZeroExtendVector(tag_double, upperFzAcc);

        int index = reversed ? i - 1 : i + 1;
        fxPtr[index] += highway::ReduceSum(tag_double, upperFxAccExt);
        fyPtr[index] += highway::ReduceSum(tag_double, upperFyAccExt);
        fzPtr[index] += highway::ReduceSum(tag_double, upperFzAccExt);
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      const auto lowerFxAcc = highway::LowerHalf(fxAcc);
      const auto lowerFyAcc = highway::LowerHalf(fyAcc);
      const auto lowerFzAcc = highway::LowerHalf(fzAcc);

      const auto upperFxAcc = highway::UpperHalf(tag_double_half, fxAcc);
      const auto upperFyAcc = highway::UpperHalf(tag_double_half, fyAcc);
      const auto upperFzAcc = highway::UpperHalf(tag_double_half, fzAcc);

      const auto fxAccCombined = lowerFxAcc + upperFxAcc;
      const auto fyAccCombined = lowerFyAcc + upperFyAcc;
      const auto fzAccCombined = lowerFzAcc + upperFzAcc;

      const auto fxAccExt = highway::ZeroExtendVector(tag_double, fxAccCombined);
      const auto fyAccExt = highway::ZeroExtendVector(tag_double, fyAccCombined);
      const auto fzAccExt = highway::ZeroExtendVector(tag_double, fzAccCombined);

      const int index = reversed ? (remainder ? 0 : i - _vecLengthDouble / 2 + 1) : i;

      const int lanes = remainder ? restI : _vecLengthDouble / 2;

      const VectorDouble oldFx = highway::LoadN(tag_double, &fxPtr[index], lanes);
      const VectorDouble oldFy = highway::LoadN(tag_double, &fyPtr[index], lanes);
      const VectorDouble oldFz = highway::LoadN(tag_double, &fzPtr[index], lanes);

      const auto newFx = oldFx + fxAccExt;
      const auto newFy = oldFy + fyAccExt;
      const auto newFz = oldFz + fzAccExt;

      highway::StoreN(newFx, tag_double, &fxPtr[index], lanes);
      highway::StoreN(newFy, tag_double, &fyPtr[index], lanes);
      highway::StoreN(newFz, tag_double, &fzPtr[index], lanes);
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      const VectorDouble oldFx =
          remainder ? highway::LoadN(tag_double, &fxPtr[i], restI) : highway::LoadU(tag_double, &fxPtr[i]);
      const VectorDouble oldFy =
          remainder ? highway::LoadN(tag_double, &fyPtr[i], restI) : highway::LoadU(tag_double, &fyPtr[i]);
      const VectorDouble oldFz =
          remainder ? highway::LoadN(tag_double, &fzPtr[i], restI) : highway::LoadU(tag_double, &fzPtr[i]);

      const VectorDouble fxNew = oldFx + fxAcc;
      const VectorDouble fyNew = oldFy + fyAcc;
      const VectorDouble fzNew = oldFz + fzAcc;

      remainder ? highway::StoreN(fxNew, tag_double, &fxPtr[i], restI) : highway::StoreU(fxNew, tag_double, &fxPtr[i]);
      remainder ? highway::StoreN(fyNew, tag_double, &fyPtr[i], restI) : highway::StoreU(fyNew, tag_double, &fyPtr[i]);
      remainder ? highway::StoreN(fzNew, tag_double, &fzPtr[i], restI) : highway::StoreU(fzNew, tag_double, &fzPtr[i]);
    }
  }

  template <bool newton3>
  inline void computeGlobals(const VectorDouble &virialSumX, const VectorDouble &virialSumY,
                             const VectorDouble &virialSumZ, const VectorDouble &uPotSum) {
    const int threadnum = autopas::autopas_get_thread_num();

    double globals[4]{highway::ReduceSum(tag_double, virialSumX), highway::ReduceSum(tag_double, virialSumY),
                      highway::ReduceSum(tag_double, virialSumZ), highway::ReduceSum(tag_double, uPotSum)};

    double factor = 1.;
    if constexpr (newton3) {
      factor = 0.5;
    }

    _aosThreadData[threadnum].virialSum[0] += globals[0] * factor;
    _aosThreadData[threadnum].virialSum[1] += globals[1] * factor;
    _aosThreadData[threadnum].virialSum[2] += globals[2] * factor;
    _aosThreadData[threadnum].uPotSum += globals[3] * factor;
  }

  template <bool reversed, bool newton3, bool remainderI, VectorizationPattern vecPattern>
  inline void handleILoopBody(
      const size_t i, const double *const __restrict xPtr1, const double *const __restrict yPtr1,
      const double *const __restrict zPtr1, const autopas::OwnershipState *const __restrict ownedStatePtr1,
      const double *const __restrict xPtr2, const double *const __restrict yPtr2, const double *const __restrict zPtr2,
      const autopas::OwnershipState *const __restrict ownedStatePtr2, double *const __restrict fxPtr1,
      double *const __restrict fyPtr1, double *const __restrict fzPtr1, double *const __restrict fxPtr2,
      double *const __restrict fyPtr2, double *const __restrict fzPtr2, const size_t *const __restrict typeIDptr1,
      const size_t *const __restrict typeIDptr2, VectorDouble &virialSumX, VectorDouble &virialSumY,
      VectorDouble &virialSumZ, VectorDouble &uPotSum, const size_t restI, const size_t jStop) {
    VectorDouble fxAcc = _zeroDouble;
    VectorDouble fyAcc = _zeroDouble;
    VectorDouble fzAcc = _zeroDouble;

    MaskDouble ownedMaskI;

    VectorDouble x1 = _zeroDouble;
    VectorDouble y1 = _zeroDouble;
    VectorDouble z1 = _zeroDouble;

    fillIRegisters<remainderI, reversed, vecPattern>(i, xPtr1, yPtr1, zPtr1, ownedStatePtr1, x1, y1, z1, ownedMaskI,
                                                     restI);

    unsigned int j = 0;
    for (; checkSecondLoopCondition<vecPattern>(jStop, j); incrementSecondLoop<vecPattern>(j)) {
      SoAKernel<newton3, remainderI, false, reversed, vecPattern>(
          i, j, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2, zPtr2, fxPtr2,
          fyPtr2, fzPtr2, &typeIDptr1[i], &typeIDptr2[j], fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ,
          uPotSum, restI, 0);
    }

    const int restJ = obtainSecondLoopRest<vecPattern>(jStop);
    if (restJ > 0) {
      SoAKernel<newton3, remainderI, true, reversed, vecPattern>(
          i, j, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2, zPtr2, fxPtr2,
          fyPtr2, fzPtr2, &typeIDptr1[i], &typeIDptr2[j], fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ,
          uPotSum, restI, restJ);
    }

    reduceAccumulatedForce<reversed, remainderI, vecPattern>(i, fxPtr1, fyPtr1, fzPtr1, fxAcc, fyAcc, fzAcc, restI);
  }

  /**
   * Templatized version of SoAFunctorSingle
   * @tparam newton3 Whether to use newton3 (TODO: compare with other functors)
   * @tparam vecPattern Vectorization Pattern for the interactions of the same list
   * @param soa
   */
  template <bool newton3, VectorizationPattern vecPattern>
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
    for (; checkFirstLoopCondition<true, vecPattern>(i, 0); decrementFirstLoop<vecPattern>(i)) {
      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");

      handleILoopBody<true, true, false, vecPattern>(i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr,
                                                     ownedStatePtr, fxPtr, fyPtr, fzPtr, fxPtr, fyPtr, fzPtr, typeIDptr,
                                                     typeIDptr, virialSumX, virialSumY, virialSumZ, uPotSum, 0, i);
    }

    if constexpr (vecPattern != VectorizationPattern::p1xVec) {
      // Rest I can't occur in 1xVec case
      const int restI = obtainFirstLoopRest<true>(i, -1);

      if (restI > 0) {
        handleILoopBody<true, true, true, vecPattern>(
            i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, fxPtr, fyPtr,
            fzPtr, typeIDptr, typeIDptr, virialSumX, virialSumY, virialSumZ, uPotSum, restI, i);
      }
    }

    if constexpr (calculateGlobals) {
      computeGlobals<newton3>(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  /**
   * Templatized version of SoAFunctorPairImpl
   * @tparam newton3
   * @tparam vecPattern
   * @param soa1
   * @param soa2
   */
  template <bool newton3, VectorizationPattern vecPattern>
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
    for (; checkFirstLoopCondition<false, vecPattern>(i, soa1.size()); incrementFirstLoop<vecPattern>(i)) {
      handleILoopBody<false, newton3, false, vecPattern>(
          i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx1Ptr, fy1Ptr, fz1Ptr, fx2Ptr,
          fy2Ptr, fz2Ptr, typeID1ptr, typeID2ptr, virialSumX, virialSumY, virialSumZ, uPotSum, 0, soa2.size());
    }

    if constexpr (vecPattern != VectorizationPattern::p1xVec) {
      // Rest I can't occur in 1xVec case
      const int restI = obtainFirstLoopRest<false>(i, soa1.size());
      if (restI > 0) {
        handleILoopBody<false, newton3, true, vecPattern>(
            i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx1Ptr, fy1Ptr, fz1Ptr, fx2Ptr,
            fy2Ptr, fz2Ptr, typeID1ptr, typeID2ptr, virialSumX, virialSumY, virialSumZ, uPotSum, restI, soa2.size());
      }
    }

    if constexpr (calculateGlobals) {
      computeGlobals<newton3>(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  template <bool remainder, VectorizationPattern vecPattern>
  inline void fillJRegisters(const size_t j, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
                             const double *const __restrict z2Ptr, const int64_t *const __restrict ownedStatePtr2,
                             VectorDouble &x2, VectorDouble &y2, VectorDouble &z2, VectorDouble &ownedStateJDouble,
                             const unsigned int rest) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      if constexpr (remainder) {
        x2 = highway::LoadN(tag_double, &x2Ptr[j], rest);
        y2 = highway::LoadN(tag_double, &y2Ptr[j], rest);
        z2 = highway::LoadN(tag_double, &z2Ptr[j], rest);

        const VectorLong ownedStateJ = highway::LoadN(tag_long, &ownedStatePtr2[j], rest);
        ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
      } else {
        x2 = highway::LoadU(tag_double, &x2Ptr[j]);
        y2 = highway::LoadU(tag_double, &y2Ptr[j]);
        z2 = highway::LoadU(tag_double, &z2Ptr[j]);

        const VectorLong ownedStateJ = highway::LoadU(tag_long, &ownedStatePtr2[j]);
        ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
      }
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      const int lanes = remainder ? rest : _vecLengthDouble / 2;

      VectorLong ownedStateJ = highway::LoadN(tag_long, &ownedStatePtr2[j], lanes);
      x2 = highway::LoadN(tag_double, &x2Ptr[j], lanes);
      y2 = highway::LoadN(tag_double, &y2Ptr[j], lanes);
      z2 = highway::LoadN(tag_double, &z2Ptr[j], lanes);

      // "broadcast" lower half to upper half
      ownedStateJ = highway::ConcatLowerLower(tag_long, ownedStateJ, ownedStateJ);
      ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
      x2 = highway::ConcatLowerLower(tag_double, x2, x2);
      y2 = highway::ConcatLowerLower(tag_double, y2, y2);
      z2 = highway::ConcatLowerLower(tag_double, z2, z2);
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      VectorLong ownedStateJ = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr2[j]));
      x2 = highway::Set(tag_double, x2Ptr[j]);
      y2 = highway::Set(tag_double, y2Ptr[j]);
      z2 = highway::Set(tag_double, z2Ptr[j]);

      if constexpr (remainder) {
        ownedStateJ = highway::ConcatLowerLower(tag_long, _zeroLong, ownedStateJ);
        x2 = highway::ConcatLowerLower(tag_double, _zeroDouble, x2);
        y2 = highway::ConcatLowerLower(tag_double, _zeroDouble, y2);
        z2 = highway::ConcatLowerLower(tag_double, _zeroDouble, z2);
      } else {
        const auto tmpOwnedJ = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr2[j + 1]));
        const auto tmpX2 = highway::Set(tag_double, x2Ptr[j + 1]);
        const auto tmpY2 = highway::Set(tag_double, y2Ptr[j + 1]);
        const auto tmpZ2 = highway::Set(tag_double, z2Ptr[j + 1]);

        ownedStateJ = highway::ConcatLowerLower(tag_long, tmpOwnedJ, ownedStateJ);
        x2 = highway::ConcatLowerLower(tag_double, tmpX2, x2);
        y2 = highway::ConcatLowerLower(tag_double, tmpY2, y2);
        z2 = highway::ConcatLowerLower(tag_double, tmpZ2, z2);
      }
      ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      VectorLong ownedStateJ = highway::Set(tag_long, static_cast<int64_t>(ownedStatePtr2[j]));
      x2 = highway::Set(tag_double, x2Ptr[j]);
      y2 = highway::Set(tag_double, y2Ptr[j]);
      z2 = highway::Set(tag_double, z2Ptr[j]);
      ownedStateJDouble = highway::ConvertTo(tag_double, ownedStateJ);
    }
  }

  template <bool remainderI, bool remainderJ, bool reversed, VectorizationPattern vecPattern>
  inline void fillPhysicsRegisters(const size_t *const typeID1Ptr, const size_t *const typeID2Ptr,
                                   VectorDouble &epsilon24s, VectorDouble &sigmaSquareds, VectorDouble &shift6s,
                                   const unsigned int j, const unsigned int restI, const unsigned int restJ) {
    HWY_ALIGN double epsilons[_vecLengthDouble] = {0.};
    HWY_ALIGN double sigmas[_vecLengthDouble] = {0.};
    HWY_ALIGN double shifts[_vecLengthDouble] = {0.};

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      for (int n = 0; n < (remainderJ ? restJ : _vecLengthDouble); ++n) {
        epsilons[n] = _PPLibrary->getMixing24Epsilon(*typeID1Ptr, *(typeID2Ptr + n));
        sigmas[n] = _PPLibrary->getMixingSigmaSquared(*typeID1Ptr, *(typeID2Ptr + n));
        if constexpr (applyShift) {
          shifts[n] = _PPLibrary->getMixingShift6(*typeID1Ptr, *(typeID2Ptr + n));
        }
      }
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      for (int i = 0; i < (remainderI ? 1 : 2); ++i) {
        for (int k = 0; k < (remainderJ ? restJ : _vecLengthDouble / 2); ++k) {
          int index = i * (_vecLengthDouble / 2) + k;
          auto typeID1 = reversed ? typeID1Ptr - i : typeID1Ptr + i;
          epsilons[index] = _PPLibrary->getMixing24Epsilon(*typeID1, *(typeID2Ptr + k));
          sigmas[index] = _PPLibrary->getMixingSigmaSquared(*typeID1, *(typeID2Ptr + k));

          if constexpr (applyShift) {
            shifts[index] = _PPLibrary->getMixingShift6(*typeID1, *(typeID2Ptr + k));
          }
        }
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      for (int i = 0; i < (remainderI ? restI : _vecLengthDouble / 2); ++i) {
        for (int k = 0; k < (remainderJ ? 1 : 2); ++k) {
          int index = i + _vecLengthDouble / 2 * k;
          auto typeID1 = reversed ? typeID1Ptr - i : typeID1Ptr + i;

          epsilons[index] = _PPLibrary->getMixing24Epsilon(*typeID1, *(typeID2Ptr + k));
          sigmas[index] = _PPLibrary->getMixingSigmaSquared(*typeID1, *(typeID2Ptr + k));

          if constexpr (applyShift) {
            shifts[index] = _PPLibrary->getMixingShift6(*typeID1, *(typeID2Ptr + k));
          }
        }
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      for (int n = 0; n < (remainderI ? restI : _vecLengthDouble); ++n) {
        auto typeID1 = reversed ? typeID1Ptr - n : typeID1Ptr + n;
        epsilons[n] = _PPLibrary->getMixing24Epsilon(*typeID1, *typeID2Ptr);
        sigmas[n] = _PPLibrary->getMixingSigmaSquared(*typeID1, *typeID2Ptr);

        if constexpr (applyShift) {
          shifts[n] = _PPLibrary->getMixingShift6(*typeID1, *typeID2Ptr);
        }
      }
    }

    epsilon24s = highway::Load(tag_double, epsilons);
    sigmaSquareds = highway::Load(tag_double, sigmas);
    if constexpr (applyShift) {
      shift6s = highway::Load(tag_double, shifts);
    }
  }

  /**
   * Actual innter kernel of the SoAFunctors
   *
   * @tparam newton3
   * @tparam remainderI
   * @tparam remainderJ
   * @tparam reversed
   * @tparam vecPattern
   * @param i
   * @param j
   * @param ownedMaskI
   * @param ownedStatePtr2
   * @param x1
   * @param y1
   * @param z1
   * @param x2Ptr
   * @param y2Ptr
   * @param z2Ptr
   * @param fx2Ptr
   * @param fy2Ptr
   * @param fz2Ptr
   * @param typeID1Ptr
   * @param typeID2Ptr
   * @param fxAcc
   * @param fyAcc
   * @param fzAcc
   * @param virialSumX
   * @param virialSumY
   * @param virialSumZ
   * @param uPotSum
   * @param restI
   * @param restJ
   */
  template <bool newton3, bool remainderI, bool remainderJ, bool reversed, VectorizationPattern vecPattern>
  inline void SoAKernel(const size_t i, const size_t j, const MaskDouble &ownedMaskI,
                        const int64_t *const __restrict ownedStatePtr2, const VectorDouble &x1, const VectorDouble &y1,
                        const VectorDouble &z1, const double *const __restrict x2Ptr,
                        const double *const __restrict y2Ptr, const double *const __restrict z2Ptr,
                        double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                        double *const __restrict fz2Ptr, const size_t *const typeID1Ptr, const size_t *const typeID2Ptr,
                        VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc, VectorDouble &virialSumX,
                        VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &uPotSum,
                        const unsigned int restI, const unsigned int restJ) {
    VectorDouble epsilon24s = _epsilon24;
    VectorDouble sigmaSquareds = _sigmaSquared;
    VectorDouble shift6s = _shift6;

    if constexpr (useMixing) {
      fillPhysicsRegisters<remainderI, remainderJ, reversed, vecPattern>(typeID1Ptr, typeID2Ptr, epsilon24s,
                                                                         sigmaSquareds, shift6s, j, restI, restJ);
    }

    VectorDouble x2;
    VectorDouble y2;
    VectorDouble z2;
    VectorDouble ownedStateJDouble;

    fillJRegisters<remainderJ, vecPattern>(j, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, x2, y2, z2, ownedStateJDouble,
                                           restJ);

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
    const auto cutoffDummyMask = highway::And(cutoffMask, dummyMask);

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
      handleNewton3Reduction<remainderJ, reversed, vecPattern>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, i, j, restJ);
    }

    if constexpr (calculateGlobals) {
      auto virialX = fx * drX;
      auto virialY = fy * drY;
      auto virialZ = fz * drZ;

      auto uPot = highway::MulAdd(epsilon24s, lj12m6, shift6s);
      auto uPotMasked = highway::IfThenElseZero(cutoffDummyMask, uPot);

      auto energyFactor = highway::IfThenElseZero(dummyMask, _oneDouble);

      if constexpr (newton3) {
        energyFactor = energyFactor + highway::IfThenElseZero(dummyMask, _oneDouble);
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

    HWY_ALIGN double x2Tmp[_vecLengthDouble] = {0.};
    HWY_ALIGN double y2Tmp[_vecLengthDouble] = {0.};
    HWY_ALIGN double z2Tmp[_vecLengthDouble] = {0.};
    HWY_ALIGN double fx2Tmp[_vecLengthDouble] = {0.};
    HWY_ALIGN double fy2Tmp[_vecLengthDouble] = {0.};
    HWY_ALIGN double fz2Tmp[_vecLengthDouble] = {0.};
    HWY_ALIGN size_t typeID2Tmp[_vecLengthDouble] = {0};
    HWY_ALIGN autopas::OwnershipState ownedStates2Tmp[_vecLengthDouble] = {autopas::OwnershipState::dummy};
    // alignas(64) std::array<double, _vecLengthDouble> x2Tmp{};
    // alignas(64) std::array<double, _vecLengthDouble> y2Tmp{};
    // alignas(64) std::array<double, _vecLengthDouble> z2Tmp{};
    // alignas(64) std::array<double, _vecLengthDouble> fx2Tmp{};
    // alignas(64) std::array<double, _vecLengthDouble> fy2Tmp{};
    // alignas(64) std::array<double, _vecLengthDouble> fz2Tmp{};
    // alignas(64) std::array<size_t, _vecLengthDouble> typeID2Tmp{};
    // alignas(64) std::array<autopas::OwnershipState, _vecLengthDouble> ownedStates2Tmp{};

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

      SoAKernel<newton3, false, false, false, VectorizationPattern::p1xVec>(
          0, 0, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStates2Tmp), x1, y1, z1, x2Tmp, y2Tmp, z2Tmp, fx2Tmp,
          fy2Tmp, fz2Tmp, &typeIDPtr[indexFirst], typeID2Tmp, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ,
          uPotSum, 0, 0);

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

      SoAKernel<newton3, false, true, false, VectorizationPattern::p1xVec>(
          0, 0, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStates2Tmp), x1, y1, z1, x2Tmp, y2Tmp, z2Tmp, fx2Tmp,
          fy2Tmp, fz2Tmp, &typeIDPtr[indexFirst], typeID2Tmp, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ,
          uPotSum, 0, rest);

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
      computeGlobals<newton3>(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
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
  /* TODO: compare with other functors */
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
  double getPotentialEnergy() {
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
    _epsilon24 = highway::Set(tag_double, epsilon24);
    _sigmaSquared = highway::Set(tag_double, sigmaSquare);
    if constexpr (applyShift) {
      _shift6 = highway::Set(tag_double, ParticlePropertiesLibrary<double, size_t>::calcShift6(
                                             epsilon24, sigmaSquare, highway::GetLane(_cutoffSquared)));
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

  /**
   * @copydoc autopas::Functor::setVecPattern()
   * Setter for the vectorization pattern to choose
   * @param vecPattern
   */
  void setVecPattern(const VectorizationPattern vecPattern) final { _vecPattern = vecPattern; }

 private:
  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, uPotSum{0.} {}

    void setZero() {
      virialSum = {0., 0., 0.};
      uPotSum = 0.;
    }

    std::array<double, 3> virialSum{};
    double uPotSum{0};

   private:
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

  // helper variables for the LJ-calculation
  const VectorDouble _zeroDouble{highway::Zero(tag_double)};
  const VectorLong _zeroLong{highway::Zero(tag_long)};
  const VectorDouble _oneDouble{highway::Set(tag_double, 1.)};
  const VectorLong _oneLong{highway::Set(tag_long, 1)};
  const VectorDouble _ownedStateDummy{highway::Zero(tag_double)};
  const VectorDouble _cutoffSquared{};
  VectorDouble _shift6{highway::Zero(tag_double)};
  VectorDouble _epsilon24{highway::Zero(tag_double)};
  VectorDouble _sigmaSquared{highway::Zero(tag_double)};

  MaskDouble equalMasks[_vecLengthDouble / 2];
  MaskDouble overlapMasks[_vecLengthDouble / 2];

  const double _cutoffSquareAoS{0.};
  double _epsilon24AoS, _sigmaSquareAoS, _shift6AoS = 0.;
  ParticlePropertiesLibrary<double, size_t> *_PPLibrary = nullptr;
  double _uPotSum{0.};
  std::array<double, 3> _virialSum;
  std::vector<AoSThreadData> _aosThreadData {};
  bool _postProcessed;
  bool _masksInitialized{false};

  VectorizationPattern _vecPattern;
};
}  // namespace mdLib