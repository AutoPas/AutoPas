/**
 * @file LJFunctorHWY.h
 *
 * @date 27.03.2024
 * @author Luis Gall
 */

#pragma once

#include <hwy/highway.h>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"

namespace mdLib {

namespace highway = hwy::HWY_NAMESPACE;
/** Highway tag for full double register */
constexpr highway::ScalableTag<double> tag_double;
/** Highway tag for full long register */
constexpr highway::ScalableTag<int64_t> tag_long;
/** Number of double values in a full register */
constexpr size_t _vecLengthDouble{highway::Lanes(tag_double)};
/** Type for a Double vector register */
using VectorDouble = decltype(highway::Zero(tag_double));
/** Type for a Long vector register */
using VectorLong = decltype(highway::Zero(tag_long));
/** Highway tag for a half-filled double register */
constexpr highway::Half<highway::DFromV<VectorDouble>> tag_double_half;
/** Type for a Double Mask */
using MaskDouble = decltype(highway::FirstN(tag_double, 1));
/** Type for a Long Mask */
using MaskLong = decltype(highway::FirstN(tag_long, 2));
/** Vectorization Pattern Type */
using VectorizationPattern = autopas::VectorizationPatternOption::Value;

/**
 * A functor to handle lennard-jones interactions between two particles (molecules)
 * This functor uses the SIMD abstraction library Google Highway to provide architecture independent vectorization
 * @tparam Particle_T The type of particle
 * @tparam applyShift Switch for the lj potential to be truncated shifted
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam relevantForTuning Whether the auto-tuner should consider this functor.
 * @tparam countFLOPs counts FLOPs and hitrate. Not implemented for this functor. Please use the AutoVec functor.*/
template <class Particle_T, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>

class LJFunctorHWY
    : public autopas::PairwiseFunctor<Particle_T, LJFunctorHWY<Particle_T, applyShift, useMixing, useNewton3,
                                                               calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = Particle_T::SoAArraysType;

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorHWY() = delete;

  /**
   * Constructor for Functor with mixing enabled/disabled. When using this functor without particlePropertiesLibrary it
   * is necessary to call setParticleProperties() to set internal constants.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit LJFunctorHWY(double cutoff, std::optional<std::reference_wrapper<ParticlePropertiesLibrary<double, size_t>>>
                                           particlePropertiesLibrary = std::nullopt)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorHWY>(cutoff),
        _cutoffSquared{highway::Set(tag_double, cutoff * cutoff)},
        _cutoffSquareAoS{cutoff * cutoff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _aosThreadData{},
        _postProcessed{false},
        _PPLibrary{particlePropertiesLibrary} {
    if (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }

    if constexpr (countFLOPs) {
      AutoPasLog(DEBUG, "Using LJFunctorHWY with countFLOPs but FLOP counting is not implemented.");
    }

    if constexpr (useMixing) {
      if (not _PPLibrary.has_value()) {
        throw std::runtime_error("Mixing is enabled but no ParticlePropertiesLibrary was provided!");
      }
    } else {
      if (_PPLibrary.has_value()) {
        throw std::runtime_error("Mixing is disabled but a ParticlePropertiesLibrary was provided!");
      }
    }
  }

  std::string getName() final { return "LJFunctorHWY"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  /**
   * Specifies whether the functor is capable of using the specified Vectorization Pattern in the SoA functor.
   * Thi functor can handle p1xVec, p2xVecDiv2, pVecDiv2x2, pVecx1
   *
   * @param vecPattern
   * @return whether the functor is capable of using the specified Vectorization Pattern
   */
  bool isVecPatternAllowed(const VectorizationPattern vecPattern) override final {
    return std::find(_vecPatternsAllowed.begin(), _vecPatternsAllowed.end(), vecPattern) != _vecPatternsAllowed.end();
  }

  /**
   * @copydoc autopas::PairwiseFunctor::AoSFunctor()
   */
  inline void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    auto sigmaSquare = _sigmaSquareAoS;
    auto epsilon24 = _epsilon24AoS;
    auto shift6 = _shift6AoS;
    if constexpr (useMixing) {
      sigmaSquare = _PPLibrary->get().getMixingSigmaSquared(i.getTypeId(), j.getTypeId());
      epsilon24 = _PPLibrary->get().getMixing24Epsilon(i.getTypeId(), j.getTypeId());
      if constexpr (applyShift) {
        shift6 = _PPLibrary->get().getMixingShift6(i.getTypeId(), j.getTypeId());
      }
    }
    const auto dr = i.getR() - j.getR();
    const double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffSquareAoS) {
      return;
    }

    const double invdr2 = 1. / dr2;
    double lj6 = sigmaSquare * invdr2;
    lj6 = lj6 * lj6 * lj6;
    const double lj12 = lj6 * lj6;
    const double lj12m6 = lj12 - lj6;
    const double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    const auto f = dr * fac;
    i.addF(f);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(f);
    }
    if (calculateGlobals) {
      const auto virial = dr * f;
      const double potentialEnergy6 = epsilon24 * lj12m6 + shift6;

      const int threadnum = autopas::autopas_get_thread_num();

      if (i.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadData[threadnum].virialSum += virial;
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        // for non-newton3 the division is in the post-processing step.
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy6;
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
        autopas::utils::ExceptionHandler::exception("Unknown VectorizationPattern!");
    }
  }

 private:
  /**
   * Checks whether the loop over i should be continued. This depends on the respective VectorizationPattern.
   *
   * @tparam reversed Whether iterating backwards over i.
   * @tparam vecPattern
   * @param i
   * @param vecEnd
   * @return true if the loop over i should be continued, else false.
   */
  template <bool reversed, VectorizationPattern vecPattern>
  constexpr bool checkFirstLoopCondition(const std::ptrdiff_t i, const long vecEnd) {
    if constexpr (reversed) {
      if constexpr (vecPattern == VectorizationPattern::p1xVec) {
        return i >= 0;
      } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
        return i >= 1;
      } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
        return i >= _vecLengthDouble / 2;
      } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
        return i >= _vecLengthDouble;
      } else {
        return false;
      }
    } else {
      if constexpr (vecPattern == VectorizationPattern::p1xVec) {
        return i < vecEnd;
      } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
        return i < (vecEnd - 1);
      } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
        // Switching to signed long to avoid integer underflows.
        // Note: The narrowing conversion of _vecLengthDouble shouldn't cause issues as it is not expected to be so
        // large.
        const long criterion = vecEnd - _vecLengthDouble / 2 + 1;
        return i < criterion;
      } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
        // Switching to signed long to avoid integer underflows.
        // Note: The narrowing conversion of _vecLengthDouble shouldn't cause issues as it is not expected to be so
        // large.
        const long criterion = vecEnd - _vecLengthDouble + 1;
        return i < criterion;
      } else {
        return false;
      }
    }
  }

  /**
   * Calculates the step size by which i is decremented in the outer loop, depending on VectorizationPattern. Should
   * only be called when iterating backwards over i.
   *
   * @tparam vecPattern
   * @param i
   */
  template <VectorizationPattern vecPattern>
  constexpr void decrementFirstLoop(std::ptrdiff_t &i) {
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

  /**
   * Calculates the step size by which i is incremented in the outer loop, depending on VectorizationPattern.
   * Should only be called when iterating forwards over i.
   *
   * @tparam vecPattern
   * @param i
   */
  template <VectorizationPattern vecPattern>
  constexpr void incrementFirstLoop(std::ptrdiff_t &i) {
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

  /**
   * Determines the number of registers filled in the case of a final, only partially filled, kernel call in the loop
   * over i.
   *
   * @tparam reversed
   * @param i
   * @param vecEnd
   * @return the number of lanes filled.
   */
  template <bool reversed>
  constexpr int obtainILoopRemainderLength(std::ptrdiff_t i, const int vecEnd) {
    return reversed ? (i < 0 ? 0 : i + 1) : vecEnd - i;
  }

  /**
   * Checks whether the inner loop over j should be continued. This depends on the respective VectorizationPattern.
   *
   * @tparam vecPattern
   * @param i
   * @param j
   * @return  true if the inner loop over j should be continued, else false.
   */
  template <VectorizationPattern vecPattern>
  constexpr bool checkSecondLoopCondition(std::ptrdiff_t i, size_t j) {
    std::ptrdiff_t limit = 0;

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      // Round i down to the next multiple of _vecLengthDouble
      limit = i - (i % _vecLengthDouble);
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      // Round i down to the next multiple of _vecLengthDouble / 2
      const std::ptrdiff_t block = _vecLengthDouble / 2;
      limit = i - (i % block);
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      // Round i down to the next multiple of 2
      limit = i - (i % 2);
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      // Rounding to a multiple of 1 is a no-op
      limit = i;
    } else {
      // Unknown vectorization pattern
      return false;
    }

    return j < static_cast<size_t>(limit);
  }

  /**
   * Calculates the step size by which j is incremented in the inner loop, depending on VectorizationPattern.
   * Should only be called when iterating forwards over i.
   *
   * @tparam vecPattern
   * @param j
   */
  template <VectorizationPattern vecPattern>
  constexpr void incrementSecondLoop(std::ptrdiff_t &j) {
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

  /**
   * Determines the number of registers filled in the case of a final, only partially filled, kernel call in the loop
   * over j.
   *
   * @tparam vecPattern
   * @param j
   * @return the number of lanes filled.
   */
  template <VectorizationPattern vecPattern>
  constexpr int obtainJLoopRemainderLength(const std::ptrdiff_t j) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      return static_cast<int>(j & (_vecLengthDouble - 1));
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      return static_cast<int>(j & (_vecLengthDouble / 2 - 1));
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      return static_cast<int>(j & (1));
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      return 0;
    } else {
      return -1;
    }
  }

  /**
   * Depending on the vectorization pattern, this function fills the registers for the loop over i. This includes the
   * positions and the ownership state.
   *
   * @tparam remainder
   * @tparam reversed
   * @tparam vecPattern
   * @param i Current loop index for i.
   * @param xPtr pointer to the x coordinate of the particle at position i.
   * @param yPtr pointer to the y coordinate of the particle at position i.
   * @param zPtr pointer to the z coordinate of the particle at position i.
   * @param ownedStatePtr pointer to the ownership state of the particle at position i.
   * @param x1 The register to be filled for x coordinates.
   * @param y1 The register to be filled for y coordinates.
   * @param z1 The register to be filled for z coordinates.
   * @param ownedMaskI The mask to be filled for the ownership state.
   * @param restI The number of lanes filled in case of a remainder loop.
   */
  template <bool remainder, bool reversed, VectorizationPattern vecPattern>
  inline void fillIRegisters(const size_t i, const double *const __restrict xPtr, const double *const __restrict yPtr,
                             const double *const __restrict zPtr,
                             const autopas::OwnershipState *const __restrict ownedStatePtr, VectorDouble &x1,
                             VectorDouble &y1, VectorDouble &z1, MaskDouble &ownedMaskI, const size_t restI) {
    VectorLong ownedStateILong = _zeroLong;

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      const auto owned = static_cast<int64_t>(ownedStatePtr[i]);
      ownedStateILong = highway::Set(tag_long, owned);

      x1 = highway::Set(tag_double, xPtr[i]);
      y1 = highway::Set(tag_double, yPtr[i]);
      z1 = highway::Set(tag_double, zPtr[i]);
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      const auto ownedFirst = static_cast<int64_t>(ownedStatePtr[i]);
      ownedStateILong = highway::Set(tag_long, ownedFirst);

      x1 = highway::Set(tag_double, xPtr[i]);
      y1 = highway::Set(tag_double, yPtr[i]);
      z1 = highway::Set(tag_double, zPtr[i]);

      VectorLong tmpOwnedI = _zeroLong;
      VectorDouble tmpX1 = _zeroDouble;
      VectorDouble tmpY1 = _zeroDouble;
      VectorDouble tmpZ1 = _zeroDouble;

      if constexpr (not remainder) {
        const auto index = reversed ? i - 1 : i + 1;
        const auto ownedSecond = static_cast<int64_t>(ownedStatePtr[index]);
        tmpOwnedI = highway::Set(tag_long, ownedSecond);
        tmpX1 = highway::Set(tag_double, xPtr[index]);
        tmpY1 = highway::Set(tag_double, yPtr[index]);
        tmpZ1 = highway::Set(tag_double, zPtr[index]);
      }

      ownedStateILong = highway::ConcatLowerLower(tag_long, tmpOwnedI, ownedStateILong);
      x1 = highway::ConcatLowerLower(tag_double, tmpX1, x1);
      y1 = highway::ConcatLowerLower(tag_double, tmpY1, y1);
      z1 = highway::ConcatLowerLower(tag_double, tmpZ1, z1);
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      const int index = reversed ? (remainder ? 0 : i - _vecLengthDouble / 2 + 1) : i;
      const int lanes = remainder ? restI : _vecLengthDouble / 2;

      ownedStateILong = highway::LoadN(tag_long, reinterpret_cast<const int64_t *>(&ownedStatePtr[index]), lanes);

      x1 = highway::LoadN(tag_double, &xPtr[index], lanes);
      y1 = highway::LoadN(tag_double, &yPtr[index], lanes);
      z1 = highway::LoadN(tag_double, &zPtr[index], lanes);

      ownedStateILong = highway::ConcatLowerLower(tag_long, ownedStateILong, ownedStateILong);
      x1 = highway::ConcatLowerLower(tag_double, x1, x1);
      y1 = highway::ConcatLowerLower(tag_double, y1, y1);
      z1 = highway::ConcatLowerLower(tag_double, z1, z1);
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      const auto index = reversed ? (remainder ? 0 : i - _vecLengthDouble + 1) : i;

      if constexpr (remainder) {
        x1 = highway::LoadN(tag_double, &xPtr[index], restI);
        y1 = highway::LoadN(tag_double, &yPtr[index], restI);
        z1 = highway::LoadN(tag_double, &zPtr[index], restI);

        ownedStateILong = highway::LoadN(tag_long, reinterpret_cast<const int64_t *>(&ownedStatePtr[index]), restI);
      } else {
        x1 = highway::LoadU(tag_double, &xPtr[index]);
        y1 = highway::LoadU(tag_double, &yPtr[index]);
        z1 = highway::LoadU(tag_double, &zPtr[index]);

        ownedStateILong = highway::LoadU(tag_long, reinterpret_cast<const int64_t *>(&ownedStatePtr[index]));
      }
    }

    MaskLong ownedMaskILong = highway::Ne(ownedStateILong, _zeroLong);

    // conert to a double mask since we perform logical operations with other double masks in the kernel.
    ownedMaskI = highway::RebindMask(tag_double, ownedMaskILong);
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
      const auto lowerFx = highway::LowerHalf(tag_double_half, fx);
      const auto lowerFy = highway::LowerHalf(tag_double_half, fy);
      const auto lowerFz = highway::LowerHalf(tag_double_half, fz);

      const auto upperFx = highway::UpperHalf(tag_double_half, fx);
      const auto upperFy = highway::UpperHalf(tag_double_half, fy);
      const auto upperFz = highway::UpperHalf(tag_double_half, fz);

      const auto fxCombined = lowerFx + upperFx;
      const auto fyCombined = lowerFy + upperFy;
      const auto fzCombined = lowerFz + upperFz;

      const int lanes = remainder ? rest : _vecLengthDouble / 2;

      const auto fx2 = highway::LoadN(tag_double_half, &fx2Ptr[j], lanes);
      const auto fy2 = highway::LoadN(tag_double_half, &fy2Ptr[j], lanes);
      const auto fz2 = highway::LoadN(tag_double_half, &fz2Ptr[j], lanes);

      const auto newFx = fx2 - fxCombined;
      const auto newFy = fy2 - fyCombined;
      const auto newFz = fz2 - fzCombined;

      highway::StoreN(newFx, tag_double_half, &fx2Ptr[j], lanes);
      highway::StoreN(newFy, tag_double_half, &fy2Ptr[j], lanes);
      highway::StoreN(newFz, tag_double_half, &fz2Ptr[j], lanes);
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      const auto lowerFx = highway::LowerHalf(tag_double_half, fx);
      const auto lowerFy = highway::LowerHalf(tag_double_half, fy);
      const auto lowerFz = highway::LowerHalf(tag_double_half, fz);

      fx2Ptr[j] -= highway::ReduceSum(tag_double_half, lowerFx);
      fy2Ptr[j] -= highway::ReduceSum(tag_double_half, lowerFy);
      fz2Ptr[j] -= highway::ReduceSum(tag_double_half, lowerFz);

      if constexpr (not remainder) {
        const auto upperFx = highway::UpperHalf(tag_double_half, fx);
        const auto upperFy = highway::UpperHalf(tag_double_half, fy);
        const auto upperFz = highway::UpperHalf(tag_double_half, fz);

        fx2Ptr[j + 1] -= highway::ReduceSum(tag_double_half, upperFx);
        fy2Ptr[j + 1] -= highway::ReduceSum(tag_double_half, upperFy);
        fz2Ptr[j + 1] -= highway::ReduceSum(tag_double_half, upperFz);
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
      const auto lowerFxAcc = highway::LowerHalf(tag_double_half, fxAcc);
      const auto lowerFyAcc = highway::LowerHalf(tag_double_half, fyAcc);
      const auto lowerFzAcc = highway::LowerHalf(tag_double_half, fzAcc);

      fxPtr[i] += highway::ReduceSum(tag_double_half, lowerFxAcc);
      fyPtr[i] += highway::ReduceSum(tag_double_half, lowerFyAcc);
      fzPtr[i] += highway::ReduceSum(tag_double_half, lowerFzAcc);

      if constexpr (not remainder) {
        const auto upperFxAcc = highway::UpperHalf(tag_double_half, fxAcc);
        const auto upperFyAcc = highway::UpperHalf(tag_double_half, fyAcc);
        const auto upperFzAcc = highway::UpperHalf(tag_double_half, fzAcc);

        const auto index = reversed ? i - 1 : i + 1;
        fxPtr[index] += highway::ReduceSum(tag_double_half, upperFxAcc);
        fyPtr[index] += highway::ReduceSum(tag_double_half, upperFyAcc);
        fzPtr[index] += highway::ReduceSum(tag_double_half, upperFzAcc);
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      const auto lowerFxAcc = highway::LowerHalf(tag_double_half, fxAcc);
      const auto lowerFyAcc = highway::LowerHalf(tag_double_half, fyAcc);
      const auto lowerFzAcc = highway::LowerHalf(tag_double_half, fzAcc);

      const auto upperFxAcc = highway::UpperHalf(tag_double_half, fxAcc);
      const auto upperFyAcc = highway::UpperHalf(tag_double_half, fyAcc);
      const auto upperFzAcc = highway::UpperHalf(tag_double_half, fzAcc);

      const auto fxAccCombined = lowerFxAcc + upperFxAcc;
      const auto fyAccCombined = lowerFyAcc + upperFyAcc;
      const auto fzAccCombined = lowerFzAcc + upperFzAcc;

      const int index = reversed ? (remainder ? 0 : i - _vecLengthDouble / 2 + 1) : i;

      const int lanes = remainder ? restI : _vecLengthDouble / 2;

      const auto oldFx = highway::LoadN(tag_double_half, &fxPtr[index], lanes);
      const auto oldFy = highway::LoadN(tag_double_half, &fyPtr[index], lanes);
      const auto oldFz = highway::LoadN(tag_double_half, &fzPtr[index], lanes);

      const auto newFx = oldFx + fxAccCombined;
      const auto newFy = oldFy + fyAccCombined;
      const auto newFz = oldFz + fzAccCombined;

      highway::StoreN(newFx, tag_double_half, &fxPtr[index], lanes);
      highway::StoreN(newFy, tag_double_half, &fyPtr[index], lanes);
      highway::StoreN(newFz, tag_double_half, &fzPtr[index], lanes);
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

    const std::array<double, 4> globals{
        highway::ReduceSum(tag_double, virialSumX), highway::ReduceSum(tag_double, virialSumY),
        highway::ReduceSum(tag_double, virialSumZ), highway::ReduceSum(tag_double, uPotSum)};

    _aosThreadData[threadnum].virialSum[0] += globals[0];
    _aosThreadData[threadnum].virialSum[1] += globals[1];
    _aosThreadData[threadnum].virialSum[2] += globals[2];
    _aosThreadData[threadnum].potentialEnergySum += globals[3];
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
      VectorDouble &virialSumZ, VectorDouble &uPotSum, const size_t restI, const size_t jVecEnd) {
    VectorDouble fxAcc = _zeroDouble;
    VectorDouble fyAcc = _zeroDouble;
    VectorDouble fzAcc = _zeroDouble;

    MaskDouble ownedMaskI;

    VectorDouble x1 = _zeroDouble;
    VectorDouble y1 = _zeroDouble;
    VectorDouble z1 = _zeroDouble;

    fillIRegisters<remainderI, reversed, vecPattern>(i, xPtr1, yPtr1, zPtr1, ownedStatePtr1, x1, y1, z1, ownedMaskI,
                                                     restI);

    std::ptrdiff_t j = 0;
    for (; checkSecondLoopCondition<vecPattern>(jVecEnd, j); incrementSecondLoop<vecPattern>(j)) {
      SoAKernel<newton3, remainderI, false, reversed, vecPattern>(
          i, j, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2, zPtr2, fxPtr2,
          fyPtr2, fzPtr2, &typeIDptr1[i], &typeIDptr2[j], fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ,
          uPotSum, restI, 0);
    }

    const int restJ = obtainJLoopRemainderLength<vecPattern>(jVecEnd);
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
   * @tparam newton3 Whether to use newton3
   * @param soa
   */
  template <bool newton3>
  inline void SoAFunctorSingleImpl(autopas::SoAView<SoAArraysType> soa) {
    if (soa.size() == 0) return;

    // obtain iterators for the various values
    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    const auto *const __restrict typeIDptr = soa.template begin<Particle_T::AttributeNames::typeId>();

    // initialize and declare vector variables
    auto virialSumX = _zeroDouble;
    auto virialSumY = _zeroDouble;
    auto virialSumZ = _zeroDouble;
    auto uPotSum = _zeroDouble;

    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(soa.size()) - 1;
         checkFirstLoopCondition<true, VectorizationPattern::p1xVec>(i, 0);
         decrementFirstLoop<VectorizationPattern::p1xVec>(i)) {
      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");

      handleILoopBody<true, true, false, VectorizationPattern::p1xVec>(
          i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, fxPtr, fyPtr, fzPtr,
          typeIDptr, typeIDptr, virialSumX, virialSumY, virialSumZ, uPotSum, 0, i);
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

    const auto *const __restrict x1Ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y1Ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z1Ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict x2Ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y2Ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z2Ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fx1Ptr = soa1.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy1Ptr = soa1.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz1Ptr = soa1.template begin<Particle_T::AttributeNames::forceZ>();
    auto *const __restrict fx2Ptr = soa2.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy2Ptr = soa2.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz2Ptr = soa2.template begin<Particle_T::AttributeNames::forceZ>();

    const auto *const __restrict typeID1ptr = soa1.template begin<Particle_T::AttributeNames::typeId>();
    const auto *const __restrict typeID2ptr = soa2.template begin<Particle_T::AttributeNames::typeId>();

    VectorDouble virialSumX = _zeroDouble;
    VectorDouble virialSumY = _zeroDouble;
    VectorDouble virialSumZ = _zeroDouble;
    VectorDouble uPotSum = _zeroDouble;

    std::ptrdiff_t i = 0;
    for (; checkFirstLoopCondition<false, vecPattern>(i, soa1.size()); incrementFirstLoop<vecPattern>(i)) {
      handleILoopBody<false, newton3, false, vecPattern>(
          i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx1Ptr, fy1Ptr, fz1Ptr, fx2Ptr,
          fy2Ptr, fz2Ptr, typeID1ptr, typeID2ptr, virialSumX, virialSumY, virialSumZ, uPotSum, 0, soa2.size());
    }

    if constexpr (vecPattern != VectorizationPattern::p1xVec) {
      // Rest I can't occur in 1xVec case
      const int restI = obtainILoopRemainderLength<false>(i, soa1.size());
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

  /**
   * Depending on the vectorization pattern, this function fills the registers for the loop over j. This includes the
   * positions and the ownership state.
   *
   * @tparam remainder
   * @tparam vecPattern
   * @param j Current loop index for j.
   * @param x2Ptr pointer to the x coordinate of the particle at position j.
   * @param y2Ptr pointer to the y coordinate of the particle at position j.
   * @param z2Ptr pointer to the z coordinate of the particle at position j.
   * @param ownedStatePtr2 pointer to the ownership state of the particle at position j.
   * @param x2 The register to be filled for x coordinates.
   * @param y2 The register to be filled for y coordinates.
   * @param z2 The register to be filled for z coordinates.
   * @param ownedMaskJ The register to be filled for the ownership state.
   * @param rest The number of lanes filled in case of a remainder loop.
   */
  template <bool remainder, VectorizationPattern vecPattern>
  inline void fillJRegisters(const size_t j, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
                             const double *const __restrict z2Ptr, const int64_t *const __restrict ownedStatePtr2,
                             VectorDouble &x2, VectorDouble &y2, VectorDouble &z2, MaskDouble &ownedMaskJ,
                             const unsigned int rest) {
    VectorLong ownedStateJLong = _zeroLong;

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      if constexpr (remainder) {
        x2 = highway::LoadN(tag_double, &x2Ptr[j], rest);
        y2 = highway::LoadN(tag_double, &y2Ptr[j], rest);
        z2 = highway::LoadN(tag_double, &z2Ptr[j], rest);

        ownedStateJLong = highway::LoadN(tag_long, &ownedStatePtr2[j], rest);
      } else {
        x2 = highway::LoadU(tag_double, &x2Ptr[j]);
        y2 = highway::LoadU(tag_double, &y2Ptr[j]);
        z2 = highway::LoadU(tag_double, &z2Ptr[j]);

        ownedStateJLong = highway::LoadU(tag_long, &ownedStatePtr2[j]);
      }
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      const int lanes = remainder ? rest : _vecLengthDouble / 2;

      VectorLong ownedStateJ = highway::LoadN(tag_long, &ownedStatePtr2[j], lanes);
      x2 = highway::LoadN(tag_double, &x2Ptr[j], lanes);
      y2 = highway::LoadN(tag_double, &y2Ptr[j], lanes);
      z2 = highway::LoadN(tag_double, &z2Ptr[j], lanes);

      // "broadcast" lower half to upper half
      ownedStateJLong = highway::ConcatLowerLower(tag_long, ownedStateJ, ownedStateJ);
      x2 = highway::ConcatLowerLower(tag_double, x2, x2);
      y2 = highway::ConcatLowerLower(tag_double, y2, y2);
      z2 = highway::ConcatLowerLower(tag_double, z2, z2);
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      VectorLong ownedStateJ = highway::Set(tag_long, ownedStatePtr2[j]);
      x2 = highway::Set(tag_double, x2Ptr[j]);
      y2 = highway::Set(tag_double, y2Ptr[j]);
      z2 = highway::Set(tag_double, z2Ptr[j]);

      if constexpr (remainder) {
        ownedStateJLong = highway::ConcatLowerLower(tag_long, _zeroLong, ownedStateJ);
        x2 = highway::ConcatLowerLower(tag_double, _zeroDouble, x2);
        y2 = highway::ConcatLowerLower(tag_double, _zeroDouble, y2);
        z2 = highway::ConcatLowerLower(tag_double, _zeroDouble, z2);
      } else {
        const auto tmpOwnedJ = highway::Set(tag_long, ownedStatePtr2[j + 1]);
        const auto tmpX2 = highway::Set(tag_double, x2Ptr[j + 1]);
        const auto tmpY2 = highway::Set(tag_double, y2Ptr[j + 1]);
        const auto tmpZ2 = highway::Set(tag_double, z2Ptr[j + 1]);

        ownedStateJLong = highway::ConcatLowerLower(tag_long, tmpOwnedJ, ownedStateJ);
        x2 = highway::ConcatLowerLower(tag_double, tmpX2, x2);
        y2 = highway::ConcatLowerLower(tag_double, tmpY2, y2);
        z2 = highway::ConcatLowerLower(tag_double, tmpZ2, z2);
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      ownedStateJLong = highway::Set(tag_long, ownedStatePtr2[j]);
      x2 = highway::Set(tag_double, x2Ptr[j]);
      y2 = highway::Set(tag_double, y2Ptr[j]);
      z2 = highway::Set(tag_double, z2Ptr[j]);
    }

    MaskLong ownedMaskJLong = highway::Ne(ownedStateJLong, _zeroLong);

    // convert to a double mask since we perform logical operations with other double masks in the kernel.
    ownedMaskJ = highway::RebindMask(tag_double, ownedMaskJLong);
  }

  template <bool remainderI, bool remainderJ, bool reversed, VectorizationPattern vecPattern>
  inline void fillPhysicsRegisters(const size_t *const typeID1Ptr, const size_t *const typeID2Ptr,
                                   VectorDouble &epsilon24s, VectorDouble &sigmaSquareds, VectorDouble &shift6s,
                                   const unsigned int j, const unsigned int restI, const unsigned int restJ) {
    HWY_ALIGN double epsilons[_vecLengthDouble] = {0.};
    HWY_ALIGN double sigmas[_vecLengthDouble] = {0.};
    HWY_ALIGN double shifts[_vecLengthDouble] = {0.};

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      for (int j = 0; j < (remainderJ ? restJ : _vecLengthDouble); ++j) {
        epsilons[j] = _PPLibrary->get().getMixing24Epsilon(*typeID1Ptr, *(typeID2Ptr + j));
        sigmas[j] = _PPLibrary->get().getMixingSigmaSquared(*typeID1Ptr, *(typeID2Ptr + j));
        if constexpr (applyShift) {
          shifts[j] = _PPLibrary->get().getMixingShift6(*typeID1Ptr, *(typeID2Ptr + j));
        }
      }
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      for (int i = 0; i < (remainderI ? 1 : 2); ++i) {
        for (int j = 0; j < (remainderJ ? restJ : _vecLengthDouble / 2); ++j) {
          const auto index = i * (_vecLengthDouble / 2) + j;
          const auto typeID1 = reversed ? typeID1Ptr - i : typeID1Ptr + i;
          epsilons[index] = _PPLibrary->get().getMixing24Epsilon(*typeID1, *(typeID2Ptr + j));
          sigmas[index] = _PPLibrary->get().getMixingSigmaSquared(*typeID1, *(typeID2Ptr + j));

          if constexpr (applyShift) {
            shifts[index] = _PPLibrary->get().getMixingShift6(*typeID1, *(typeID2Ptr + j));
          }
        }
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      for (int i = 0; i < (remainderI ? restI : _vecLengthDouble / 2); ++i) {
        for (int j = 0; j < (remainderJ ? 1 : 2); ++j) {
          const auto index = i + _vecLengthDouble / 2 * j;
          const auto typeID1 = reversed ? typeID1Ptr - i : typeID1Ptr + i;

          epsilons[index] = _PPLibrary->get().getMixing24Epsilon(*typeID1, *(typeID2Ptr + j));
          sigmas[index] = _PPLibrary->get().getMixingSigmaSquared(*typeID1, *(typeID2Ptr + j));

          if constexpr (applyShift) {
            shifts[index] = _PPLibrary->get().getMixingShift6(*typeID1, *(typeID2Ptr + j));
          }
        }
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      for (int i = 0; i < (remainderI ? restI : _vecLengthDouble); ++i) {
        auto typeID1 = reversed ? typeID1Ptr - i : typeID1Ptr + i;
        epsilons[i] = _PPLibrary->get().getMixing24Epsilon(*typeID1, *typeID2Ptr);
        sigmas[i] = _PPLibrary->get().getMixingSigmaSquared(*typeID1, *typeID2Ptr);

        if constexpr (applyShift) {
          shifts[i] = _PPLibrary->get().getMixingShift6(*typeID1, *typeID2Ptr);
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
   * Actual inner kernel of the SoAFunctors
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
    MaskDouble ownedMaskJ;

    fillJRegisters<remainderJ, vecPattern>(j, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, x2, y2, z2, ownedMaskJ, restJ);

    // distance calculations
    const auto drX = x1 - x2;
    const auto drY = y1 - y2;
    const auto drZ = z1 - z2;

    const auto drX2 = drX * drX;
    const auto drY2 = drY * drY;
    const auto drZ2 = drZ * drZ;

    const auto dr2 = drX2 + drY2 + drZ2;

    const auto dummyMask = highway::And(ownedMaskI, ownedMaskJ);
    const auto cutoffDummyMask = highway::MaskedLe(dummyMask, dr2, _cutoffSquared);

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
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
      return;
    }

    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    const auto *const __restrict typeIDPtr = soa.template begin<Particle_T::AttributeNames::typeId>();

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

    size_t j = 0;
    const size_t vecEnd = (neighborList.size() / _vecLengthDouble) * _vecLengthDouble;

    for (; j < vecEnd; j += _vecLengthDouble) {
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
    return std::array<typename Particle_T::AttributeNames, 9>{Particle_T::AttributeNames::id,
                                                              Particle_T::AttributeNames::posX,
                                                              Particle_T::AttributeNames::posY,
                                                              Particle_T::AttributeNames::posZ,
                                                              Particle_T::AttributeNames::forceX,
                                                              Particle_T::AttributeNames::forceY,
                                                              Particle_T::AttributeNames::forceZ,
                                                              Particle_T::AttributeNames::typeId,
                                                              Particle_T::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle_T::AttributeNames, 6>{
        Particle_T::AttributeNames::id,     Particle_T::AttributeNames::posX,
        Particle_T::AttributeNames::posY,   Particle_T::AttributeNames::posZ,
        Particle_T::AttributeNames::typeId, Particle_T::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 3>{
        Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY, Particle_T::AttributeNames::forceZ};
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
   * Accumulates global values, e.g. upot and virial.
   * @param newton3
   */
  void endTraversal(const bool newton3) final {
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
      _potentialEnergySum *= 0.5;
      _virialSum *= 0.5;

      // We have always calculated 6*potentialEnergy, so we divide by 6 here!
      _potentialEnergySum /= 6.;
      _postProcessed = true;

      AutoPasLog(DEBUG, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(DEBUG, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
    }
  }

  /**
   * Get the potential Energy
   * @return the potential Energy
   */
  double getPotentialEnergy() const {
    if (not calculateGlobals) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Trying to get upot even though calculateGlobals is false. If you want this functor to calculate global "
          "values, please specify calculateGlobals to be true.");
    }
    if (not _postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException("Cannot get upot, because endTraversal was not called.");
    }
    return _potentialEnergySum;
  }

  /**
   * Get the virial
   * @return the virial
   */
  double getVirial() const {
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
  void setParticleProperties(const double epsilon24, const double sigmaSquare) {
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
   */
  void setVecPattern(const VectorizationPattern vecPattern) final { _vecPattern = vecPattern; }

 private:
  /**
   * This class stores internal data of each thread, make sure that this data has proper size, i.e. k*64 Bytes!
   */
  class AoSThreadData {
   public:
    AoSThreadData() : virialSum{0., 0., 0.}, potentialEnergySum{0.} {}

    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }

    std::array<double, 3> virialSum{};
    double potentialEnergySum{};

   private:
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };
  static_assert(sizeof(AoSThreadData) % 64 == 0, "AoSThreadData has wrong size");

  // helper variables for the LJ-calculation used in the kernel.
  // vector register of doubles containing only zeros.
  const VectorDouble _zeroDouble{highway::Zero(tag_double)};
  // vector register of long integers containing only zeros.
  const VectorLong _zeroLong{highway::Zero(tag_long)};
  // vector register of doubles containing only ones.
  const VectorDouble _oneDouble{highway::Set(tag_double, 1.)};
  // vector register of long integers containing only ones.
  const VectorLong _oneLong{highway::Set(tag_long, 1)};
  // vector register of doubles containing dummy ownership state.
  const VectorDouble _ownedStateDummy{highway::Zero(tag_double)};
  // vector register to hold the squared cutoff in all lanes.
  const VectorDouble _cutoffSquared{};
  // vector register to hold the _hift6 values.
  VectorDouble _shift6{highway::Zero(tag_double)};
  // vector register to hold the epsilon24 values.
  VectorDouble _epsilon24{highway::Zero(tag_double)};
  // vector register to hold the sigmaSquared values.
  VectorDouble _sigmaSquared{highway::Zero(tag_double)};

  // cutoff squared used in the AoS functor.
  const double _cutoffSquareAoS{0.};
  // epsilon, sigma and shift6 used in the AoS functor.
  double _epsilon24AoS{0.}, _sigmaSquareAoS{0.}, _shift6AoS{0.};

  // optional to hold a reference to the ParticlePropertiesLibrary. If a ParticlePropertiesLibrary is not used the
  // optional is empty.
  std::optional<std::reference_wrapper<ParticlePropertiesLibrary<double, size_t>>> _PPLibrary;

  // accumulators for the globals (potential energy and virial).
  double _potentialEnergySum{0.};
  std::array<double, 3> _virialSum;
  std::vector<AoSThreadData> _aosThreadData{};

  // flag to indicate wether post processing has been performed.
  bool _postProcessed;

  // The Vectorization Pattern currently used in the SoA functor.
  VectorizationPattern _vecPattern;

  // Vectorization Pattern that the functor can handle.
  static constexpr std::array<VectorizationPattern, 4> _vecPatternsAllowed = {
      VectorizationPattern::p1xVec, VectorizationPattern::p2xVecDiv2, VectorizationPattern::pVecDiv2x2,
      VectorizationPattern::pVecx1};
};
}  // namespace mdLib