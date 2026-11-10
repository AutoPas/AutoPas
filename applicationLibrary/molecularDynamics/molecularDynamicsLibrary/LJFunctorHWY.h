/**
 * @file LJFunctorHWY.h
 *
 * @date 27.03.2024
 * @author Luis Gall
 */

#pragma once

#include <hwy/highway.h>

#include <algorithm>
#include <array>
#include <optional>
#include <vector>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/optRef.h"

namespace mdLib {

namespace highway = hwy::HWY_NAMESPACE;
/** Highway tag for full double register */
constexpr highway::ScalableTag<double> tag_double;
/** Highway tag for full long register */
constexpr highway::ScalableTag<int64_t> tag_long;
/**
 * Number of double values in a full register. Not constexpr on all architectures, so cannot be assumed, but, when
 * possible (e.g. on x86), we provide this constexpr so compilers may make optimizations.*/
HWY_LANES_CONSTEXPR inline size_t _vecLengthDouble{highway::Lanes(tag_double)};
/**
 * Upper bound on the number of double values in a full register. Always constexpr.
 */
constexpr size_t _maxVecLengthDouble{highway::MaxLanes(tag_double)};
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

  static constexpr auto AttrSqrt24Epsilon =
      static_cast<typename Particle_T::AttributeNames>(std::tuple_size_v<SoAArraysType> + 0);
  static constexpr auto AttrHalfSigma =
      static_cast<typename Particle_T::AttributeNames>(std::tuple_size_v<SoAArraysType> + 1);

 public:
  /**
   * Deleted default constructor
   */
  LJFunctorHWY() = delete;

  /**
   * Hook for SoALoader to load particle attributes or functor-computed mixing data.
   */
  template <typename Particle_T::AttributeNames attribute>
  inline auto loadSoAAttribute(const Particle_T &particle) const {
    if constexpr (attribute == AttrSqrt24Epsilon) {
      if constexpr (useMixing) {
        return _PPLibrary->get().getSqrt24Epsilon(particle.getTypeId());
      } else {
        return 0.0;
      }
    } else if constexpr (attribute == AttrHalfSigma) {
      if constexpr (useMixing) {
        return _PPLibrary->get().getHalfSigma(particle.getTypeId());
      } else {
        return 0.0;
      }
    } else {
      return particle.template get<attribute>();
    }
  }

  /**
   * Hook for SoALoader (non-const particle overload for ptr).
   */
  template <typename Particle_T::AttributeNames attribute>
  inline auto loadSoAAttribute(Particle_T &particle) const {
    if constexpr (attribute == AttrSqrt24Epsilon) {
      if constexpr (useMixing) {
        return _PPLibrary->get().getSqrt24Epsilon(particle.getTypeId());
      } else {
        return 0.0;
      }
    } else if constexpr (attribute == AttrHalfSigma) {
      if constexpr (useMixing) {
        return _PPLibrary->get().getHalfSigma(particle.getTypeId());
      } else {
        return 0.0;
      }
    } else {
      return particle.template get<attribute>();
    }
  }

  /**
   * Constructor for Functor with mixing enabled/disabled. When using this functor without particlePropertiesLibrary it
   * is necessary to call setParticleProperties() to set internal constants.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit LJFunctorHWY(double cutoff, std::optional<std::reference_wrapper<ParticlePropertiesLibrary<double, size_t>>>
                                           particlePropertiesLibrary = std::nullopt)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorHWY>(cutoff),
        _cutoffSquareAoS{cutoff * cutoff},
        _invCutoffSquareAoS{1.0 / (cutoff * cutoff)},
        _PPLibrary{particlePropertiesLibrary} {
    if constexpr (calculateGlobals) {
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
   * This functor can handle p1xVec, p2xVecDiv2, pVecDiv2x2, pVecx1
   *
   * @param vecPattern
   * @return whether the functor is capable of using the specified Vectorization Pattern
   */
  bool isVecPatternAllowed(const VectorizationPattern vecPattern) final {
    return std::ranges::find(_vecPatternsAllowed, vecPattern) != _vecPatternsAllowed.end();
  }

  /**
   * @copydoc autopas::PairwiseFunctor::AoSFunctor()
   */
  inline void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;
    if (i.isDummy() or j.isDummy()) {
      return;
    }
    const auto dr = i.getR() - j.getR();
    const double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffSquareAoS) {
      return;
    }

    auto sigmaSquare = _sigmaSquareAoS;
    auto epsilon24 = _epsilon24AoS;
    [[maybe_unused]] auto shift6 = _shift6AoS;
    if constexpr (useMixing) {
      const auto &mixingData = _PPLibrary->get().getLJMixingData(i.getTypeId(), j.getTypeId());
      sigmaSquare = mixingData.sigmaSquared;
      epsilon24 = mixingData.epsilon24;
      if constexpr (applyShift and calculateGlobals) {
        shift6 = mixingData.shift6;
      }
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
      j.subF(f);
    }
    if constexpr (calculateGlobals) {
      const auto virial = dr * f;
      const double potentialEnergy6 = epsilon24 * lj12m6 + shift6;

      const int threadID = autopas::autopas_get_thread_num();

      if (i.isOwned()) {
        _aosThreadData[threadID].potentialEnergySum += potentialEnergy6;
        _aosThreadData[threadID].virialSum += virial;
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        // for non-newton3 the division is in the post-processing step.
        _aosThreadData[threadID].potentialEnergySum += potentialEnergy6;
        _aosThreadData[threadID].virialSum += virial;
      }
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always do a newton3 like traversal of the SoA using only the 1xVec Pattern.
   */
  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, [[maybe_unused]] const bool newton3) final {
    const size_t soaSize = soa.size();
    if (soaSize <= 1) return;

    // obtain iterators for the various values
    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    const double *sqrtEpsilonPtr = nullptr;
    const double *halfSigmaPtr = nullptr;
    if constexpr (useMixing) {
      sqrtEpsilonPtr = soa.template begin<AttrSqrt24Epsilon>();
      halfSigmaPtr = soa.template begin<AttrHalfSigma>();
    }

    if (soaSize <= 7) {
      scalarSoAKernel<true, true>(soaSize, soaSize, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr,
                                  sqrtEpsilonPtr, halfSigmaPtr, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr,
                                  sqrtEpsilonPtr, halfSigmaPtr);
      return;
    }

    if constexpr (useMixing) {
      if (soaIsHomogeneous(sqrtEpsilonPtr, halfSigmaPtr, soaSize)) {
        runSingleLoop<true>(soaSize, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, sqrtEpsilonPtr,
                            halfSigmaPtr);
      } else {
        runSingleLoop<false>(soaSize, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, sqrtEpsilonPtr,
                             halfSigmaPtr);
      }
    } else {
      runSingleLoop<true>(soaSize, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, nullptr, nullptr);
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
   */
  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                             bool newton3) final {
    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }
    SoAFunctorPairDispatch<false>(soa1, soa2, newton3);
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPairSorted()
   */
  inline void SoAFunctorPairSorted(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                   const autopas::SoASortingData &sortingData, bool newton3) final {
    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }
    SoAFunctorPairDispatch<true>(soa1, soa2, newton3, sortingData);
  }

 private:
  /**
   * Number of i-particles processed per outer-loop step for the given vectorization pattern.
   * @tparam vecPattern
   */
  template <VectorizationPattern vecPattern>
  HWY_INLINE constexpr static size_t iStepSize() {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      return 1;
    }
    if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      return 2;
    }
    if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      return _vecLengthDouble / 2;
    }
    if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      return _vecLengthDouble;
    }
    autopas::utils::ExceptionHandler::exception("Unknown VectorizationPattern!");
    return {};
  }

  /**
   * Number of j-particles processed per inner-loop step for the given vectorization pattern.
   * @tparam vecPattern
   */
  template <VectorizationPattern vecPattern>
  HWY_INLINE constexpr static size_t jStepSize() {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      return _vecLengthDouble;
    }
    if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      return _vecLengthDouble / 2;
    }
    if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      return 2;
    }
    if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      return 1;
    }
    autopas::utils::ExceptionHandler::exception("Unknown VectorizationPattern!");
    return {};
  }

  /**
   * Fast check for whether all particles in an SoA share identical mixing parameters.
   */
  HWY_INLINE static bool soaIsHomogeneous(const double *const __restrict sqrtEpsilonPtr,
                                          const double *const __restrict halfSigmaPtr, const size_t soaSize) {
    if (soaSize <= 1) return true;
    const double firstSig = halfSigmaPtr[0];
    const double firstEps = sqrtEpsilonPtr[0];
    if (halfSigmaPtr[1] != firstSig or sqrtEpsilonPtr[1] != firstEps) {
      return false;
    }
    const VectorDouble vFirstSig = highway::Set(tag_double, firstSig);
    const VectorDouble vFirstEps = highway::Set(tag_double, firstEps);
    size_t k = 0;
    for (; k + _vecLengthDouble <= soaSize; k += _vecLengthDouble) {
      const auto vSig = highway::LoadU(tag_double, &halfSigmaPtr[k]);
      const auto vEps = highway::LoadU(tag_double, &sqrtEpsilonPtr[k]);
      const auto eqMask = highway::And(highway::Eq(vSig, vFirstSig), highway::Eq(vEps, vFirstEps));
      if (not highway::AllTrue(tag_double, eqMask)) {
        return false;
      }
    }
    for (; k < soaSize; ++k) {
      if (halfSigmaPtr[k] != firstSig or sqrtEpsilonPtr[k] != firstEps) {
        return false;
      }
    }
    return true;
  }

  /**
   * Computes scalar Lennard-Jones force and 6*potentialEnergy for small-N fallback.
   */
  HWY_INLINE void computeScalarLJInteraction(const double dr2, const double dx, const double dy, const double dz,
                                             [[maybe_unused]] const double sqrtEps1,
                                             [[maybe_unused]] const double halfSig1,
                                             [[maybe_unused]] const double sqrtEps2,
                                             [[maybe_unused]] const double halfSig2, double &fx, double &fy, double &fz,
                                             [[maybe_unused]] double &pot6) const {
    double sigmaSquared = _sigmaSquareAoS;
    double epsilon24 = _epsilon24AoS;
    [[maybe_unused]] double shift6 = _shift6AoS;
    if constexpr (useMixing) {
      const double sig = halfSig1 + halfSig2;
      sigmaSquared = sig * sig;
      epsilon24 = sqrtEps1 * sqrtEps2;
      if constexpr (calculateGlobals and applyShift) {
        shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquareAoS);
      }
    }
    const double invdr2 = 1.0 / dr2;
    const double lj2 = sigmaSquared * invdr2;
    const double lj6 = lj2 * lj2 * lj2;
    const double lj12 = lj6 * lj6;
    const double lj12m6 = lj12 - lj6;
    const double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    fx = dx * fac;
    fy = dy * fac;
    fz = dz * fac;
    if constexpr (calculateGlobals) {
      pot6 = epsilon24 * lj12m6 + shift6;
    }
  }

  /**
   * Scalar fallback kernel for very small particle counts (N <= 7), avoiding SIMD remainder overhead and horizontal
   * reduction trees when only a few lanes would be active.
   */
  template <bool isSingleCell, bool newton3>
  HWY_INLINE void scalarSoAKernel(const size_t soaSize1, const size_t soaSize2, const double *const __restrict x1Ptr,
                                  const double *const __restrict y1Ptr, const double *const __restrict z1Ptr,
                                  const autopas::OwnershipState *const __restrict ownedStatePtr1,
                                  double *const __restrict fx1Ptr, double *const __restrict fy1Ptr,
                                  double *const __restrict fz1Ptr, const double *const __restrict sqrtEpsilon1Ptr,
                                  const double *const __restrict halfSigma1Ptr, const double *const __restrict x2Ptr,
                                  const double *const __restrict y2Ptr, const double *const __restrict z2Ptr,
                                  const autopas::OwnershipState *const __restrict ownedStatePtr2,
                                  double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                                  double *const __restrict fz2Ptr, const double *const __restrict sqrtEpsilon2Ptr,
                                  const double *const __restrict halfSigma2Ptr) {
    [[maybe_unused]] double virialX = 0., virialY = 0., virialZ = 0., potentialEnergy = 0.;

    std::array<double, 8> fx2Acc{}, fy2Acc{}, fz2Acc{};

    constexpr size_t startI = isSingleCell ? 1 : 0;
    for (size_t i = startI; i < soaSize1; ++i) {
      if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
        continue;
      }
      const double xi = x1Ptr[i], yi = y1Ptr[i], zi = z1Ptr[i];
      double sqrtEps1 = 0., halfSig1 = 0.;
      if constexpr (useMixing) {
        sqrtEps1 = sqrtEpsilon1Ptr[i];
        halfSig1 = halfSigma1Ptr[i];
      }
      [[maybe_unused]] const bool isOwnedI = (ownedStatePtr1[i] == autopas::OwnershipState::owned);

      double fxi = 0., fyi = 0., fzi = 0.;
      const size_t endJ = isSingleCell ? i : soaSize2;
      for (size_t j = 0; j < endJ; ++j) {
        if (ownedStatePtr2[j] == autopas::OwnershipState::dummy) {
          continue;
        }
        const double dx = xi - x2Ptr[j];
        const double dy = yi - y2Ptr[j];
        const double dz = zi - z2Ptr[j];
        const double dr2 = dx * dx + dy * dy + dz * dz;

        if (dr2 <= _cutoffSquareAoS) {
          double fx = 0., fy = 0., fz = 0., pot6 = 0.;
          computeScalarLJInteraction(dr2, dx, dy, dz, sqrtEps1, halfSig1, useMixing ? sqrtEpsilon2Ptr[j] : 0.,
                                     useMixing ? halfSigma2Ptr[j] : 0., fx, fy, fz, pot6);
          fxi += fx;
          fyi += fy;
          fzi += fz;
          if constexpr (newton3) {
            fx2Acc[j] -= fx;
            fy2Acc[j] -= fy;
            fz2Acc[j] -= fz;
          }
          if constexpr (calculateGlobals) {
            const bool isOwnedJ = newton3 and (ownedStatePtr2[j] == autopas::OwnershipState::owned);
            const double energyFactor = (isOwnedI ? 1.0 : 0.0) + (isOwnedJ ? 1.0 : 0.0);
            potentialEnergy += energyFactor * pot6;
            virialX += energyFactor * fx * dx;
            virialY += energyFactor * fy * dy;
            virialZ += energyFactor * fz * dz;
          }
        }
      }
      if constexpr (isSingleCell) {
        fx2Acc[i] += fxi;
        fy2Acc[i] += fyi;
        fz2Acc[i] += fzi;
      } else {
        fx1Ptr[i] += fxi;
        fy1Ptr[i] += fyi;
        fz1Ptr[i] += fzi;
      }
    }

    if constexpr (isSingleCell) {
      for (size_t i = 0; i < soaSize1; ++i) {
        fx1Ptr[i] += fx2Acc[i];
        fy1Ptr[i] += fy2Acc[i];
        fz1Ptr[i] += fz2Acc[i];
      }
    } else if constexpr (newton3) {
      for (size_t j = 0; j < soaSize2; ++j) {
        fx2Ptr[j] += fx2Acc[j];
        fy2Ptr[j] += fy2Acc[j];
        fz2Ptr[j] += fz2Acc[j];
      }
    }

    if constexpr (calculateGlobals) {
      const int threadID = autopas::autopas_get_thread_num();
      _aosThreadData[threadID].potentialEnergySum += potentialEnergy;
      _aosThreadData[threadID].virialSum[0] += virialX;
      _aosThreadData[threadID].virialSum[1] += virialY;
      _aosThreadData[threadID].virialSum[2] += virialZ;
    }
  }

  /**
   * Depending on the vectorization pattern, loads/broadcasts one scalar array (positions, ownership state, or mixing
   * properties) into a full vector register. This is the shared primitive behind fillIRegisters and fillJRegisters:
   * every VectorizationPattern reduces, per side, to one of two idioms - a "paired" idiom
   * (broadcast the current element, and, unless in the remainder case, its neighbor, then combine both halves) or a
   * "contiguous" idiom (load a run of elements, self-combining to fill the register if it is only a half-load).
   * Which idiom applies to which side is fixed by vecPattern: e.g. for p2xVecDiv2 the i-side is paired (2 wide) and
   * the j-side is contiguous (Vec/2 wide); for pVecDiv2x2 it is the other way around.
   *
   * @tparam isISide Whether this call fills a register for the i-side (true) or the j-side (false) of vecPattern.
   * @tparam remainder Whether this is the remainder-loop case, i.e. fewer than a full step of elements are valid.
   * @tparam reversed Whether the i-loop iterates backwards. Only affects the i-side's neighbor direction in the
   * paired idiom; ignored for the j-side, which always takes its neighbor forwards.
   * @tparam vecPattern
   * @param tag Highway tag selecting the element type (e.g. tag_double or tag_long), and thus the vector type.
   * @param basePtr Pointer to the scalar array, unoffset (i.e. callers pass e.g. xPtr, not &xPtr[i]).
   * @param idx The current loop index (i or j) into basePtr.
   * @param rest The number of valid lanes in the remainder case.
   */
  template <bool isISide, bool remainder, bool reversed, VectorizationPattern vecPattern, class TagType, class PtrType>
  HWY_INLINE static auto loadRegister(TagType tag, const PtrType *const __restrict basePtr, const size_t idx,
                                      const size_t rest) {
    if constexpr (isISide) {
      if constexpr (vecPattern == VectorizationPattern::p1xVec) {
        return highway::Set(tag, basePtr[idx]);
      } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
        const auto cur = highway::Set(tag, basePtr[idx]);
        auto nbr = highway::Zero(tag);
        if constexpr (not remainder) {
          const auto nbrIdx = reversed ? idx - 1 : idx + 1;
          nbr = highway::Set(tag, basePtr[nbrIdx]);
        }
        return highway::ConcatLowerLower(tag, nbr, cur);
      } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
        const size_t windowStart = reversed ? (remainder ? 0 : idx - _vecLengthDouble / 2 + 1) : idx;
        const size_t lanes = remainder ? rest : _vecLengthDouble / 2;
        const auto v = highway::LoadN(tag, &basePtr[windowStart], lanes);
        return highway::ConcatLowerLower(tag, v, v);
      } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
        const size_t windowStart = reversed ? (remainder ? 0 : idx - _vecLengthDouble + 1) : idx;
        if constexpr (remainder) {
          return highway::LoadN(tag, &basePtr[windowStart], rest);
        } else {
          return highway::LoadU(tag, &basePtr[windowStart]);
        }
      }
    } else {
      if constexpr (vecPattern == VectorizationPattern::p1xVec) {
        if constexpr (remainder) {
          return highway::LoadN(tag, &basePtr[idx], rest);
        } else {
          return highway::LoadU(tag, &basePtr[idx]);
        }
      } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
        const size_t lanes = remainder ? rest : _vecLengthDouble / 2;
        const auto v = highway::LoadN(tag, &basePtr[idx], lanes);
        return highway::ConcatLowerLower(tag, v, v);
      } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
        const auto cur = highway::Set(tag, basePtr[idx]);
        if constexpr (remainder) {
          return highway::ConcatLowerLower(tag, highway::Zero(tag), cur);
        } else {
          const auto nbr = highway::Set(tag, basePtr[idx + 1]);
          return highway::ConcatLowerLower(tag, nbr, cur);
        }
      } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
        return highway::Set(tag, basePtr[idx]);
      }
    }
  }

  /**
   * Recursively folds a vector in half down to 2 lanes by adding its lower and upper halves.
   */
  template <class D>
  HWY_INLINE static auto foldTo2Lanes(D d, const highway::VFromD<D> &v) {
    if constexpr (highway::MaxLanes(d) == 8) {
      const highway::Half<D> d4;
      const auto v4 = highway::Add(highway::LowerHalf(d4, v), highway::UpperHalf(d4, v));
      const highway::Half<decltype(d4)> d2;
      return highway::Add(highway::LowerHalf(d2, v4), highway::UpperHalf(d2, v4));
    } else if constexpr (highway::MaxLanes(d) == 4) {
      const highway::Half<D> d2;
      return highway::Add(highway::LowerHalf(d2, v), highway::UpperHalf(d2, v));
    } else {
      return v;
    }
  }

  /**
   * Reduces a single vector to a scalar using a halving tree down to 2 lanes before calling highway::ReduceSum.
   */
  template <class D>
  HWY_INLINE static double reduceSumHalving(D d, const highway::VFromD<D> &v) {
    if constexpr (highway::MaxLanes(d) <= 2) {
      return highway::ReduceSum(d, v);
    } else {
      const auto v2 = foldTo2Lanes(d, v);
      const highway::DFromV<decltype(v2)> d2;
      return highway::ReduceSum(d2, v2);
    }
  }

  /**
   * Reduces two vectors simultaneously by interleaving their even/odd lanes and folding to 2 lanes,
   * adding the results to output scalar accumulators.
   */
  template <class D>
  HWY_INLINE static void reduceAdd2_InterleaveFirst(D d, const highway::VFromD<D> &va, const highway::VFromD<D> &vb,
                                                    double &outA, double &outB) {
    if constexpr (highway::MaxLanes(d) == 1) {
      outA += highway::GetLane(va);
      outB += highway::GetLane(vb);
    } else {
      const auto ab = highway::Add(highway::InterleaveLower(d, va, vb), highway::InterleaveUpper(d, va, vb));
      const auto ab2 = foldTo2Lanes(d, ab);
      outA += highway::GetLane(ab2);
      outB += highway::ExtractLane(ab2, 1);
    }
  }

  /**
   * Reduces two vectors simultaneously by interleaving their even/odd lanes and folding to 2 lanes,
   * subtracting the results from output scalar accumulators.
   */
  template <class D>
  HWY_INLINE static void reduceSub2_InterleaveFirst(D d, const highway::VFromD<D> &va, const highway::VFromD<D> &vb,
                                                    double &outA, double &outB) {
    if constexpr (highway::MaxLanes(d) == 1) {
      outA -= highway::GetLane(va);
      outB -= highway::GetLane(vb);
    } else {
      const auto ab = highway::Add(highway::InterleaveLower(d, va, vb), highway::InterleaveUpper(d, va, vb));
      const auto ab2 = foldTo2Lanes(d, ab);
      outA -= highway::GetLane(ab2);
      outB -= highway::ExtractLane(ab2, 1);
    }
  }

  /**
   * Reduces 3 vectors (x, y, z) into scalar accumulators by interleaving (x, y) and halving z.
   */
  template <class D>
  HWY_INLINE static void reduceAdd3_InterleaveFirst(D d, const highway::VFromD<D> &fx, const highway::VFromD<D> &fy,
                                                    const highway::VFromD<D> &fz, double &outX, double &outY,
                                                    double &outZ) {
    reduceAdd2_InterleaveFirst(d, fx, fy, outX, outY);
    outZ += reduceSumHalving(d, fz);
  }

  /**
   * Reduces 3 vectors (x, y, z) into scalar accumulators by interleaving (x, y) and halving z,
   * subtracting the results from output scalar accumulators.
   */
  template <class D>
  HWY_INLINE static void reduceSub3_InterleaveFirst(D d, const highway::VFromD<D> &fx, const highway::VFromD<D> &fy,
                                                    const highway::VFromD<D> &fz, double &outX, double &outY,
                                                    double &outZ) {
    reduceSub2_InterleaveFirst(d, fx, fy, outX, outY);
    outZ -= reduceSumHalving(d, fz);
  }

  /**
   * Reduces 4 vectors into scalar accumulators using two parallel interleaved pairs.
   */
  template <class D>
  HWY_INLINE static void reduceAdd4_Interleaved(D d, const highway::VFromD<D> &v0, const highway::VFromD<D> &v1,
                                                const highway::VFromD<D> &v2, const highway::VFromD<D> &v3,
                                                double &out0, double &out1, double &out2, double &out3) {
    reduceAdd2_InterleaveFirst(d, v0, v1, out0, out1);
    reduceAdd2_InterleaveFirst(d, v2, v3, out2, out3);
  }

  /**
   * Pairs lower and upper half-vectors for 3 dimensions into 2-lane vectors and adds them directly
   * to memory at 2 contiguous particle positions using 128-bit vector stores.
   */
  template <class DHalf>
  HWY_INLINE static void reduceAdd6Half_Vec2Store(DHalf dh, const highway::VFromD<DHalf> &lfx,
                                                  const highway::VFromD<DHalf> &lfy, const highway::VFromD<DHalf> &lfz,
                                                  const highway::VFromD<DHalf> &ufx, const highway::VFromD<DHalf> &ufy,
                                                  const highway::VFromD<DHalf> &ufz, double *const __restrict fxPtr,
                                                  double *const __restrict fyPtr, double *const __restrict fzPtr) {
    if constexpr (highway::MaxLanes(dh) == 1) {
      fxPtr[0] += highway::GetLane(lfx);
      fyPtr[0] += highway::GetLane(lfy);
      fzPtr[0] += highway::GetLane(lfz);
      fxPtr[1] += highway::GetLane(ufx);
      fyPtr[1] += highway::GetLane(ufy);
      fzPtr[1] += highway::GetLane(ufz);
    } else {
      const auto fxLU = highway::Add(highway::InterleaveLower(dh, lfx, ufx), highway::InterleaveUpper(dh, lfx, ufx));
      const auto fyLU = highway::Add(highway::InterleaveLower(dh, lfy, ufy), highway::InterleaveUpper(dh, lfy, ufy));
      const auto fzLU = highway::Add(highway::InterleaveLower(dh, lfz, ufz), highway::InterleaveUpper(dh, lfz, ufz));
      const auto fx2 = foldTo2Lanes(dh, fxLU);
      const auto fy2 = foldTo2Lanes(dh, fyLU);
      const auto fz2 = foldTo2Lanes(dh, fzLU);
      const highway::DFromV<decltype(fx2)> d2;

      highway::StoreU(highway::Add(highway::LoadU(d2, fxPtr), fx2), d2, fxPtr);
      highway::StoreU(highway::Add(highway::LoadU(d2, fyPtr), fy2), d2, fyPtr);
      highway::StoreU(highway::Add(highway::LoadU(d2, fzPtr), fz2), d2, fzPtr);
    }
  }

  /**
   * Pairs lower and upper half-vectors for 3 dimensions into 2-lane vectors and subtracts them directly
   * from memory at 2 contiguous particle positions using 128-bit vector stores.
   */
  template <class DHalf>
  HWY_INLINE static void reduceSub6Half_Vec2Store(DHalf dh, const highway::VFromD<DHalf> &lfx,
                                                  const highway::VFromD<DHalf> &lfy, const highway::VFromD<DHalf> &lfz,
                                                  const highway::VFromD<DHalf> &ufx, const highway::VFromD<DHalf> &ufy,
                                                  const highway::VFromD<DHalf> &ufz, double *const __restrict fxPtr,
                                                  double *const __restrict fyPtr, double *const __restrict fzPtr) {
    if constexpr (highway::MaxLanes(dh) == 1) {
      fxPtr[0] -= highway::GetLane(lfx);
      fyPtr[0] -= highway::GetLane(lfy);
      fzPtr[0] -= highway::GetLane(lfz);
      fxPtr[1] -= highway::GetLane(ufx);
      fyPtr[1] -= highway::GetLane(ufy);
      fzPtr[1] -= highway::GetLane(ufz);
    } else {
      const auto fxLU = highway::Add(highway::InterleaveLower(dh, lfx, ufx), highway::InterleaveUpper(dh, lfx, ufx));
      const auto fyLU = highway::Add(highway::InterleaveLower(dh, lfy, ufy), highway::InterleaveUpper(dh, lfy, ufy));
      const auto fzLU = highway::Add(highway::InterleaveLower(dh, lfz, ufz), highway::InterleaveUpper(dh, lfz, ufz));
      const auto fx2 = foldTo2Lanes(dh, fxLU);
      const auto fy2 = foldTo2Lanes(dh, fyLU);
      const auto fz2 = foldTo2Lanes(dh, fzLU);
      const highway::DFromV<decltype(fx2)> d2;

      highway::StoreU(highway::Sub(highway::LoadU(d2, fxPtr), fx2), d2, fxPtr);
      highway::StoreU(highway::Sub(highway::LoadU(d2, fyPtr), fy2), d2, fyPtr);
      highway::StoreU(highway::Sub(highway::LoadU(d2, fzPtr), fz2), d2, fzPtr);
    }
  }

  template <bool remainder, VectorizationPattern vecPattern>
  HWY_INLINE static void handleNewton3Reduction(const VectorDouble &fx, const VectorDouble &fy, const VectorDouble &fz,
                                                double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                                                double *const __restrict fz2Ptr, const size_t j, const size_t rest) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      const VectorDouble fx2 =
          remainder ? highway::LoadN(tag_double, &fx2Ptr[j], rest) : highway::LoadU(tag_double, &fx2Ptr[j]);
      const VectorDouble fy2 =
          remainder ? highway::LoadN(tag_double, &fy2Ptr[j], rest) : highway::LoadU(tag_double, &fy2Ptr[j]);
      const VectorDouble fz2 =
          remainder ? highway::LoadN(tag_double, &fz2Ptr[j], rest) : highway::LoadU(tag_double, &fz2Ptr[j]);

      const VectorDouble fx2New = highway::Sub(fx2, fx);
      const VectorDouble fy2New = highway::Sub(fy2, fy);
      const VectorDouble fz2New = highway::Sub(fz2, fz);

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

      const auto fxCombined = highway::Add(lowerFx, upperFx);
      const auto fyCombined = highway::Add(lowerFy, upperFy);
      const auto fzCombined = highway::Add(lowerFz, upperFz);

      const size_t lanes = remainder ? rest : _vecLengthDouble / 2;

      const auto fx2 =
          remainder ? highway::LoadN(tag_double_half, &fx2Ptr[j], lanes) : highway::LoadU(tag_double_half, &fx2Ptr[j]);
      const auto fy2 =
          remainder ? highway::LoadN(tag_double_half, &fy2Ptr[j], lanes) : highway::LoadU(tag_double_half, &fy2Ptr[j]);
      const auto fz2 =
          remainder ? highway::LoadN(tag_double_half, &fz2Ptr[j], lanes) : highway::LoadU(tag_double_half, &fz2Ptr[j]);

      const auto newFx = highway::Sub(fx2, fxCombined);
      const auto newFy = highway::Sub(fy2, fyCombined);
      const auto newFz = highway::Sub(fz2, fzCombined);

      if constexpr (remainder) {
        highway::StoreN(newFx, tag_double_half, &fx2Ptr[j], lanes);
        highway::StoreN(newFy, tag_double_half, &fy2Ptr[j], lanes);
        highway::StoreN(newFz, tag_double_half, &fz2Ptr[j], lanes);
      } else {
        highway::StoreU(newFx, tag_double_half, &fx2Ptr[j]);
        highway::StoreU(newFy, tag_double_half, &fy2Ptr[j]);
        highway::StoreU(newFz, tag_double_half, &fz2Ptr[j]);
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      const auto lowerFx = highway::LowerHalf(tag_double_half, fx);
      const auto lowerFy = highway::LowerHalf(tag_double_half, fy);
      const auto lowerFz = highway::LowerHalf(tag_double_half, fz);

      if constexpr (remainder) {
        reduceSub3_InterleaveFirst(tag_double_half, lowerFx, lowerFy, lowerFz, fx2Ptr[j], fy2Ptr[j], fz2Ptr[j]);
      } else {
        const auto upperFx = highway::UpperHalf(tag_double_half, fx);
        const auto upperFy = highway::UpperHalf(tag_double_half, fy);
        const auto upperFz = highway::UpperHalf(tag_double_half, fz);

        reduceSub6Half_Vec2Store(tag_double_half, lowerFx, lowerFy, lowerFz, upperFx, upperFy, upperFz, &fx2Ptr[j],
                                 &fy2Ptr[j], &fz2Ptr[j]);
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      reduceSub3_InterleaveFirst(tag_double, fx, fy, fz, fx2Ptr[j], fy2Ptr[j], fz2Ptr[j]);
    }
  }

  template <bool reversed, bool remainder, VectorizationPattern vecPattern>
  HWY_INLINE static void reduceAccumulatedForce(const size_t i, double *const __restrict fxPtr,
                                                double *const __restrict fyPtr, double *const __restrict fzPtr,
                                                const VectorDouble &fxAcc, const VectorDouble &fyAcc,
                                                const VectorDouble &fzAcc, const size_t restI) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      reduceAdd3_InterleaveFirst(tag_double, fxAcc, fyAcc, fzAcc, fxPtr[i], fyPtr[i], fzPtr[i]);
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      const auto lowerFxAcc = highway::LowerHalf(tag_double_half, fxAcc);
      const auto lowerFyAcc = highway::LowerHalf(tag_double_half, fyAcc);
      const auto lowerFzAcc = highway::LowerHalf(tag_double_half, fzAcc);

      if constexpr (remainder) {
        reduceAdd3_InterleaveFirst(tag_double_half, lowerFxAcc, lowerFyAcc, lowerFzAcc, fxPtr[i], fyPtr[i], fzPtr[i]);
      } else {
        const auto upperFxAcc = highway::UpperHalf(tag_double_half, fxAcc);
        const auto upperFyAcc = highway::UpperHalf(tag_double_half, fyAcc);
        const auto upperFzAcc = highway::UpperHalf(tag_double_half, fzAcc);

        if constexpr (reversed) {
          reduceAdd6Half_Vec2Store(tag_double_half, upperFxAcc, upperFyAcc, upperFzAcc, lowerFxAcc, lowerFyAcc,
                                   lowerFzAcc, &fxPtr[i - 1], &fyPtr[i - 1], &fzPtr[i - 1]);
        } else {
          reduceAdd6Half_Vec2Store(tag_double_half, lowerFxAcc, lowerFyAcc, lowerFzAcc, upperFxAcc, upperFyAcc,
                                   upperFzAcc, &fxPtr[i], &fyPtr[i], &fzPtr[i]);
        }
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      const auto lowerFxAcc = highway::LowerHalf(tag_double_half, fxAcc);
      const auto lowerFyAcc = highway::LowerHalf(tag_double_half, fyAcc);
      const auto lowerFzAcc = highway::LowerHalf(tag_double_half, fzAcc);

      const auto upperFxAcc = highway::UpperHalf(tag_double_half, fxAcc);
      const auto upperFyAcc = highway::UpperHalf(tag_double_half, fyAcc);
      const auto upperFzAcc = highway::UpperHalf(tag_double_half, fzAcc);

      const auto fxAccCombined = highway::Add(lowerFxAcc, upperFxAcc);
      const auto fyAccCombined = highway::Add(lowerFyAcc, upperFyAcc);
      const auto fzAccCombined = highway::Add(lowerFzAcc, upperFzAcc);

      const size_t index = reversed ? (remainder ? 0 : i - _vecLengthDouble / 2 + 1) : i;
      const size_t lanes = remainder ? restI : _vecLengthDouble / 2;

      const auto oldFx = remainder ? highway::LoadN(tag_double_half, &fxPtr[index], lanes)
                                   : highway::LoadU(tag_double_half, &fxPtr[index]);
      const auto oldFy = remainder ? highway::LoadN(tag_double_half, &fyPtr[index], lanes)
                                   : highway::LoadU(tag_double_half, &fyPtr[index]);
      const auto oldFz = remainder ? highway::LoadN(tag_double_half, &fzPtr[index], lanes)
                                   : highway::LoadU(tag_double_half, &fzPtr[index]);

      const auto newFx = highway::Add(oldFx, fxAccCombined);
      const auto newFy = highway::Add(oldFy, fyAccCombined);
      const auto newFz = highway::Add(oldFz, fzAccCombined);

      if constexpr (remainder) {
        highway::StoreN(newFx, tag_double_half, &fxPtr[index], lanes);
        highway::StoreN(newFy, tag_double_half, &fyPtr[index], lanes);
        highway::StoreN(newFz, tag_double_half, &fzPtr[index], lanes);
      } else {
        highway::StoreU(newFx, tag_double_half, &fxPtr[index]);
        highway::StoreU(newFy, tag_double_half, &fyPtr[index]);
        highway::StoreU(newFz, tag_double_half, &fzPtr[index]);
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      const size_t index = reversed ? (remainder ? 0 : i - _vecLengthDouble + 1) : i;
      const VectorDouble oldFx =
          remainder ? highway::LoadN(tag_double, &fxPtr[index], restI) : highway::LoadU(tag_double, &fxPtr[index]);
      const VectorDouble oldFy =
          remainder ? highway::LoadN(tag_double, &fyPtr[index], restI) : highway::LoadU(tag_double, &fyPtr[index]);
      const VectorDouble oldFz =
          remainder ? highway::LoadN(tag_double, &fzPtr[index], restI) : highway::LoadU(tag_double, &fzPtr[index]);

      const VectorDouble fxNew = highway::Add(oldFx, fxAcc);
      const VectorDouble fyNew = highway::Add(oldFy, fyAcc);
      const VectorDouble fzNew = highway::Add(oldFz, fzAcc);

      remainder ? highway::StoreN(fxNew, tag_double, &fxPtr[index], restI)
                : highway::StoreU(fxNew, tag_double, &fxPtr[index]);
      remainder ? highway::StoreN(fyNew, tag_double, &fyPtr[index], restI)
                : highway::StoreU(fyNew, tag_double, &fyPtr[index]);
      remainder ? highway::StoreN(fzNew, tag_double, &fzPtr[index], restI)
                : highway::StoreU(fzNew, tag_double, &fzPtr[index]);
    }
  }

  HWY_INLINE void computeGlobals(const VectorDouble &virialSumX, const VectorDouble &virialSumY,
                                 const VectorDouble &virialSumZ, const VectorDouble &uPotSum) {
    const int threadID = autopas::autopas_get_thread_num();

    reduceAdd4_Interleaved(tag_double, virialSumX, virialSumY, virialSumZ, uPotSum,
                           _aosThreadData[threadID].virialSum[0], _aosThreadData[threadID].virialSum[1],
                           _aosThreadData[threadID].virialSum[2], _aosThreadData[threadID].potentialEnergySum);
  }

  /**
   * Helper function that handles an iteration of the outer-loop (the loop that changes "i").
   * @tparam reversed Whether iterating backwards over the outer-loop.
   * @tparam newton3
   * @tparam remainderI True iff the "i" vector register is not fully filled, i.e. there are less than _vecLength
   * elements of the SoA to iterate over.
   * @tparam vecPattern
   * @tparam homogeneousParams Whether all interacting pairs in this call share constant LJ parameters.
   */
  template <bool reversed, bool newton3, bool remainderI, VectorizationPattern vecPattern,
            bool homogeneousParams = not useMixing>
  HWY_INLINE void handleILoopBody(
      const size_t i, const double *const __restrict xPtr1, const double *const __restrict yPtr1,
      const double *const __restrict zPtr1, const autopas::OwnershipState *const __restrict ownedStatePtr1,
      const double *const __restrict xPtr2, const double *const __restrict yPtr2, const double *const __restrict zPtr2,
      const autopas::OwnershipState *const __restrict ownedStatePtr2, double *const __restrict fxPtr1,
      double *const __restrict fyPtr1, double *const __restrict fzPtr1, double *const __restrict fxPtr2,
      double *const __restrict fyPtr2, double *const __restrict fzPtr2, const double *const __restrict sqrtEpsilonPtr1,
      const double *const __restrict halfSigmaPtr1, const double *const __restrict sqrtEpsilonPtr2,
      const double *const __restrict halfSigmaPtr2, VectorDouble &virialSumX, VectorDouble &virialSumY,
      VectorDouble &virialSumZ, VectorDouble &uPotSum, const size_t restI, const size_t jVecStart,
      const size_t jVecEnd) const {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) [[unlikely]] {
        return;
      }
    }

    MaskDouble notDummyMaskI;
    [[maybe_unused]] VectorDouble energyFactorI = highway::Zero(tag_double);
    [[maybe_unused]] VectorDouble energyFactorIPlus1 = highway::Zero(tag_double);

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      if constexpr (calculateGlobals) {
        const double baseFactor = (ownedStatePtr1[i] == autopas::OwnershipState::owned) ? 1.0 : 0.0;
        energyFactorI = highway::Set(tag_double, baseFactor);
        if constexpr (newton3) {
          energyFactorIPlus1 = highway::Set(tag_double, baseFactor + 1.0);
        }
      }
    } else {
      const VectorLong ownedStateI = loadRegister<true, remainderI, reversed, vecPattern>(
          tag_long, reinterpret_cast<const int64_t *>(ownedStatePtr1), i, restI);
      notDummyMaskI = highway::RebindMask(tag_double, highway::Ne(ownedStateI, highway::Zero(tag_long)));
      if (highway::AllFalse(tag_double, notDummyMaskI)) [[unlikely]] {
        return;
      }
      if constexpr (calculateGlobals) {
        const auto ownedVal = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
        const auto ownedMaskI = highway::RebindMask(tag_double, highway::Eq(ownedStateI, ownedVal));
        energyFactorI = highway::MaskedSet(tag_double, ownedMaskI, 1.0);
        if constexpr (newton3) {
          energyFactorIPlus1 = highway::Add(energyFactorI, highway::Set(tag_double, 1.0));
        }
      }
    }

    const VectorDouble x1 = loadRegister<true, remainderI, reversed, vecPattern>(tag_double, xPtr1, i, restI);
    const VectorDouble y1 = loadRegister<true, remainderI, reversed, vecPattern>(tag_double, yPtr1, i, restI);
    const VectorDouble z1 = loadRegister<true, remainderI, reversed, vecPattern>(tag_double, zPtr1, i, restI);

    VectorDouble sqrtEpsilon1;
    VectorDouble sigmaHalf1;
    if constexpr (homogeneousParams) {
      if constexpr (useMixing) {
        const double sig = halfSigmaPtr1[0] + halfSigmaPtr2[0];
        sqrtEpsilon1 = highway::Set(tag_double, sqrtEpsilonPtr1[0] * sqrtEpsilonPtr2[0]);
        sigmaHalf1 = highway::Set(tag_double, sig * sig);
      } else {
        sqrtEpsilon1 = highway::Set(tag_double, _epsilon24AoS);
        sigmaHalf1 = highway::Set(tag_double, _sigmaSquareAoS);
      }
    } else {
      sqrtEpsilon1 = loadRegister<true, remainderI, reversed, vecPattern>(tag_double, sqrtEpsilonPtr1, i, restI);
      sigmaHalf1 = loadRegister<true, remainderI, reversed, vecPattern>(tag_double, halfSigmaPtr1, i, restI);
    }

    VectorDouble fxAcc = highway::Zero(tag_double);
    VectorDouble fyAcc = highway::Zero(tag_double);
    VectorDouble fzAcc = highway::Zero(tag_double);

    constexpr size_t jStep = jStepSize<vecPattern>();
    const size_t restJ = jVecEnd & (jStep - 1);
    const size_t jVecEndFull = jVecEnd - restJ;
    size_t j = jVecStart;

    for (; j < jVecEndFull; j += jStep) {
      SoAKernel<newton3, false, reversed, vecPattern, homogeneousParams>(
          j, notDummyMaskI, energyFactorI, energyFactorIPlus1, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1,
          y1, z1, xPtr2, yPtr2, zPtr2, fxPtr2, fyPtr2, fzPtr2, sqrtEpsilon1, sigmaHalf1, sqrtEpsilonPtr2, halfSigmaPtr2,
          fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, 0);
    }

    if (restJ > 0) {
      SoAKernel<newton3, true, reversed, vecPattern, homogeneousParams>(
          j, notDummyMaskI, energyFactorI, energyFactorIPlus1, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1,
          y1, z1, xPtr2, yPtr2, zPtr2, fxPtr2, fyPtr2, fzPtr2, sqrtEpsilon1, sigmaHalf1, sqrtEpsilonPtr2, halfSigmaPtr2,
          fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum, restJ);
    }

    reduceAccumulatedForce<reversed, remainderI, vecPattern>(i, fxPtr1, fyPtr1, fzPtr1, fxAcc, fyAcc, fzAcc, restI);
  }

  /**
   * Dispatches SoAFunctorPair or SoAFunctorPairSorted to scalar small-N fallback or SIMD VectorizationPattern.
   */
  template <bool sorted>
  inline void SoAFunctorPairDispatch(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                     const bool newton3,
                                     autopas::utils::optRef<const autopas::SoASortingData> sortingData = std::nullopt) {
    const size_t n1 = soa1.size();
    const size_t n2 = soa2.size();

    if constexpr (useMixing) {
      if (n1 <= 7 and n2 <= 7) {
        if constexpr (sorted) {
          if (sortingData.has_value()) {
            const auto &sd = sortingData->get();
            if (sd.startI >= sd.endI) {
              return;
            }
          }
        }
        const auto *const __restrict x1Ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
        const auto *const __restrict y1Ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
        const auto *const __restrict z1Ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
        const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle_T::AttributeNames::ownershipState>();
        auto *const __restrict fx1Ptr = soa1.template begin<Particle_T::AttributeNames::forceX>();
        auto *const __restrict fy1Ptr = soa1.template begin<Particle_T::AttributeNames::forceY>();
        auto *const __restrict fz1Ptr = soa1.template begin<Particle_T::AttributeNames::forceZ>();
        const auto *const __restrict sqrtEpsilon1Ptr = soa1.template begin<AttrSqrt24Epsilon>();
        const auto *const __restrict halfSigma1Ptr = soa1.template begin<AttrHalfSigma>();

        const auto *const __restrict x2Ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
        const auto *const __restrict y2Ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
        const auto *const __restrict z2Ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();
        const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle_T::AttributeNames::ownershipState>();
        auto *const __restrict fx2Ptr = soa2.template begin<Particle_T::AttributeNames::forceX>();
        auto *const __restrict fy2Ptr = soa2.template begin<Particle_T::AttributeNames::forceY>();
        auto *const __restrict fz2Ptr = soa2.template begin<Particle_T::AttributeNames::forceZ>();
        const auto *const __restrict sqrtEpsilon2Ptr = soa2.template begin<AttrSqrt24Epsilon>();
        const auto *const __restrict halfSigma2Ptr = soa2.template begin<AttrHalfSigma>();

        if (newton3) {
          scalarSoAKernel<false, true>(n1, n2, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, fx1Ptr, fy1Ptr, fz1Ptr,
                                       sqrtEpsilon1Ptr, halfSigma1Ptr, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr,
                                       fy2Ptr, fz2Ptr, sqrtEpsilon2Ptr, halfSigma2Ptr);
        } else {
          scalarSoAKernel<false, false>(n1, n2, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, fx1Ptr, fy1Ptr, fz1Ptr,
                                        sqrtEpsilon1Ptr, halfSigma1Ptr, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr,
                                        fy2Ptr, fz2Ptr, sqrtEpsilon2Ptr, halfSigma2Ptr);
        }
        return;
      }
    }

    switch (_vecPattern) {
      case VectorizationPattern::p1xVec: {
        if (newton3) {
          SoAFunctorPairImpl<true, sorted, VectorizationPattern::p1xVec>(soa1, soa2, sortingData);
        } else {
          SoAFunctorPairImpl<false, sorted, VectorizationPattern::p1xVec>(soa1, soa2, sortingData);
        }
        break;
      }
      case VectorizationPattern::p2xVecDiv2: {
        if (newton3) {
          SoAFunctorPairImpl<true, sorted, VectorizationPattern::p2xVecDiv2>(soa1, soa2, sortingData);
        } else {
          SoAFunctorPairImpl<false, sorted, VectorizationPattern::p2xVecDiv2>(soa1, soa2, sortingData);
        }
        break;
      }
      case VectorizationPattern::pVecDiv2x2: {
        if (newton3) {
          SoAFunctorPairImpl<true, sorted, VectorizationPattern::pVecDiv2x2>(soa1, soa2, sortingData);
        } else {
          SoAFunctorPairImpl<false, sorted, VectorizationPattern::pVecDiv2x2>(soa1, soa2, sortingData);
        }
        break;
      }
      case VectorizationPattern::pVecx1: {
        if (newton3) {
          SoAFunctorPairImpl<true, sorted, VectorizationPattern::pVecx1>(soa1, soa2, sortingData);
        } else {
          SoAFunctorPairImpl<false, sorted, VectorizationPattern::pVecx1>(soa1, soa2, sortingData);
        }
        break;
      }
      default:
        autopas::utils::ExceptionHandler::exception("Unknown VectorizationPattern!");
    }
  }

  /**
   * Templatized implementation of SoAFunctorPair.
   *
   * When sorted=false the original SoA arrays are iterated directly.
   * When sorted=true the iteration is pruned via the supplied sortingData.
   *
   * @tparam newton3
   * @tparam sorted Whether to use pruned index data in sortingData.
   * @tparam vecPattern Vectorization pattern.
   * @param soa1
   * @param soa2
   * @param sortingData Pruning data needed for sorting optimization.
   */
  template <bool newton3, bool sorted, VectorizationPattern vecPattern>
  inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                 autopas::utils::optRef<const autopas::SoASortingData> sortingData = std::nullopt) {
    const size_t n1 = soa1.size();
    const size_t n2 = soa2.size();
    if (n1 == 0 || n2 == 0) {
      return;
    }

    const double *sqrtEpsilon1Ptr = nullptr;
    const double *halfSigma1Ptr = nullptr;
    const double *sqrtEpsilon2Ptr = nullptr;
    const double *halfSigma2Ptr = nullptr;
    if constexpr (useMixing) {
      sqrtEpsilon1Ptr = soa1.template begin<AttrSqrt24Epsilon>();
      halfSigma1Ptr = soa1.template begin<AttrHalfSigma>();
      sqrtEpsilon2Ptr = soa2.template begin<AttrSqrt24Epsilon>();
      halfSigma2Ptr = soa2.template begin<AttrHalfSigma>();

      if (soaIsHomogeneous(sqrtEpsilon1Ptr, halfSigma1Ptr, n1) and
          soaIsHomogeneous(sqrtEpsilon2Ptr, halfSigma2Ptr, n2)) {
        runPairLoop<newton3, sorted, vecPattern, true>(soa1, soa2, sortingData, sqrtEpsilon1Ptr, halfSigma1Ptr,
                                                       sqrtEpsilon2Ptr, halfSigma2Ptr);
        return;
      }
      runPairLoop<newton3, sorted, vecPattern, false>(soa1, soa2, sortingData, sqrtEpsilon1Ptr, halfSigma1Ptr,
                                                      sqrtEpsilon2Ptr, halfSigma2Ptr);
    } else {
      runPairLoop<newton3, sorted, vecPattern, true>(soa1, soa2, sortingData, nullptr, nullptr, nullptr, nullptr);
    }
  }

  template <bool homogeneousParams>
  inline void runSingleLoop(const size_t soaSize, const double *const __restrict xPtr,
                            const double *const __restrict yPtr, const double *const __restrict zPtr,
                            const autopas::OwnershipState *const __restrict ownedStatePtr,
                            double *const __restrict fxPtr, double *const __restrict fyPtr,
                            double *const __restrict fzPtr, const double *const __restrict sqrtEpsilonPtr,
                            const double *const __restrict halfSigmaPtr) {
    auto virialSumX = highway::Zero(tag_double);
    auto virialSumY = highway::Zero(tag_double);
    auto virialSumZ = highway::Zero(tag_double);
    auto uPotSum = highway::Zero(tag_double);

    static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                  "OwnershipStates underlying type should be int64_t!");

    for (size_t i = soaSize - 1; i >= 1; i -= 1) {
      handleILoopBody<true, true, false, VectorizationPattern::p1xVec, homogeneousParams>(
          i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, fxPtr, fyPtr, fzPtr,
          sqrtEpsilonPtr, halfSigmaPtr, sqrtEpsilonPtr, halfSigmaPtr, virialSumX, virialSumY, virialSumZ, uPotSum, 0, 0,
          i);
    }

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  template <bool newton3, bool sorted, VectorizationPattern vecPattern, bool homogeneousParams>
  inline void runPairLoop(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                          autopas::utils::optRef<const autopas::SoASortingData> sortingData,
                          const double *const __restrict sqrtEpsilon1Ptr, const double *const __restrict halfSigma1Ptr,
                          const double *const __restrict sqrtEpsilon2Ptr,
                          const double *const __restrict halfSigma2Ptr) {
    const size_t n1 = soa1.size();
    const size_t n2 = soa2.size();

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

    const size_t startI = sorted && sortingData.has_value() ? sortingData->get().startI : 0;
    const size_t endI = sorted && sortingData.has_value() ? sortingData->get().endI : n1;

    VectorDouble virialSumX = highway::Zero(tag_double);
    VectorDouble virialSumY = highway::Zero(tag_double);
    VectorDouble virialSumZ = highway::Zero(tag_double);
    VectorDouble uPotSum = highway::Zero(tag_double);

    constexpr size_t iStep = iStepSize<vecPattern>();
    constexpr size_t jStep = jStepSize<vecPattern>();

    size_t i = startI;
    for (; i + iStep <= endI; i += iStep) {
      size_t jVecEnd{};
      size_t jVecStart = 0;
      if constexpr (sorted) {
        const auto &sd = sortingData->get();
        jVecEnd = sd.maxIndex[i + iStep - 1];
        jVecStart = sd.minIndex[i];
        if (jVecStart >= jVecEnd) {
          continue;
        }
        jVecStart = jVecStart - (jVecStart % jStep);
      } else {
        jVecEnd = n2;
      }
      handleILoopBody<false, newton3, false, vecPattern, homogeneousParams>(
          i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx1Ptr, fy1Ptr, fz1Ptr, fx2Ptr,
          fy2Ptr, fz2Ptr, sqrtEpsilon1Ptr, halfSigma1Ptr, sqrtEpsilon2Ptr, halfSigma2Ptr, virialSumX, virialSumY,
          virialSumZ, uPotSum, 0, jVecStart, jVecEnd);
    }
    if constexpr (vecPattern != VectorizationPattern::p1xVec) {
      const size_t restI = endI - i;
      if (restI > 0) {
        size_t jVecEnd = n2;
        size_t jVecStart = 0;
        if constexpr (sorted) {
          const auto &sd = sortingData->get();
          jVecEnd = sd.maxIndex[i + restI - 1];
          jVecStart = sd.minIndex[i];
          if (jVecStart < jVecEnd) {
            jVecStart = jVecStart - (jVecStart % jStep);
          }
        }
        if (jVecStart < jVecEnd) {
          handleILoopBody<false, newton3, true, vecPattern, homogeneousParams>(
              i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx1Ptr, fy1Ptr, fz1Ptr,
              fx2Ptr, fy2Ptr, fz2Ptr, sqrtEpsilon1Ptr, halfSigma1Ptr, sqrtEpsilon2Ptr, halfSigma2Ptr, virialSumX,
              virialSumY, virialSumZ, uPotSum, restI, jVecStart, jVecEnd);
        }
      }
    }

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  /**
   * Actual inner kernel of the SoAFunctors.
   */
  template <bool newton3, bool remainderJ, bool reversed, VectorizationPattern vecPattern,
            bool homogeneousParams = not useMixing>
  HWY_INLINE void SoAKernel(
      const size_t j, [[maybe_unused]] const MaskDouble &notDummyMaskI,
      [[maybe_unused]] const VectorDouble &energyFactorI, [[maybe_unused]] const VectorDouble &energyFactorIPlus1,
      const int64_t *const __restrict ownedStatePtr2, const VectorDouble &x1, const VectorDouble &y1,
      const VectorDouble &z1, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
      const double *const __restrict z2Ptr, double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
      double *const __restrict fz2Ptr, const VectorDouble &sqrtEpsilon1, const VectorDouble &sigmaHalf1,
      [[maybe_unused]] const double *const __restrict sqrtEpsilonPtr2,
      [[maybe_unused]] const double *const __restrict halfSigmaPtr2, VectorDouble &fxAcc, VectorDouble &fyAcc,
      VectorDouble &fzAcc, [[maybe_unused]] VectorDouble &virialSumX, [[maybe_unused]] VectorDouble &virialSumY,
      [[maybe_unused]] VectorDouble &virialSumZ, [[maybe_unused]] VectorDouble &uPotSum, const size_t restJ) const {
    const VectorDouble x2 = loadRegister<false, remainderJ, reversed, vecPattern>(tag_double, x2Ptr, j, restJ);
    const VectorDouble y2 = loadRegister<false, remainderJ, reversed, vecPattern>(tag_double, y2Ptr, j, restJ);
    const VectorDouble z2 = loadRegister<false, remainderJ, reversed, vecPattern>(tag_double, z2Ptr, j, restJ);

    // distance calculations using FMA
    const auto drX = highway::Sub(x1, x2);
    const auto drY = highway::Sub(y1, y2);
    const auto drZ = highway::Sub(z1, z2);
    const auto dr2 = highway::MulAdd(drZ, drZ, highway::MulAdd(drY, drY, highway::Mul(drX, drX)));

    const VectorDouble cutoffSquaredVec = highway::Set(tag_double, _cutoffSquareAoS);

    MaskDouble cutoffMaskI;
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      cutoffMaskI = highway::Le(dr2, cutoffSquaredVec);
    } else {
      cutoffMaskI = highway::MaskedLe(notDummyMaskI, dr2, cutoffSquaredVec);
    }

    if (highway::AllFalse(tag_double, cutoffMaskI)) {
      return;
    }

    const VectorLong ownedStateJ =
        loadRegister<false, remainderJ, reversed, vecPattern>(tag_long, ownedStatePtr2, j, restJ);
    const auto notDummyMaskJ = highway::RebindMask(tag_double, highway::Ne(ownedStateJ, highway::Zero(tag_long)));
    const auto cutoffDummyMask = highway::And(cutoffMaskI, notDummyMaskJ);

    if (highway::AllFalse(tag_double, cutoffDummyMask)) [[unlikely]] {
      return;
    }

    VectorDouble epsilon24s;
    VectorDouble sigmaSquareds;

    if constexpr (homogeneousParams) {
      epsilon24s = sqrtEpsilon1;
      sigmaSquareds = sigmaHalf1;
    } else {
      const VectorDouble sqrtEpsilon2 =
          loadRegister<false, remainderJ, reversed, vecPattern>(tag_double, sqrtEpsilonPtr2, j, restJ);
      const VectorDouble sigmaHalf2 =
          loadRegister<false, remainderJ, reversed, vecPattern>(tag_double, halfSigmaPtr2, j, restJ);
      epsilon24s = highway::Mul(sqrtEpsilon1, sqrtEpsilon2);
      const auto sigma = highway::Add(sigmaHalf1, sigmaHalf2);
      sigmaSquareds = highway::Mul(sigma, sigma);
    }

    // compute LJ Potential with overlapped dependency chains and FMA
    const auto one = highway::Set(tag_double, 1.0);
    const auto invDr2 = highway::Div(one, dr2);
    const auto lj2 = highway::Mul(sigmaSquareds, invDr2);
    const auto eps24InvDr2 = highway::Mul(epsilon24s, invDr2);
    const auto lj4 = highway::Mul(lj2, lj2);
    const auto lj6 = highway::Mul(lj2, lj4);
    const auto twoLj6m1 = highway::MulSub(lj6, highway::Set(tag_double, 2.0), one);
    const auto lj6Eps24InvDr2 = highway::Mul(lj6, eps24InvDr2);
    const auto fac = highway::Mul(twoLj6m1, lj6Eps24InvDr2);

    const auto facMasked = highway::IfThenElseZero(cutoffDummyMask, fac);

    if constexpr (not newton3 and not calculateGlobals) {
      fxAcc = highway::MulAdd(drX, facMasked, fxAcc);
      fyAcc = highway::MulAdd(drY, facMasked, fyAcc);
      fzAcc = highway::MulAdd(drZ, facMasked, fzAcc);
    } else {
      VectorDouble fx = highway::Mul(drX, facMasked);
      VectorDouble fy = highway::Mul(drY, facMasked);
      VectorDouble fz = highway::Mul(drZ, facMasked);

      fxAcc = highway::Add(fxAcc, fx);
      fyAcc = highway::Add(fyAcc, fy);
      fzAcc = highway::Add(fzAcc, fz);

      if constexpr (calculateGlobals) {
        VectorDouble shift6s;
        if constexpr (applyShift) {
          if constexpr (not useMixing) {
            shift6s = highway::Set(tag_double, _shift6AoS);
          } else {
            const VectorDouble invCutoffSquared = highway::Set(tag_double, _invCutoffSquareAoS);
            const auto tmpShift1 = highway::Mul(sigmaSquareds, invCutoffSquared);
            const auto tmpShift2 = highway::Mul(highway::Mul(tmpShift1, tmpShift1), tmpShift1);
            shift6s = highway::Mul(epsilon24s, highway::Sub(tmpShift2, highway::Mul(tmpShift2, tmpShift2)));
          }
        } else {
          shift6s = highway::Zero(tag_double);
        }

        const auto virialX = highway::Mul(fx, drX);
        const auto virialY = highway::Mul(fy, drY);
        const auto virialZ = highway::Mul(fz, drZ);

        const auto lj12m6 = highway::MulSub(lj6, lj6, lj6);
        const auto uPot = highway::MulAdd(epsilon24s, lj12m6, shift6s);
        const auto uPotMasked = highway::IfThenElseZero(cutoffDummyMask, uPot);

        VectorDouble energyFactor;
        if constexpr (newton3) {
          const auto owned = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
          const auto ownedMaskJ = highway::RebindMask(tag_double, highway::Eq(ownedStateJ, owned));
          energyFactor = highway::IfThenElse(ownedMaskJ, energyFactorIPlus1, energyFactorI);
        } else {
          energyFactor = energyFactorI;
        }

        uPotSum = highway::MulAdd(energyFactor, uPotMasked, uPotSum);
        virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
        virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
        virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
      }

      if constexpr (newton3) {
        handleNewton3Reduction<remainderJ, vecPattern>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, restJ);
      }
    }
  }

  template <bool newton3, bool remainder = false, bool homogeneousParams = not useMixing>
  HWY_INLINE void SoAKernelVerlet(
      const size_t j, const VectorDouble &cutoffSquared, [[maybe_unused]] const VectorDouble &shiftOrInvCutoffSquared,
      [[maybe_unused]] const VectorDouble &energyFactorI, const VectorDouble &x1, const VectorDouble &y1,
      const VectorDouble &z1, const double *const __restrict xPtr, const double *const __restrict yPtr,
      const double *const __restrict zPtr, const int64_t *const __restrict ownedStatePtr,
      double *const __restrict fxPtr, double *const __restrict fyPtr, double *const __restrict fzPtr,
      [[maybe_unused]] const VectorDouble &sqrtEpsilon1, [[maybe_unused]] const VectorDouble &sigmaHalf1,
      [[maybe_unused]] const double *const __restrict sqrtEpsilonPtr,
      [[maybe_unused]] const double *const __restrict sigmaHalfPtr, const size_t *const __restrict neighborList,
      VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc, [[maybe_unused]] VectorDouble &virialSumX,
      [[maybe_unused]] VectorDouble &virialSumY, [[maybe_unused]] VectorDouble &virialSumZ,
      [[maybe_unused]] VectorDouble &uPotSum, [[maybe_unused]] const size_t rest = 0) const {
    [[maybe_unused]] MaskDouble restMaskDouble;
    [[maybe_unused]] MaskLong restMaskLong;
    VectorLong indices;

    if constexpr (remainder) {
      restMaskDouble = highway::FirstN(tag_double, rest);
      restMaskLong = highway::FirstN(tag_long, rest);
      indices = highway::LoadN(tag_long, reinterpret_cast<const int64_t *>(neighborList + j), rest);
    } else {
      indices = highway::LoadU(tag_long, reinterpret_cast<const int64_t *>(neighborList + j));
    }

    VectorDouble x2;
    VectorDouble y2;
    VectorDouble z2;
    VectorLong ownedState2;
    VectorDouble sqrtEpsilon2 = highway::Zero(tag_double);
    VectorDouble sigmaHalf2 = highway::Zero(tag_double);

    if constexpr (remainder) {
      x2 = highway::MaskedGatherIndex(restMaskDouble, tag_double, xPtr, indices);
      y2 = highway::MaskedGatherIndex(restMaskDouble, tag_double, yPtr, indices);
      z2 = highway::MaskedGatherIndex(restMaskDouble, tag_double, zPtr, indices);
      ownedState2 = highway::MaskedGatherIndex(restMaskLong, tag_long, ownedStatePtr, indices);
      if constexpr (useMixing and not homogeneousParams) {
        sqrtEpsilon2 = highway::MaskedGatherIndex(restMaskDouble, tag_double, sqrtEpsilonPtr, indices);
        sigmaHalf2 = highway::MaskedGatherIndex(restMaskDouble, tag_double, sigmaHalfPtr, indices);
      }
    } else {
      x2 = highway::GatherIndex(tag_double, xPtr, indices);
      y2 = highway::GatherIndex(tag_double, yPtr, indices);
      z2 = highway::GatherIndex(tag_double, zPtr, indices);
      ownedState2 = highway::GatherIndex(tag_long, ownedStatePtr, indices);
      if constexpr (useMixing and not homogeneousParams) {
        sigmaHalf2 = highway::GatherIndex(tag_double, sigmaHalfPtr, indices);
        HWY_ALIGN std::array<double, _maxVecLengthDouble> sqrtEpsArr;
        for (size_t k = 0; k < _vecLengthDouble; ++k) {
          sqrtEpsArr[k] = sqrtEpsilonPtr[neighborList[j + k]];
        }
        sqrtEpsilon2 = highway::Load(tag_double, sqrtEpsArr.data());
      }
    }

    const auto drX = highway::Sub(x1, x2);
    const auto drY = highway::Sub(y1, y2);
    const auto drZ = highway::Sub(z1, z2);
    const auto dr2 = highway::MulAdd(drZ, drZ, highway::MulAdd(drY, drY, highway::Mul(drX, drX)));

    VectorDouble epsilon24s;
    VectorDouble sigmaSquareds;

    if constexpr (homogeneousParams) {
      epsilon24s = sqrtEpsilon1;
      sigmaSquareds = sigmaHalf1;
    } else {
      epsilon24s = highway::Mul(sqrtEpsilon1, sqrtEpsilon2);
      const auto sigma = highway::Add(sigmaHalf1, sigmaHalf2);
      sigmaSquareds = highway::Mul(sigma, sigma);
    }

    const auto notDummyMask2 = highway::RebindMask(tag_double, highway::Ne(ownedState2, highway::Zero(tag_long)));
    const auto cutoffDummyMask = highway::MaskedLe(notDummyMask2, dr2, cutoffSquared);

    if constexpr (remainder or not useMixing) {
      if (highway::AllFalse(tag_double, cutoffDummyMask)) {
        return;
      }
    }

    const auto one = highway::Set(tag_double, 1.0);
    VectorDouble invDr2 = highway::Div(one, dr2);

    [[maybe_unused]] VectorDouble fx2Old;
    [[maybe_unused]] VectorDouble fy2Old;
    [[maybe_unused]] VectorDouble fz2Old;
    if constexpr (newton3 and not remainder and (homogeneousParams or not calculateGlobals)) {
      fx2Old = highway::GatherIndex(tag_double, fxPtr, indices);
      fy2Old = highway::GatherIndex(tag_double, fyPtr, indices);
      fz2Old = highway::GatherIndex(tag_double, fzPtr, indices);
      if constexpr (homogeneousParams) {
        asm("" : "+v"(invDr2), "+v"(fx2Old), "+v"(fy2Old), "+v"(fz2Old));
      }
    }

    const auto lj2 = highway::Mul(sigmaSquareds, invDr2);
    const auto eps24InvDr2 = highway::Mul(epsilon24s, invDr2);
    const auto lj4 = highway::Mul(lj2, lj2);
    const auto lj6 = highway::Mul(lj2, lj4);
    const auto twoLj6m1 = highway::MulSub(lj6, highway::Set(tag_double, 2.0), one);
    const auto lj6Eps24InvDr2 = highway::Mul(lj6, eps24InvDr2);
    const auto fac = highway::Mul(twoLj6m1, lj6Eps24InvDr2);

    const auto facMasked = highway::IfThenElseZero(cutoffDummyMask, fac);

    if constexpr (not newton3 and not calculateGlobals) {
      fxAcc = highway::MulAdd(drX, facMasked, fxAcc);
      fyAcc = highway::MulAdd(drY, facMasked, fyAcc);
      fzAcc = highway::MulAdd(drZ, facMasked, fzAcc);
    } else {
      const VectorDouble fx = highway::Mul(drX, facMasked);
      const VectorDouble fy = highway::Mul(drY, facMasked);
      const VectorDouble fz = highway::Mul(drZ, facMasked);

      fxAcc = highway::Add(fxAcc, fx);
      fyAcc = highway::Add(fyAcc, fy);
      fzAcc = highway::Add(fzAcc, fz);

      if constexpr (newton3) {
        if constexpr (remainder) {
          fx2Old = highway::MaskedGatherIndex(restMaskDouble, tag_double, fxPtr, indices);
          fy2Old = highway::MaskedGatherIndex(restMaskDouble, tag_double, fyPtr, indices);
          fz2Old = highway::MaskedGatherIndex(restMaskDouble, tag_double, fzPtr, indices);
        } else if constexpr (not homogeneousParams and calculateGlobals) {
          fx2Old = highway::GatherIndex(tag_double, fxPtr, indices);
          fy2Old = highway::GatherIndex(tag_double, fyPtr, indices);
          fz2Old = highway::GatherIndex(tag_double, fzPtr, indices);
        }

        const auto fx2New = highway::Sub(fx2Old, fx);
        const auto fy2New = highway::Sub(fy2Old, fy);
        const auto fz2New = highway::Sub(fz2Old, fz);

        if constexpr (remainder) {
          highway::MaskedScatterIndex(fx2New, restMaskDouble, tag_double, fxPtr, indices);
          highway::MaskedScatterIndex(fy2New, restMaskDouble, tag_double, fyPtr, indices);
          highway::MaskedScatterIndex(fz2New, restMaskDouble, tag_double, fzPtr, indices);
        } else {
          highway::ScatterIndex(fx2New, tag_double, fxPtr, indices);
          highway::ScatterIndex(fy2New, tag_double, fyPtr, indices);
          highway::ScatterIndex(fz2New, tag_double, fzPtr, indices);
        }
      }

      if constexpr (calculateGlobals) {
        const auto virialX = highway::Mul(fx, drX);
        const auto virialY = highway::Mul(fy, drY);
        const auto virialZ = highway::Mul(fz, drZ);

        const auto lj12m6 = highway::MulSub(lj6, lj6, lj6);
        VectorDouble uPot;
        if constexpr (applyShift) {
          VectorDouble shift6s;
          if constexpr (useMixing) {
            const auto tmpShift1 = highway::Mul(sigmaSquareds, shiftOrInvCutoffSquared);
            const auto tmpShift2 = highway::Mul(highway::Mul(tmpShift1, tmpShift1), tmpShift1);
            shift6s = highway::Mul(epsilon24s, highway::Sub(tmpShift2, highway::Mul(tmpShift2, tmpShift2)));
          } else {
            shift6s = shiftOrInvCutoffSquared;
          }
          uPot = highway::MulAdd(epsilon24s, lj12m6, shift6s);
        } else {
          uPot = highway::Mul(epsilon24s, lj12m6);
        }
        const auto uPotMasked = highway::IfThenElseZero(cutoffDummyMask, uPot);

        VectorDouble energyFactor;
        if constexpr (newton3) {
          const auto owned = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
          const auto ownedMaskJ = highway::RebindMask(tag_double, highway::Eq(ownedState2, owned));
          energyFactor = highway::MaskedAddOr(energyFactorI, ownedMaskJ, energyFactorI, one);
        } else {
          energyFactor = energyFactorI;
        }

        uPotSum = highway::MulAdd(energyFactor, uPotMasked, uPotSum);
        virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
        virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
        virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
      }
    }
  }

 public:
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorVerlet()
   */
  inline void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                               const std::span<const size_t> neighborList, const bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

 private:
  template <bool newton3>
  HWY_NOINLINE void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                                         std::span<const size_t> neighborList) {
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto stateI = ownedStatePtr[indexFirst];
    if (stateI == autopas::OwnershipState::dummy) {
      return;
    }

    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    const double *sqrtEpsilonPtr = nullptr;
    const double *halfSigmaPtr = nullptr;
    if constexpr (useMixing) {
      sqrtEpsilonPtr = soa.template begin<AttrSqrt24Epsilon>();
      halfSigmaPtr = soa.template begin<AttrHalfSigma>();
    }

    const size_t neighborListSize = neighborList.size();
    const size_t rest = neighborListSize & (_vecLengthDouble - 1);
    const size_t vecEnd = neighborListSize - rest;
    const size_t *const __restrict nlData = neighborList.data();
    size_t j = 0;

    if (vecEnd > 0) {
      VectorDouble sqrtEpsilon1;
      VectorDouble sigmaHalf1;
      if constexpr (useMixing) {
        sqrtEpsilon1 = highway::Set(tag_double, sqrtEpsilonPtr[indexFirst]);
        sigmaHalf1 = highway::Set(tag_double, halfSigmaPtr[indexFirst]);
      } else {
        sqrtEpsilon1 = highway::Set(tag_double, _epsilon24AoS);
        sigmaHalf1 = highway::Set(tag_double, _sigmaSquareAoS);
      }

      const VectorDouble cutoffSquared = highway::Set(tag_double, _cutoffSquareAoS);
      [[maybe_unused]] VectorDouble shiftOrInvCutoffSquared = highway::Zero(tag_double);
      if constexpr (calculateGlobals and applyShift) {
        shiftOrInvCutoffSquared = highway::Set(tag_double, useMixing ? _invCutoffSquareAoS : _shift6AoS);
      }

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

      [[maybe_unused]] VectorDouble energyFactorI = highway::Zero(tag_double);
      if constexpr (calculateGlobals) {
        const double baseFactor = (stateI == autopas::OwnershipState::owned) ? 1.0 : 0.0;
        energyFactorI = highway::Set(tag_double, baseFactor);
      }

      for (; j < vecEnd; j += _vecLengthDouble) {
        SoAKernelVerlet<newton3, false>(j, cutoffSquared, shiftOrInvCutoffSquared, energyFactorI, x1, y1, z1, xPtr,
                                        yPtr, zPtr, reinterpret_cast<const int64_t *>(ownedStatePtr), fxPtr, fyPtr,
                                        fzPtr, sqrtEpsilon1, sigmaHalf1, sqrtEpsilonPtr, halfSigmaPtr, nlData, fxAcc,
                                        fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, uPotSum);
      }

      reduceAdd3_InterleaveFirst(tag_double, fxAcc, fyAcc, fzAcc, fxPtr[indexFirst], fyPtr[indexFirst],
                                 fzPtr[indexFirst]);

      if constexpr (calculateGlobals) {
        computeGlobals(virialSumX, virialSumY, virialSumZ, uPotSum);
      }
    }

    if (j < neighborListSize) {
      const double xi = xPtr[indexFirst];
      const double yi = yPtr[indexFirst];
      const double zi = zPtr[indexFirst];
      double sqrtEpsI = 0.;
      double halfSigI = 0.;
      if constexpr (useMixing) {
        sqrtEpsI = sqrtEpsilonPtr[indexFirst];
        halfSigI = halfSigmaPtr[indexFirst];
      }
      [[maybe_unused]] const double baseFactor =
          (calculateGlobals and stateI == autopas::OwnershipState::owned) ? 1.0 : 0.0;
      const double cutoffSq = _cutoffSquareAoS;
      double fxi = 0., fyi = 0., fzi = 0.;
      [[maybe_unused]] double virX = 0., virY = 0., virZ = 0., pot6Sum = 0.;
      for (; j < neighborListSize; ++j) {
        const size_t idxJ = nlData[j];
        const auto stateJ = ownedStatePtr[idxJ];
        if (stateJ == autopas::OwnershipState::dummy) {
          continue;
        }
        const double dx = xi - xPtr[idxJ];
        const double dy = yi - yPtr[idxJ];
        const double dz = zi - zPtr[idxJ];
        const double dr2 = dx * dx + dy * dy + dz * dz;
        if (dr2 <= cutoffSq) {
          double fx = 0., fy = 0., fz = 0., pot6 = 0.;
          computeScalarLJInteraction(dr2, dx, dy, dz, sqrtEpsI, halfSigI, useMixing ? sqrtEpsilonPtr[idxJ] : 0.,
                                     useMixing ? halfSigmaPtr[idxJ] : 0., fx, fy, fz, pot6);
          fxi += fx;
          fyi += fy;
          fzi += fz;
          if constexpr (newton3) {
            fxPtr[idxJ] -= fx;
            fyPtr[idxJ] -= fy;
            fzPtr[idxJ] -= fz;
          }
          if constexpr (calculateGlobals) {
            const double factor = baseFactor + ((newton3 and stateJ == autopas::OwnershipState::owned) ? 1.0 : 0.0);
            pot6Sum += factor * pot6;
            virX += factor * fx * dx;
            virY += factor * fy * dy;
            virZ += factor * fz * dz;
          }
        }
      }
      fxPtr[indexFirst] += fxi;
      fyPtr[indexFirst] += fyi;
      fzPtr[indexFirst] += fzi;
      if constexpr (calculateGlobals) {
        const int threadID = autopas::autopas_get_thread_num();
        _aosThreadData[threadID].potentialEnergySum += pot6Sum;
        _aosThreadData[threadID].virialSum[0] += virX;
        _aosThreadData[threadID].virialSum[1] += virY;
        _aosThreadData[threadID].virialSum[2] += virZ;
      }
    }
  }

 public:
  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    if constexpr (useMixing) {
      return std::array<typename Particle_T::AttributeNames, 10>{Particle_T::AttributeNames::id,
                                                                 Particle_T::AttributeNames::posX,
                                                                 Particle_T::AttributeNames::posY,
                                                                 Particle_T::AttributeNames::posZ,
                                                                 Particle_T::AttributeNames::forceX,
                                                                 Particle_T::AttributeNames::forceY,
                                                                 Particle_T::AttributeNames::forceZ,
                                                                 AttrSqrt24Epsilon,
                                                                 AttrHalfSigma,
                                                                 Particle_T::AttributeNames::ownershipState};
    } else {
      return std::array<typename Particle_T::AttributeNames, 8>{
          Particle_T::AttributeNames::id,     Particle_T::AttributeNames::posX,
          Particle_T::AttributeNames::posY,   Particle_T::AttributeNames::posZ,
          Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY,
          Particle_T::AttributeNames::forceZ, Particle_T::AttributeNames::ownershipState};
    }
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    if constexpr (useMixing) {
      return std::array<typename Particle_T::AttributeNames, 7>{Particle_T::AttributeNames::id,
                                                                Particle_T::AttributeNames::posX,
                                                                Particle_T::AttributeNames::posY,
                                                                Particle_T::AttributeNames::posZ,
                                                                AttrSqrt24Epsilon,
                                                                AttrHalfSigma,
                                                                Particle_T::AttributeNames::ownershipState};
    } else {
      return std::array<typename Particle_T::AttributeNames, 5>{
          Particle_T::AttributeNames::id, Particle_T::AttributeNames::posX, Particle_T::AttributeNames::posY,
          Particle_T::AttributeNames::posZ, Particle_T::AttributeNames::ownershipState};
    }
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle_T::AttributeNames, 3>{
        Particle_T::AttributeNames::forceX, Particle_T::AttributeNames::forceY, Particle_T::AttributeNames::forceZ};
  }

  /**
   * @copydoc autopas::Functor::supportsSoASorting
   */
  static constexpr bool supportsSoASorting = true;

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

    if constexpr (calculateGlobals) {
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
    if constexpr (not calculateGlobals) {
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
    if constexpr (not calculateGlobals) {
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

  // cutoff squared used in the AoS functor.
  const double _cutoffSquareAoS{0.};
  const double _invCutoffSquareAoS{0.};
  // epsilon, sigma and shift6 used in the AoS functor.
  double _epsilon24AoS{0.}, _sigmaSquareAoS{0.}, _shift6AoS{0.};

  // optional to hold a reference to the ParticlePropertiesLibrary. If a ParticlePropertiesLibrary is not used the
  // optional is empty. Kept only to validate the useMixing/PPL-presence invariant in the constructor (for
  // API/behavioral compatibility with other functors); the mixing math itself no longer looks anything up here:
  // epsilon and sigma are read directly from the interacting particles/SoA columns and mixed inline.
  std::optional<std::reference_wrapper<ParticlePropertiesLibrary<double, size_t>>> _PPLibrary;

  // accumulators for the globals (potential energy and virial).
  double _potentialEnergySum{0.};
  std::array<double, 3> _virialSum{0., 0., 0.};
  std::vector<AoSThreadData> _aosThreadData{};

  // flag to indicate whether post-processing has been performed.
  bool _postProcessed{false};

  // The Vectorization Pattern currently used in the SoA functor.
  VectorizationPattern _vecPattern{};

  // Vectorization Pattern that the functor can handle.
  static constexpr std::array<VectorizationPattern, 4> _vecPatternsAllowed = {
      VectorizationPattern::p1xVec, VectorizationPattern::p2xVecDiv2, VectorizationPattern::pVecDiv2x2,
      VectorizationPattern::pVecx1};
};
}  // namespace mdLib
