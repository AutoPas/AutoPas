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
#include <cmath>
#include <functional>
#include <optional>
#include <vector>

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
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
 * A functor to handle Lennard-Jones interactions between two particles (molecules).
 * This functor uses the SIMD abstraction library Google Highway to provide architecture independent vectorization.
 *
 * Key Architectural Principles:
 * 1. Zero Particle Bloat: Particles only store standard attributes (positions, typeId, ownershipState).
 * 2. 1D On-The-Fly Parameter Projection: For separable Lorentz-Berthelot mixing, per-type parameters
 *    (halfSigma, sqrtEpsilon) are projected into a contiguous 64-byte aligned buffer at cell entry.
 * 3. Multi-Regime Hardware Specialization:
 *    - Regime 1: Homogeneous Single-Type (Broadcast once into vector registers, zero parameter loads).
 *    - Regime 2: Small-N Direct Scalar Shortcut (N <= 8 in SoAFunctorSingle, avoids SIMD masking & reductions).
 *    - Regime 3: 1D Contiguous SIMD Streaming (N > 8, multi-pattern contiguous loads + vector add/mul, 0 STLF stalls).
 *
 * @tparam Particle_T The type of particle
 * @tparam applyShift Switch for the LJ potential to be truncated shifted
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate. Not implemented for this functor. Please use the AutoVec functor.
 * @tparam relevantForTuning Whether the auto-tuner should consider this functor.
 */
template <class Particle_T, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctorHWY
    : public autopas::PairwiseFunctor<Particle_T, LJFunctorHWY<Particle_T, applyShift, useMixing, useNewton3,
                                                               calculateGlobals, countFLOPs, relevantForTuning>> {
  using SoAArraysType = typename Particle_T::SoAArraysType;

  template <size_t StackCapacity = 512>
  struct MixingBuffers {
    // Aligned buffers to make sure we can use hwy::Load instead of hwy::LoadU
    alignas(64) std::array<double, StackCapacity> halfSigmaStack;
    alignas(64) std::array<double, StackCapacity> sqrtEpsStack;
    std::vector<double, autopas::AlignedAllocator<double, 64>> halfSigmaHeap;
    std::vector<double, autopas::AlignedAllocator<double, 64>> sqrtEpsHeap;

    void getMixingPointers(const size_t *__restrict typeIDs, const size_t soaSize,
                           const std::vector<double> &halfSigmaPerType, const std::vector<double> &sqrtEpsilonPerType,
                           double *&halfSigmaPtr, double *&sqrtEpsilonPtr) {
      if (soaSize <= StackCapacity) {
        halfSigmaPtr = halfSigmaStack.data();
        sqrtEpsilonPtr = sqrtEpsStack.data();
      } else {
        halfSigmaHeap.resize(soaSize);
        sqrtEpsHeap.resize(soaSize);
        halfSigmaPtr = halfSigmaHeap.data();
        sqrtEpsilonPtr = sqrtEpsHeap.data();
      }
#pragma omp simd
      for (size_t k = 0; k < soaSize; ++k) {
        const auto typeID = typeIDs[k];
        halfSigmaPtr[k] = halfSigmaPerType[typeID];
        sqrtEpsilonPtr[k] = sqrtEpsilonPerType[typeID];
      }
    }
  };

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
        _cutoffSquareAoS{cutoff * cutoff},
        _invCutoffSquareAoS(1. / (cutoff * cutoff)),
        _cutoffSquaredVec{highway::Set(tag_double, cutoff * cutoff)} {
    if constexpr (calculateGlobals) {
      _aosThreadData.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      AutoPasLog(DEBUG, "Using LJFunctorHWY with countFLOPs but FLOP counting is not implemented.");
    }

    if constexpr (useMixing) {
      if (not particlePropertiesLibrary.has_value()) {
        throw std::runtime_error("Mixing is enabled but no ParticlePropertiesLibrary was provided!");
      }
      initFromPPL(particlePropertiesLibrary.value().get());
    } else {
      if (particlePropertiesLibrary.has_value()) {
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
    auto sigmaSquare = _sigmaSquareAoS;
    auto epsilon24 = _epsilon24AoS;
    auto shift6 = _shift6AoS;
    if constexpr (useMixing) {
      const auto typeI = i.getTypeId();
      const auto typeJ = j.getTypeId();
      const auto sigma = _halfSigmaPerType[typeI] + _halfSigmaPerType[typeJ];
      sigmaSquare = sigma * sigma;
      epsilon24 = 24.0 * _sqrtEpsilonPerType[typeI] * _sqrtEpsilonPerType[typeJ];
      if constexpr (applyShift) {
        shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, _cutoffSquareAoS);
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
      j.subF(f);
    }
    if constexpr (calculateGlobals) {
      const auto virial = dr * f;
      const double potentialEnergy6 = epsilon24 * lj12m6 + shift6;

      const int threadnum = autopas::autopas_get_thread_num();

      if (i.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadData[threadnum].virialSum += virial;
      }
      if (newton3 and j.isOwned()) {
        _aosThreadData[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadData[threadnum].virialSum += virial;
      }
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always do a newton3 like traversal of the SoA using only the 1xVec Pattern.
   */
  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, const bool newton3) final {
    const auto soaSize = soa.size();
    if (soaSize <= 1) return;

    // obtain pointers for particle attributes
    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    const auto *const __restrict typeIDPtr = soa.template begin<Particle_T::AttributeNames::typeId>();

    // -------------------------------------------------------------------------
    // Scalar computations for low particle counts
    // -------------------------------------------------------------------------
    if (soaSize <= 7) {
      const int threadID = autopas::autopas_get_thread_num();
      for (size_t i = 1; i < soaSize; ++i) {
        if (ownedStatePtr[i] == autopas::OwnershipState::dummy) continue;
        const double xi = xPtr[i], yi = yPtr[i], zi = zPtr[i];
        double halfSigI = 0., sqrtEpsI = 0.;
        if constexpr (useMixing) {
          const size_t ti = typeIDPtr[i];
          halfSigI = _halfSigmaPerType[ti];
          sqrtEpsI = _sqrtEpsilonPerType[ti];
        }
        double fxi = 0., fyi = 0., fzi = 0.;
        for (size_t j = 0; j < i; ++j) {
          if (ownedStatePtr[j] == autopas::OwnershipState::dummy) continue;
          const double dx = xi - xPtr[j];
          const double dy = yi - yPtr[j];
          const double dz = zi - zPtr[j];
          const double dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 <= _cutoffSquareAoS) {
            double sigSq = _sigmaSquareAoS;
            double eps24 = _epsilon24AoS;
            double shift6 = _shift6AoS;
            if constexpr (useMixing) {
              const size_t tj = typeIDPtr[j];
              const double sig = halfSigI + _halfSigmaPerType[tj];
              sigSq = sig * sig;
              eps24 = 24.0 * sqrtEpsI * _sqrtEpsilonPerType[tj];
              if constexpr (applyShift) {
                shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(eps24, sigSq, _cutoffSquareAoS);
              }
            }
            const double invdr2 = 1.0 / dr2;
            const double lj2 = sigSq * invdr2;
            const double lj4 = lj2 * lj2;
            const double lj6 = lj4 * lj2;
            const double lj12 = lj6 * lj6;
            const double lj12m6 = lj12 - lj6;
            const double fac = eps24 * (lj12 + lj12m6) * invdr2;
            const double fx = dx * fac, fy = dy * fac, fz = dz * fac;
            fxi += fx;
            fyi += fy;
            fzi += fz;
            fxPtr[j] -= fx;
            fyPtr[j] -= fy;
            fzPtr[j] -= fz;
            if constexpr (calculateGlobals) {
              const double pot6 = eps24 * lj12m6 + shift6;
              if (ownedStatePtr[i] == autopas::OwnershipState::owned) {
                _aosThreadData[threadID].potentialEnergySum += pot6;
                _aosThreadData[threadID].virialSum[0] += fx * dx;
                _aosThreadData[threadID].virialSum[1] += fy * dy;
                _aosThreadData[threadID].virialSum[2] += fz * dz;
              }
              if (ownedStatePtr[j] == autopas::OwnershipState::owned) {
                _aosThreadData[threadID].potentialEnergySum += pot6;
                _aosThreadData[threadID].virialSum[0] += fx * dx;
                _aosThreadData[threadID].virialSum[1] += fy * dy;
                _aosThreadData[threadID].virialSum[2] += fz * dz;
              }
            }
          }
        }
        fxPtr[i] += fxi;
        fyPtr[i] += fyi;
        fzPtr[i] += fzi;
      }
      return;
    }

    // initialize vector accumulators
    auto virialSumX = highway::Zero(tag_double);
    auto virialSumY = highway::Zero(tag_double);
    auto virialSumZ = highway::Zero(tag_double);
    auto uPotSum = highway::Zero(tag_double);

    // Check if all particles in the cell are of the same type
    const bool isHomogeneous = soaIsHomogeneous(typeIDPtr, soaSize);

    auto runSingleLoop = [&](auto &ljParams) {
      for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(soaSize) - 1; i > 0; i -= 1) {
        handleILoopBody<true, true, false, VectorizationPattern::p1xVec>(
            i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, fxPtr, fyPtr,
            fzPtr, ljParams, virialSumX, virialSumY, virialSumZ, uPotSum, 0, 0, i);
      }
    };

    // -------------------------------------------------------------------------
    // Fast path for a homogeneous cell
    // -------------------------------------------------------------------------
    if (isHomogeneous) {
      const size_t firstType = useMixing ? typeIDPtr[0] : 0;
      double sigmaSquare = _sigmaSquareAoS;
      double epsilon24 = _epsilon24AoS;
      double shift6 = _shift6AoS;
      if constexpr (useMixing) {
        const double sigma = 2.0 * _halfSigmaPerType[firstType];
        sigmaSquare = sigma * sigma;
        epsilon24 = 24.0 * _sqrtEpsilonPerType[firstType] * _sqrtEpsilonPerType[firstType];
        if constexpr (applyShift) {
          shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquare, _cutoffSquareAoS);
        }
      }
      HomogeneousParameters homogeneousParams{sigmaSquare, epsilon24, shift6};
      runSingleLoop(homogeneousParams);
    } else {
      // -------------------------------------------------------------------------
      // Cell consists of multiple particle types
      // -------------------------------------------------------------------------
      MixingBuffers<512> mixingBuffers;
      double *halfSigPtr = nullptr;
      double *sqrtEpsPtr = nullptr;
      mixingBuffers.getMixingPointers(typeIDPtr, soaSize, _halfSigmaPerType, _sqrtEpsilonPerType, halfSigPtr,
                                      sqrtEpsPtr);

      HeterogeneousParameters heterogeneousParams{halfSigPtr, sqrtEpsPtr, halfSigPtr, sqrtEpsPtr, _invCutoffSquareAoS};
      runSingleLoop(heterogeneousParams);
    }

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
   */
  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                             const bool newton3) final {
    if (newton3) {
      SoAFunctorPairDispatch<true, false>(soa1, soa2);
    } else {
      SoAFunctorPairDispatch<false, false>(soa1, soa2);
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPairSorted()
   */
  inline void SoAFunctorPairSorted(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                   const autopas::SoASortingData &sortingData, const bool newton3) final {
    if (newton3) {
      SoAFunctorPairDispatch<true, true>(soa1, soa2, sortingData);
    } else {
      SoAFunctorPairDispatch<false, true>(soa1, soa2, sortingData);
    }
  }

 private:
  template <bool newton3, bool sorted>
  inline void SoAFunctorPairDispatch(
      autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
      const autopas::utils::optRef<const autopas::SoASortingData> sortingData = std::nullopt) {
    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }

    switch (_vecPattern) {
      case VectorizationPattern::p1xVec: {
        SoAFunctorPairImpl<newton3, sorted, VectorizationPattern::p1xVec>(soa1, soa2, sortingData);
        break;
      }
      case VectorizationPattern::p2xVecDiv2: {
        SoAFunctorPairImpl<newton3, sorted, VectorizationPattern::p2xVecDiv2>(soa1, soa2, sortingData);
        break;
      }
      case VectorizationPattern::pVecDiv2x2: {
        SoAFunctorPairImpl<newton3, sorted, VectorizationPattern::pVecDiv2x2>(soa1, soa2, sortingData);
        break;
      }
      case VectorizationPattern::pVecx1: {
        SoAFunctorPairImpl<newton3, sorted, VectorizationPattern::pVecx1>(soa1, soa2, sortingData);
        break;
      }
      default: {
        autopas::utils::ExceptionHandler::exception("Unknown VectorizationPattern!");
      }
    }
  }

  void initFromPPL(const ParticlePropertiesLibrary<double, size_t> &ppl) {
    const size_t nTypes = ppl.getNumberRegisteredSiteTypes();
    _halfSigmaPerType.resize(nTypes);
    _sqrtEpsilonPerType.resize(nTypes);
    for (size_t t = 0; t < nTypes; ++t) {
      _halfSigmaPerType[t] = ppl.getSigma(t) * 0.5;
      _sqrtEpsilonPerType[t] = std::sqrt(ppl.getEpsilon(t));
    }
    if (nTypes > 0) {
      _sigmaSquareAoS = 4.0 * _halfSigmaPerType[0] * _halfSigmaPerType[0];
      _epsilon24AoS = 24.0 * _sqrtEpsilonPerType[0] * _sqrtEpsilonPerType[0];
      if constexpr (applyShift) {
        _shift6AoS = ppl.getMixingShift6(0, 0);
      }
    }
  }

  /**
   * Checks whether all particles in an SoA buffer share the same typeId.
   */
  static inline bool soaIsHomogeneous(const size_t *const __restrict typeIDPtr, const size_t soaSize) {
    if constexpr (useMixing) {
      const size_t firstType = typeIDPtr[0];
      for (size_t k = 1; k < soaSize; ++k) {
        if (typeIDPtr[k] != firstType) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Number of i-particles processed per outer-loop step for the given vectorization pattern.
   * @tparam vecPattern
   */
  template <VectorizationPattern vecPattern>
  static size_t iStepSize() {
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
  static size_t jStepSize() {
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
   * Loads a single attribute array for the outer i-loop according to the specified VectorizationPattern.
   */
  template <bool remainder, bool reversed, VectorizationPattern vecPattern, typename Tag, typename ValueType>
  static inline auto loadIRegister(Tag tag, const ValueType *const __restrict dataPtr, const size_t i,
                                   const size_t restI) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      return highway::Set(tag, dataPtr[i]);
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      const auto lowerHalf = highway::Set(tag, dataPtr[i]);
      auto upperHalf = highway::Zero(tag);
      if constexpr (not remainder) {
        const size_t secondIndex = reversed ? i - 1 : i + 1;
        upperHalf = highway::Set(tag, dataPtr[secondIndex]);
      }
      return highway::ConcatLowerLower(tag, upperHalf, lowerHalf);
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      const size_t startIndex = reversed ? (remainder ? 0 : i - _vecLengthDouble / 2 + 1) : i;
      const size_t numLanes = remainder ? restI : _vecLengthDouble / 2;
      const auto loadedHalf = highway::LoadN(tag, &dataPtr[startIndex], numLanes);
      return highway::ConcatLowerLower(tag, loadedHalf, loadedHalf);
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      const size_t startIndex = reversed ? (remainder ? 0 : i - _vecLengthDouble + 1) : i;
      if constexpr (remainder) {
        return highway::LoadN(tag, &dataPtr[startIndex], restI);
      } else {
        return highway::LoadU(tag, &dataPtr[startIndex]);
      }
    }
  }

  /**
   * Fills the registers for the loop over i. This includes positions and ownership state.
   */
  template <bool remainder, bool reversed, VectorizationPattern vecPattern>
  static void fillIRegisters(const size_t i, const double *const __restrict xPtr, const double *const __restrict yPtr,
                             const double *const __restrict zPtr,
                             const autopas::OwnershipState *const __restrict ownedStatePtr, VectorDouble &x1,
                             VectorDouble &y1, VectorDouble &z1, MaskDouble &notDummyMaskI,
                             [[maybe_unused]] MaskDouble &ownedMaskI, const size_t restI) {
    x1 = loadIRegister<remainder, reversed, vecPattern>(tag_double, xPtr, i, restI);
    y1 = loadIRegister<remainder, reversed, vecPattern>(tag_double, yPtr, i, restI);
    z1 = loadIRegister<remainder, reversed, vecPattern>(tag_double, zPtr, i, restI);

    const auto *const __restrict ownedStateIntPtr = reinterpret_cast<const int64_t *>(ownedStatePtr);
    const VectorLong ownedStateILong =
        loadIRegister<remainder, reversed, vecPattern>(tag_long, ownedStateIntPtr, i, restI);
    const MaskLong notDummyMaskILong = highway::Ne(ownedStateILong, highway::Zero(tag_long));
    notDummyMaskI = highway::RebindMask(tag_double, notDummyMaskILong);
    if constexpr (calculateGlobals) {
      const VectorLong ownedStateVal = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
      ownedMaskI = highway::RebindMask(tag_double, highway::Eq(ownedStateILong, ownedStateVal));
    }
  }

  /**
   * Fills physics parameter registers for i (halfSigma and sqrtEpsilon) matching the vector pattern layout.
   */
  template <bool remainder, bool reversed, VectorizationPattern vecPattern>
  static void fillIPhysicsRegisters(const size_t i, const double *const __restrict halfSigPtr,
                                    const double *const __restrict sqrtEpsPtr, VectorDouble &halfSigI,
                                    VectorDouble &sqrtEpsI, const size_t restI) {
    halfSigI = loadIRegister<remainder, reversed, vecPattern>(tag_double, halfSigPtr, i, restI);
    sqrtEpsI = loadIRegister<remainder, reversed, vecPattern>(tag_double, sqrtEpsPtr, i, restI);
  }

  /**
   * Loads a single attribute array for the inner j-loop according to the specified VectorizationPattern.
   */
  template <bool remainder, VectorizationPattern vecPattern, class Tag, typename ValueType>
  static inline auto loadJRegister(Tag tag, const ValueType *const __restrict dataPtr, const size_t j,
                                   const size_t restJ) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      if constexpr (remainder) {
        return highway::LoadN(tag, &dataPtr[j], restJ);
      } else {
        return highway::LoadU(tag, &dataPtr[j]);
      }
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      const size_t numLanes = remainder ? restJ : _vecLengthDouble / 2;
      const auto loadedHalf = highway::LoadN(tag, &dataPtr[j], numLanes);
      return highway::ConcatLowerLower(tag, loadedHalf, loadedHalf);
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      const auto lowerHalf = highway::Set(tag, dataPtr[j]);
      if constexpr (remainder) {
        return highway::ConcatLowerLower(tag, highway::Zero(tag), lowerHalf);
      } else {
        const auto upperHalf = highway::Set(tag, dataPtr[j + 1]);
        return highway::ConcatLowerLower(tag, upperHalf, lowerHalf);
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      return highway::Set(tag, dataPtr[j]);
    }
  }

  /**
   * Fills the registers for the loop over j. This includes positions and ownership state.
   */
  template <bool remainder, VectorizationPattern vecPattern>
  static void fillJRegisters(const size_t j, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
                             const double *const __restrict z2Ptr, const int64_t *const __restrict ownedStatePtr2,
                             VectorDouble &x2, VectorDouble &y2, VectorDouble &z2, MaskDouble &notDummyMaskJ,
                             [[maybe_unused]] MaskDouble &ownedMaskJ, const unsigned int restJ) {
    x2 = loadJRegister<remainder, vecPattern>(tag_double, x2Ptr, j, restJ);
    y2 = loadJRegister<remainder, vecPattern>(tag_double, y2Ptr, j, restJ);
    z2 = loadJRegister<remainder, vecPattern>(tag_double, z2Ptr, j, restJ);

    const VectorLong ownedStateJLong = loadJRegister<remainder, vecPattern>(tag_long, ownedStatePtr2, j, restJ);
    const MaskLong notDummyMaskJLong = highway::Ne(ownedStateJLong, highway::Zero(tag_long));
    notDummyMaskJ = highway::RebindMask(tag_double, notDummyMaskJLong);
    if constexpr (calculateGlobals) {
      const VectorLong ownedStateVal = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
      ownedMaskJ = highway::RebindMask(tag_double, highway::Eq(ownedStateJLong, ownedStateVal));
    }
  }

  /**
   * Fills physics parameter registers for j (halfSigma and sqrtEpsilon) matching the vector pattern layout.
   */
  template <bool remainder, VectorizationPattern vecPattern>
  static void fillJPhysicsRegisters(const size_t j, const double *const __restrict halfSigPtr,
                                    const double *const __restrict sqrtEpsPtr, VectorDouble &halfSigJ,
                                    VectorDouble &sqrtEpsJ, const unsigned int restJ) {
    halfSigJ = loadJRegister<remainder, vecPattern>(tag_double, halfSigPtr, j, restJ);
    sqrtEpsJ = loadJRegister<remainder, vecPattern>(tag_double, sqrtEpsPtr, j, restJ);
  }

  template <bool remainder, bool subtract, typename Tag>
  static inline void updateForceArrays(const Tag tag, const size_t offset, double *const __restrict fxPtr,
                                       double *const __restrict fyPtr, double *const __restrict fzPtr,
                                       const highway::VFromD<Tag> &dfx, const highway::VFromD<Tag> &dfy,
                                       const highway::VFromD<Tag> &dfz, const size_t lanes) {
    const auto oldFx = remainder ? highway::LoadN(tag, &fxPtr[offset], lanes) : highway::LoadU(tag, &fxPtr[offset]);
    const auto oldFy = remainder ? highway::LoadN(tag, &fyPtr[offset], lanes) : highway::LoadU(tag, &fyPtr[offset]);
    const auto oldFz = remainder ? highway::LoadN(tag, &fzPtr[offset], lanes) : highway::LoadU(tag, &fzPtr[offset]);

    const auto newFx = subtract ? highway::Sub(oldFx, dfx) : highway::Add(oldFx, dfx);
    const auto newFy = subtract ? highway::Sub(oldFy, dfy) : highway::Add(oldFy, dfy);
    const auto newFz = subtract ? highway::Sub(oldFz, dfz) : highway::Add(oldFz, dfz);

    if constexpr (remainder) {
      highway::StoreN(newFx, tag, &fxPtr[offset], lanes);
      highway::StoreN(newFy, tag, &fyPtr[offset], lanes);
      highway::StoreN(newFz, tag, &fzPtr[offset], lanes);
    } else {
      highway::StoreU(newFx, tag, &fxPtr[offset]);
      highway::StoreU(newFy, tag, &fyPtr[offset]);
      highway::StoreU(newFz, tag, &fzPtr[offset]);
    }
  }

  template <bool remainder, VectorizationPattern vecPattern>
  static void handleNewton3Reduction(const VectorDouble &fx, const VectorDouble &fy, const VectorDouble &fz,
                                     double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                                     double *const __restrict fz2Ptr, const size_t j, const size_t rest) {
    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      updateForceArrays<remainder, true>(tag_double, j, fx2Ptr, fy2Ptr, fz2Ptr, fx, fy, fz, rest);
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
      updateForceArrays<remainder, true>(tag_double_half, j, fx2Ptr, fy2Ptr, fz2Ptr, fxCombined, fyCombined, fzCombined,
                                         lanes);
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
  static void reduceAccumulatedForce(const size_t i, double *const __restrict fxPtr, double *const __restrict fyPtr,
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

      const auto fxAccCombined = highway::Add(lowerFxAcc, upperFxAcc);
      const auto fyAccCombined = highway::Add(lowerFyAcc, upperFyAcc);
      const auto fzAccCombined = highway::Add(lowerFzAcc, upperFzAcc);

      const size_t index = reversed ? (remainder ? 0 : i - _vecLengthDouble / 2 + 1) : i;
      const size_t lanes = remainder ? restI : _vecLengthDouble / 2;

      updateForceArrays<remainder, false>(tag_double_half, index, fxPtr, fyPtr, fzPtr, fxAccCombined, fyAccCombined,
                                          fzAccCombined, lanes);
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      const size_t index = reversed ? (remainder ? 0 : i - _vecLengthDouble + 1) : i;
      updateForceArrays<remainder, false>(tag_double, index, fxPtr, fyPtr, fzPtr, fxAcc, fyAcc, fzAcc, restI);
    }
  }

  inline void computeGlobals(const VectorDouble &virialSumX, const VectorDouble &virialSumY,
                             const VectorDouble &virialSumZ, const VectorDouble &uPotSum) {
    const int threadnum = autopas::autopas_get_thread_num();

    _aosThreadData[threadnum].virialSum[0] += highway::ReduceSum(tag_double, virialSumX);
    _aosThreadData[threadnum].virialSum[1] += highway::ReduceSum(tag_double, virialSumY);
    _aosThreadData[threadnum].virialSum[2] += highway::ReduceSum(tag_double, virialSumZ);
    _aosThreadData[threadnum].potentialEnergySum += highway::ReduceSum(tag_double, uPotSum);
  }

  struct LJParameters {
    VectorDouble sigmaSquare;
    VectorDouble epsilon24;
    VectorDouble shift6;
  };

  static inline VectorDouble computeShift6SIMD(const VectorDouble &sigmaSquareds, const VectorDouble &epsilon24s,
                                               const VectorDouble &invCutoffSquared) {
    const auto sig2Cut2 = highway::Mul(sigmaSquareds, invCutoffSquared);
    const auto sig4Cut4 = highway::Mul(sig2Cut2, sig2Cut2);
    const auto sig6Cut6 = highway::Mul(sig4Cut4, sig2Cut2);
    const auto sig12Cut12 = highway::Mul(sig6Cut6, sig6Cut6);
    return highway::Mul(epsilon24s, highway::Sub(sig6Cut6, sig12Cut12));
  }

  struct HomogeneousParameters {
    const VectorDouble sigmaSquareConst;
    const VectorDouble epsilon24Const;
    const VectorDouble shift6Const;

    HomogeneousParameters(const double sigmaSquaredConst, const double epsilon24Const, const double shift6Const)
        : sigmaSquareConst{highway::Set(tag_double, sigmaSquaredConst)},
          epsilon24Const{highway::Set(tag_double, epsilon24Const)},
          shift6Const{highway::Set(tag_double, shift6Const)} {}

    template <bool remainderI, bool reversed, VectorizationPattern vecPattern>
    inline void prepareI([[maybe_unused]] const size_t i, [[maybe_unused]] const size_t restI) const {}

    template <bool remainderJ, VectorizationPattern vecPattern>
    inline LJParameters getJ([[maybe_unused]] const size_t j, [[maybe_unused]] const size_t restJ) const {
      return {sigmaSquareConst, epsilon24Const, shift6Const};
    }
  };

  struct HeterogeneousParameters {
    const double *const __restrict halfSigma1Ptr;
    const double *const __restrict sqrtEpsilon1Ptr;
    const double *const __restrict halfSigma2Ptr;
    const double *const __restrict sqrtEpsilon2Ptr;
    VectorDouble invCutoffSquaredVec;

    VectorDouble halfSigmaI{highway::Zero(tag_double)};
    VectorDouble sqrtEpsilon24I{highway::Zero(tag_double)};

    HeterogeneousParameters(const double *const __restrict halfSig1, const double *const __restrict sqrtEps1,
                            const double *const __restrict halfSig2, const double *const __restrict sqrtEps2,
                            const double invCutoffSquared)
        : halfSigma1Ptr{halfSig1},
          sqrtEpsilon1Ptr{sqrtEps1},
          halfSigma2Ptr{halfSig2},
          sqrtEpsilon2Ptr{sqrtEps2},
          invCutoffSquaredVec{highway::Set(tag_double, invCutoffSquared)} {}

    template <bool remainderI, bool reversed, VectorizationPattern vecPattern>
    inline void prepareI(const size_t i, const size_t restI) {
      fillIPhysicsRegisters<remainderI, reversed, vecPattern>(i, halfSigma1Ptr, sqrtEpsilon1Ptr, halfSigmaI,
                                                              sqrtEpsilon24I, restI);
      sqrtEpsilon24I = highway::Mul(highway::Set(tag_double, 24.0), sqrtEpsilon24I);
    }

    template <bool remainderJ, VectorizationPattern vecPattern>
    inline LJParameters getJ(const size_t j, const size_t restJ) const {
      VectorDouble halfSigJ;
      VectorDouble sqrtEpsJ;
      fillJPhysicsRegisters<remainderJ, vecPattern>(j, halfSigma2Ptr, sqrtEpsilon2Ptr, halfSigJ, sqrtEpsJ, restJ);
      const auto sigma = highway::Add(halfSigmaI, halfSigJ);
      const auto sigSq = highway::Mul(sigma, sigma);
      const auto eps24 = highway::Mul(sqrtEpsilon24I, sqrtEpsJ);
      VectorDouble shift6 = highway::Zero(tag_double);
      if constexpr (calculateGlobals and applyShift) {
        shift6 = computeShift6SIMD(sigSq, eps24, invCutoffSquaredVec);
      }
      return {sigSq, eps24, shift6};
    }
  };

  struct HomogeneousVerletParameters {
    VectorDouble sigSqConst;
    VectorDouble eps24Const;
    VectorDouble shift6Const;

    HomogeneousVerletParameters(const double sigmaSquaredConst, const double epsilon24Const, const double shift6Const)
        : sigSqConst{highway::Set(tag_double, sigmaSquaredConst)},
          eps24Const{highway::Set(tag_double, epsilon24Const)},
          shift6Const{highway::Set(tag_double, shift6Const)} {}

    template <bool remainder>
    inline LJParameters getJ([[maybe_unused]] const VectorLong &indices,
                             [[maybe_unused]] const MaskDouble &restMaskDouble) const {
      return {sigSqConst, eps24Const, shift6Const};
    }
  };

  struct HeterogeneousVerletParameters {
    VectorDouble halfSigmaI;
    VectorDouble sqrtEpsilon24I;
    const int64_t *const __restrict typeIDPtr;
    const double *const __restrict halfSigmaPerTypePtr;
    const double *const __restrict sqrtEpsilonPerTypePtr;
    const VectorDouble invCutoffSquaredVec;

    HeterogeneousVerletParameters(const double halfSigmaIScalar, const double sqrtEpsilonIScalar,
                                  const int64_t *const __restrict typeIDPtr,
                                  const double *const __restrict halfSigmaPerType,
                                  const double *const __restrict sqrtEpsilonPerType, const double invCutoffSquared)
        : halfSigmaI{highway::Set(tag_double, halfSigmaIScalar)},
          sqrtEpsilon24I{highway::Set(tag_double, 24.0 * sqrtEpsilonIScalar)},
          typeIDPtr{typeIDPtr},
          halfSigmaPerTypePtr{halfSigmaPerType},
          sqrtEpsilonPerTypePtr{sqrtEpsilonPerType},
          invCutoffSquaredVec{highway::Set(tag_double, invCutoffSquared)} {}

    template <bool remainder>
    inline LJParameters getJ(const VectorLong &indices, const MaskDouble &restMaskDouble) const {
      VectorLong typeJ;
      if constexpr (remainder) {
        const MaskLong restMaskLong = highway::RebindMask(tag_long, restMaskDouble);
        typeJ = highway::MaskedGatherIndex(restMaskLong, tag_long, typeIDPtr, indices);
      } else {
        typeJ = highway::GatherIndex(tag_long, typeIDPtr, indices);
      }
      const VectorDouble halfSigmaJ = highway::GatherIndex(tag_double, halfSigmaPerTypePtr, typeJ);
      const VectorDouble sqrtEpsilonJ = highway::GatherIndex(tag_double, sqrtEpsilonPerTypePtr, typeJ);

      const auto sigma = highway::Add(halfSigmaI, halfSigmaJ);
      const auto sigSq = highway::Mul(sigma, sigma);
      const auto eps24 = highway::Mul(sqrtEpsilon24I, sqrtEpsilonJ);
      VectorDouble shift6 = highway::Zero(tag_double);
      if constexpr (calculateGlobals and applyShift) {
        shift6 = computeShift6SIMD(sigSq, eps24, invCutoffSquaredVec);
      }
      return {sigSq, eps24, shift6};
    }
  };

  /**
   * Outer loop body handler for SoAFunctorSingle and SoAFunctorPair.
   */
  template <bool reversed, bool newton3, bool remainderI, VectorizationPattern vecPattern, class LJParametersType>
  inline void handleILoopBody(
      const size_t i, const double *const __restrict xPtr1, const double *const __restrict yPtr1,
      const double *const __restrict zPtr1, const autopas::OwnershipState *const __restrict ownedStatePtr1,
      const double *const __restrict xPtr2, const double *const __restrict yPtr2, const double *const __restrict zPtr2,
      const autopas::OwnershipState *const __restrict ownedStatePtr2, double *const __restrict fxPtr1,
      double *const __restrict fyPtr1, double *const __restrict fzPtr1, double *const __restrict fxPtr2,
      double *const __restrict fyPtr2, double *const __restrict fzPtr2, LJParametersType &ljParameters,
      VectorDouble &virialSumX, VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &uPotSum,
      const size_t restI, const size_t jVecStart, const size_t jVecEnd) {
    VectorDouble fxAcc = highway::Zero(tag_double);
    VectorDouble fyAcc = highway::Zero(tag_double);
    VectorDouble fzAcc = highway::Zero(tag_double);

    MaskDouble notDummyMaskI;
    [[maybe_unused]] MaskDouble ownedMaskI;
    VectorDouble x1 = highway::Zero(tag_double);
    VectorDouble y1 = highway::Zero(tag_double);
    VectorDouble z1 = highway::Zero(tag_double);

    fillIRegisters<remainderI, reversed, vecPattern>(i, xPtr1, yPtr1, zPtr1, ownedStatePtr1, x1, y1, z1, notDummyMaskI,
                                                     ownedMaskI, restI);

    ljParameters.template prepareI<remainderI, reversed, vecPattern>(i, restI);

    const size_t jStep = jStepSize<vecPattern>();
    const size_t restJ = jVecEnd & (jStep - 1);
    const size_t jMainEnd = jVecEnd - restJ;
    size_t j = jVecStart;
    for (; j < jMainEnd; j += jStep) {
      SoAKernel<newton3, false, vecPattern>(j, notDummyMaskI, ownedMaskI,
                                            reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2,
                                            zPtr2, fxPtr2, fyPtr2, fzPtr2, ljParameters, fxAcc, fyAcc, fzAcc,
                                            virialSumX, virialSumY, virialSumZ, uPotSum, 0);
    }

    if (restJ > 0) {
      SoAKernel<newton3, true, vecPattern>(j, notDummyMaskI, ownedMaskI,
                                           reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2,
                                           zPtr2, fxPtr2, fyPtr2, fzPtr2, ljParameters, fxAcc, fyAcc, fzAcc, virialSumX,
                                           virialSumY, virialSumZ, uPotSum, restJ);
    }

    reduceAccumulatedForce<reversed, remainderI, vecPattern>(i, fxPtr1, fyPtr1, fzPtr1, fxAcc, fyAcc, fzAcc, restI);
  }

  /**
   * Templatized implementation of SoAFunctorPair.
   */
  template <bool newton3, bool sorted, VectorizationPattern vecPattern>
  inline void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                 autopas::utils::optRef<const autopas::SoASortingData> sortingData = std::nullopt) {
    const auto soaSize1 = soa1.size();
    const auto soaSize2 = soa2.size();

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
    const auto *const __restrict typeID1Ptr = soa1.template begin<Particle_T::AttributeNames::typeId>();
    const auto *const __restrict typeID2Ptr = soa2.template begin<Particle_T::AttributeNames::typeId>();

    // -------------------------------------------------------------------------
    // Scalar computations for low particle counts
    // -------------------------------------------------------------------------
    if (soaSize1 <= 7 and soaSize2 <= 7) {
      const int threadID = autopas::autopas_get_thread_num();
      for (size_t i = 0; i < soaSize1; ++i) {
        if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) continue;
        const double xi = x1Ptr[i], yi = y1Ptr[i], zi = z1Ptr[i];
        double halfSigI = 0., sqrtEpsI = 0.;
        if constexpr (useMixing) {
          const size_t ti = typeID1Ptr[i];
          halfSigI = _halfSigmaPerType[ti];
          sqrtEpsI = _sqrtEpsilonPerType[ti];
        }
        double fxi = 0., fyi = 0., fzi = 0.;
        for (size_t j = 0; j < soaSize2; ++j) {
          if (ownedStatePtr2[j] == autopas::OwnershipState::dummy) continue;
          const double dx = xi - x2Ptr[j];
          const double dy = yi - y2Ptr[j];
          const double dz = zi - z2Ptr[j];
          const double dr2 = dx * dx + dy * dy + dz * dz;
          if (dr2 <= _cutoffSquareAoS) {
            double sigSq = _sigmaSquareAoS;
            double eps24 = _epsilon24AoS;
            double shift6 = _shift6AoS;
            if constexpr (useMixing) {
              const size_t tj = typeID2Ptr[j];
              const double sig = halfSigI + _halfSigmaPerType[tj];
              sigSq = sig * sig;
              eps24 = 24.0 * sqrtEpsI * _sqrtEpsilonPerType[tj];
              if constexpr (applyShift) {
                shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(eps24, sigSq, _cutoffSquareAoS);
              }
            }
            const double invdr2 = 1.0 / dr2;
            const double lj2 = sigSq * invdr2;
            const double lj4 = lj2 * lj2;
            const double lj6 = lj4 * lj2;
            const double lj12 = lj6 * lj6;
            const double lj12m6 = lj12 - lj6;
            const double fac = eps24 * (lj12 + lj12m6) * invdr2;
            const double fx = dx * fac, fy = dy * fac, fz = dz * fac;
            fxi += fx;
            fyi += fy;
            fzi += fz;
            if constexpr (newton3) {
              fx2Ptr[j] -= fx;
              fy2Ptr[j] -= fy;
              fz2Ptr[j] -= fz;
            }
            if constexpr (calculateGlobals) {
              const double pot6 = eps24 * lj12m6 + shift6;
              if (ownedStatePtr1[i] == autopas::OwnershipState::owned) {
                _aosThreadData[threadID].potentialEnergySum += pot6;
                _aosThreadData[threadID].virialSum[0] += fx * dx;
                _aosThreadData[threadID].virialSum[1] += fy * dy;
                _aosThreadData[threadID].virialSum[2] += fz * dz;
              }
              if constexpr (newton3) {
                if (ownedStatePtr2[j] == autopas::OwnershipState::owned) {
                  _aosThreadData[threadID].potentialEnergySum += pot6;
                  _aosThreadData[threadID].virialSum[0] += fx * dx;
                  _aosThreadData[threadID].virialSum[1] += fy * dy;
                  _aosThreadData[threadID].virialSum[2] += fz * dz;
                }
              }
            }
          }
        }
        fx1Ptr[i] += fxi;
        fy1Ptr[i] += fyi;
        fz1Ptr[i] += fzi;
      }
      return;
    }

    const std::ptrdiff_t startI =
        sorted && sortingData.has_value() ? static_cast<std::ptrdiff_t>(sortingData->get().startI) : 0;
    const std::ptrdiff_t endI = sorted && sortingData.has_value() ? static_cast<std::ptrdiff_t>(sortingData->get().endI)
                                                                  : static_cast<std::ptrdiff_t>(soaSize1);

    VectorDouble virialSumX = highway::Zero(tag_double);
    VectorDouble virialSumY = highway::Zero(tag_double);
    VectorDouble virialSumZ = highway::Zero(tag_double);
    VectorDouble uPotSum = highway::Zero(tag_double);

    const size_t iStep = iStepSize<vecPattern>();
    const size_t jStep = jStepSize<vecPattern>();

    const bool isHomogeneous = soaIsHomogeneous(typeID1Ptr, soaSize1) and soaIsHomogeneous(typeID2Ptr, soaSize2);

    auto runPairLoop = [&](auto &ljParameters) {
      std::ptrdiff_t i = startI;
      for (; i + static_cast<std::ptrdiff_t>(iStep) <= endI; i += static_cast<std::ptrdiff_t>(iStep)) {
        size_t jVecEnd{};
        size_t jVecStart = 0;
        if constexpr (sorted) {
          const auto &sd = sortingData->get();
          jVecEnd = sd.maxIndex[i + iStep - 1];
          jVecStart = sd.minIndex[i];
          if (jVecStart >= jVecEnd) {
            continue;
          }
          jVecStart &= ~(jStep - 1);
          if constexpr (vecPattern == VectorizationPattern::p1xVec) {
            if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
              continue;
            }
          }
        } else {
          jVecEnd = soaSize2;
        }
        handleILoopBody<false, newton3, false, vecPattern>(
            i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx1Ptr, fy1Ptr, fz1Ptr, fx2Ptr,
            fy2Ptr, fz2Ptr, ljParameters, virialSumX, virialSumY, virialSumZ, uPotSum, 0, jVecStart, jVecEnd);
      }
      if constexpr (vecPattern != VectorizationPattern::p1xVec) {
        const std::ptrdiff_t restI = endI - i;
        if (restI > 0) {
          size_t jVecEnd = soaSize2;
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
            handleILoopBody<false, newton3, true, vecPattern>(
                i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx1Ptr, fy1Ptr, fz1Ptr,
                fx2Ptr, fy2Ptr, fz2Ptr, ljParameters, virialSumX, virialSumY, virialSumZ, uPotSum,
                static_cast<size_t>(restI), jVecStart, jVecEnd);
          }
        }
      }
    };

    if (isHomogeneous) {
      const size_t firstType1 = useMixing ? typeID1Ptr[0] : 0;
      const size_t firstType2 = useMixing ? typeID2Ptr[0] : 0;
      double sigSq = _sigmaSquareAoS;
      double eps24 = _epsilon24AoS;
      double shift6 = _shift6AoS;
      if constexpr (useMixing) {
        const double sig = _halfSigmaPerType[firstType1] + _halfSigmaPerType[firstType2];
        sigSq = sig * sig;
        eps24 = 24.0 * _sqrtEpsilonPerType[firstType1] * _sqrtEpsilonPerType[firstType2];
        if constexpr (applyShift) {
          shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(eps24, sigSq, _cutoffSquareAoS);
        }
      }
      HomogeneousParameters homogeneousParams{sigSq, eps24, shift6};
      runPairLoop(homogeneousParams);
    } else {
      MixingBuffers<512> mixingBuffers1;
      MixingBuffers<512> mixingBuffers2;
      double *halfSig1Ptr = nullptr, *sqrtEps1Ptr = nullptr;
      double *halfSig2Ptr = nullptr, *sqrtEps2Ptr = nullptr;
      mixingBuffers1.getMixingPointers(typeID1Ptr, soaSize1, _halfSigmaPerType, _sqrtEpsilonPerType, halfSig1Ptr,
                                       sqrtEps1Ptr);
      mixingBuffers2.getMixingPointers(typeID2Ptr, soaSize2, _halfSigmaPerType, _sqrtEpsilonPerType, halfSig2Ptr,
                                       sqrtEps2Ptr);

      HeterogeneousParameters heterogeneousParams{halfSig1Ptr, sqrtEps1Ptr, halfSig2Ptr, sqrtEps2Ptr,
                                                  _invCutoffSquareAoS};
      runPairLoop(heterogeneousParams);
    }

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, uPotSum);
    }
  }

  /**
   * Actual inner kernel of the SoAFunctors
   */
  template <bool newton3, bool remainderJ, VectorizationPattern vecPattern, typename LJParametersType>
  inline void SoAKernel(const size_t j, const MaskDouble &notDummyMaskI, [[maybe_unused]] const MaskDouble &ownedMaskI,
                        const int64_t *const __restrict ownedStatePtr2, const VectorDouble &x1, const VectorDouble &y1,
                        const VectorDouble &z1, const double *const __restrict x2Ptr,
                        const double *const __restrict y2Ptr, const double *const __restrict z2Ptr,
                        double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                        double *const __restrict fz2Ptr, LJParametersType &ljParams, VectorDouble &fxAcc,
                        VectorDouble &fyAcc, VectorDouble &fzAcc, VectorDouble &virialSumX, VectorDouble &virialSumY,
                        VectorDouble &virialSumZ, VectorDouble &uPotSum, const unsigned int restJ) {
    VectorDouble x2;
    VectorDouble y2;
    VectorDouble z2;
    MaskDouble notDummyMaskJ;
    [[maybe_unused]] MaskDouble ownedMaskJ;

    fillJRegisters<remainderJ, vecPattern>(j, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, x2, y2, z2, notDummyMaskJ,
                                           ownedMaskJ, restJ);

    // distance calculations
    const auto drX = highway::Sub(x1, x2);
    const auto drY = highway::Sub(y1, y2);
    const auto drZ = highway::Sub(z1, z2);

    const auto dr2 = highway::MulAdd(drZ, drZ, highway::MulAdd(drY, drY, highway::Mul(drX, drX)));

    const auto dummyMask = highway::And(notDummyMaskI, notDummyMaskJ);
    const auto cutoffDummyMask = highway::MaskedLe(dummyMask, dr2, _cutoffSquaredVec);

    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
      return;
    }

    const auto [sigmaSquareds, epsilon24s, shift6s] = ljParams.template getJ<remainderJ, vecPattern>(j, restJ);

    // compute LJ Potential
    const auto invDr2 = highway::Div(highway::Set(tag_double, 1.0), dr2);
    const auto lj2 = highway::Mul(sigmaSquareds, invDr2);
    const auto lj4 = highway::Mul(lj2, lj2);
    const auto lj6 = highway::Mul(lj2, lj4);
    const auto lj12 = highway::Mul(lj6, lj6);
    const auto lj12m6 = highway::Sub(lj12, lj6);
    const auto lj12m6alj12 = highway::Add(lj12m6, lj12);
    const auto lj12m6alj12e = highway::Mul(lj12m6alj12, epsilon24s);
    const auto fac = highway::Mul(lj12m6alj12e, invDr2);

    const auto facMasked = highway::IfThenElseZero(cutoffDummyMask, fac);

    const VectorDouble fx = highway::Mul(drX, facMasked);
    const VectorDouble fy = highway::Mul(drY, facMasked);
    const VectorDouble fz = highway::Mul(drZ, facMasked);

    fxAcc = highway::Add(fxAcc, fx);
    fyAcc = highway::Add(fyAcc, fy);
    fzAcc = highway::Add(fzAcc, fz);

    if constexpr (newton3) {
      handleNewton3Reduction<remainderJ, vecPattern>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, restJ);
    }

    if constexpr (calculateGlobals) {
      auto virialX = highway::Mul(fx, drX);
      auto virialY = highway::Mul(fy, drY);
      auto virialZ = highway::Mul(fz, drZ);

      auto uPot = highway::MulAdd(epsilon24s, lj12m6, shift6s);
      auto uPotMasked = highway::IfThenElseZero(cutoffDummyMask, uPot);

      auto energyFactor = highway::MaskedSet(tag_double, ownedMaskI, 1.0);

      if constexpr (newton3) {
        energyFactor = highway::Add(energyFactor, highway::MaskedSet(tag_double, ownedMaskJ, 1.0));
      }

      uPotSum = highway::MulAdd(energyFactor, uPotMasked, uPotSum);
      virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
      virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
      virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
    }
  }

  template <bool newton3, bool remainder = false, typename LJParametersType>
  inline void SoAKernelVerlet(const size_t i, const size_t j, const MaskDouble &notDummyMaskI,
                              [[maybe_unused]] const MaskDouble &ownedMaskI, const VectorDouble &x1,
                              const VectorDouble &y1, const VectorDouble &z1, const double *const __restrict xPtr,
                              const double *const __restrict yPtr, const double *const __restrict zPtr,
                              const int64_t *const __restrict ownedStatePtr, double *const __restrict fxPtr,
                              double *const __restrict fyPtr, double *const __restrict fzPtr,
                              const size_t *const __restrict neighborList, const LJParametersType &ljParams,
                              VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc, VectorDouble &virialSumX,
                              VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &uPotSum,
                              [[maybe_unused]] const size_t rest = 0) const {
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

    if constexpr (remainder) {
      x2 = highway::MaskedGatherIndex(restMaskDouble, tag_double, xPtr, indices);
      y2 = highway::MaskedGatherIndex(restMaskDouble, tag_double, yPtr, indices);
      z2 = highway::MaskedGatherIndex(restMaskDouble, tag_double, zPtr, indices);
      ownedState2 = highway::MaskedGatherIndex(restMaskLong, tag_long, ownedStatePtr, indices);
    } else {
      x2 = highway::GatherIndex(tag_double, xPtr, indices);
      y2 = highway::GatherIndex(tag_double, yPtr, indices);
      z2 = highway::GatherIndex(tag_double, zPtr, indices);
      ownedState2 = highway::GatherIndex(tag_long, ownedStatePtr, indices);
    }

    const MaskLong notDummyMaskJLong = highway::Ne(ownedState2, highway::Zero(tag_long));
    const MaskDouble notDummyMaskJ = highway::RebindMask(tag_double, notDummyMaskJLong);

    const auto drX = highway::Sub(x1, x2);
    const auto drY = highway::Sub(y1, y2);
    const auto drZ = highway::Sub(z1, z2);

    const auto dr2 = highway::MulAdd(drZ, drZ, highway::MulAdd(drY, drY, highway::Mul(drX, drX)));

    const auto dummyMask = highway::And(notDummyMaskI, notDummyMaskJ);
    const auto cutoffDummyMask = highway::MaskedLe(dummyMask, dr2, _cutoffSquaredVec);

    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
      return;
    }
    const auto [sigmaSquareds, epsilon24s, shift6s] = ljParams.template getJ<remainder>(indices, restMaskDouble);

    const auto invDr2 = highway::Div(highway::Set(tag_double, 1.0), dr2);
    const auto lj2 = highway::Mul(sigmaSquareds, invDr2);
    const auto lj4 = highway::Mul(lj2, lj2);
    const auto lj6 = highway::Mul(lj2, lj4);
    const auto lj12 = highway::Mul(lj6, lj6);
    const auto lj12m6 = highway::Sub(lj12, lj6);
    const auto lj12m6alj12 = highway::Add(lj12m6, lj12);
    const auto lj12m6alj12e = highway::Mul(lj12m6alj12, epsilon24s);
    const auto fac = highway::Mul(lj12m6alj12e, invDr2);

    const auto facMasked = highway::IfThenElseZero(cutoffDummyMask, fac);

    const VectorDouble fx = highway::Mul(drX, facMasked);
    const VectorDouble fy = highway::Mul(drY, facMasked);
    const VectorDouble fz = highway::Mul(drZ, facMasked);

    fxAcc = highway::Add(fxAcc, fx);
    fyAcc = highway::Add(fyAcc, fy);
    fzAcc = highway::Add(fzAcc, fz);

    if constexpr (newton3) {
      if constexpr (remainder) {
        const auto fx2Old = highway::MaskedGatherIndex(restMaskDouble, tag_double, fxPtr, indices);
        const auto fy2Old = highway::MaskedGatherIndex(restMaskDouble, tag_double, fyPtr, indices);
        const auto fz2Old = highway::MaskedGatherIndex(restMaskDouble, tag_double, fzPtr, indices);

        const auto fx2New = highway::Sub(fx2Old, fx);
        const auto fy2New = highway::Sub(fy2Old, fy);
        const auto fz2New = highway::Sub(fz2Old, fz);

        highway::MaskedScatterIndex(fx2New, restMaskDouble, tag_double, fxPtr, indices);
        highway::MaskedScatterIndex(fy2New, restMaskDouble, tag_double, fyPtr, indices);
        highway::MaskedScatterIndex(fz2New, restMaskDouble, tag_double, fzPtr, indices);
      } else {
        const auto fx2Old = highway::GatherIndex(tag_double, fxPtr, indices);
        const auto fy2Old = highway::GatherIndex(tag_double, fyPtr, indices);
        const auto fz2Old = highway::GatherIndex(tag_double, fzPtr, indices);

        const auto fx2New = highway::Sub(fx2Old, fx);
        const auto fy2New = highway::Sub(fy2Old, fy);
        const auto fz2New = highway::Sub(fz2Old, fz);

        highway::ScatterIndex(fx2New, tag_double, fxPtr, indices);
        highway::ScatterIndex(fy2New, tag_double, fyPtr, indices);
        highway::ScatterIndex(fz2New, tag_double, fzPtr, indices);
      }
    }

    if constexpr (calculateGlobals) {
      auto virialX = highway::Mul(fx, drX);
      auto virialY = highway::Mul(fy, drY);
      auto virialZ = highway::Mul(fz, drZ);

      auto uPot = highway::MulAdd(epsilon24s, lj12m6, shift6s);
      auto uPotMasked = highway::IfThenElseZero(cutoffDummyMask, uPot);

      auto energyFactor = highway::MaskedSet(tag_double, ownedMaskI, 1.0);

      if constexpr (newton3) {
        const VectorLong ownedStateVal = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
        const MaskDouble ownedMaskJ = highway::RebindMask(tag_double, highway::Eq(ownedState2, ownedStateVal));
        energyFactor = highway::Add(energyFactor, highway::MaskedSet(tag_double, ownedMaskJ, 1.0));
      }

      uPotSum = highway::MulAdd(energyFactor, uPotMasked, uPotSum);
      virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
      virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
      virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
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
  inline void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                                   std::span<const size_t> neighborList) {
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
    const auto ownedI = static_cast<int64_t>(ownedStatePtr[indexFirst]);
    const VectorDouble ownedStateI = highway::Set(tag_double, static_cast<double>(ownedI));
    const MaskDouble notDummyMaskI = highway::Ne(ownedStateI, highway::Zero(tag_double));
    [[maybe_unused]] MaskDouble ownedMaskI;
    if constexpr (calculateGlobals) {
      const auto ownedVal = static_cast<double>(static_cast<int64_t>(autopas::OwnershipState::owned));
      ownedMaskI = highway::Eq(ownedStateI, highway::Set(tag_double, ownedVal));
    }

    const bool isHomogeneous = not useMixing or _halfSigmaPerType.size() <= 1;

    auto runVerletLoop = [&](const auto &ljParameters) {
      size_t j = 0;
      const size_t neighborListSize = neighborList.size();
      const size_t vecEnd = (neighborListSize / _vecLengthDouble) * _vecLengthDouble;

      for (; j < vecEnd; j += _vecLengthDouble) {
        SoAKernelVerlet<newton3, false>(indexFirst, j, notDummyMaskI, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr,
                                        reinterpret_cast<const int64_t *>(ownedStatePtr), fxPtr, fyPtr, fzPtr,
                                        neighborList.data(), ljParameters, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                        virialSumZ, uPotSum);
      }

      const size_t rest = neighborListSize & (_vecLengthDouble - 1);
      if (rest > 0) {
        SoAKernelVerlet<newton3, true>(indexFirst, j, notDummyMaskI, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr,
                                       reinterpret_cast<const int64_t *>(ownedStatePtr), fxPtr, fyPtr, fzPtr,
                                       neighborList.data(), ljParameters, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                       virialSumZ, uPotSum, rest);
      }
    };

    if (isHomogeneous) {
      const size_t firstType = useMixing ? typeIDPtr[indexFirst] : 0;
      double sigSq = _sigmaSquareAoS;
      double eps24 = _epsilon24AoS;
      double shift6 = _shift6AoS;
      if constexpr (useMixing) {
        const double sig = 2.0 * _halfSigmaPerType[firstType];
        sigSq = sig * sig;
        eps24 = 24.0 * _sqrtEpsilonPerType[firstType] * _sqrtEpsilonPerType[firstType];
        if constexpr (applyShift) {
          shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(eps24, sigSq, _cutoffSquareAoS);
        }
      }
      HomogeneousVerletParameters homogeneousVerletParams{sigSq, eps24, shift6};
      runVerletLoop(homogeneousVerletParams);
    } else {
      const size_t typeI = typeIDPtr[indexFirst];
      HeterogeneousVerletParameters heterogeneousVerletParams{
          _halfSigmaPerType[typeI], _sqrtEpsilonPerType[typeI], reinterpret_cast<const int64_t *>(typeIDPtr),
          _halfSigmaPerType.data(), _sqrtEpsilonPerType.data(), _invCutoffSquareAoS};
      runVerletLoop(heterogeneousVerletParams);
    }

    fxPtr[indexFirst] += highway::ReduceSum(tag_double, fxAcc);
    fyPtr[indexFirst] += highway::ReduceSum(tag_double, fyAcc);
    fzPtr[indexFirst] += highway::ReduceSum(tag_double, fzAcc);

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, uPotSum);
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
   * @copydoc autopas::Functor::supportsSoASorting
   */
  static constexpr bool supportsSoASorting = true;

  /**
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
      _potentialEnergySum *= 0.5;
      _virialSum *= 0.5;

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
  [[nodiscard]] double getPotentialEnergy() const {
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
  [[nodiscard]] double getVirial() const {
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
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param epsilon24
   * @param sigmaSquare
   */
  void setParticleProperties(const double epsilon24, const double sigmaSquare) {
    _epsilon24AoS = epsilon24;
    _sigmaSquareAoS = sigmaSquare;
    _halfSigmaPerType = {std::sqrt(sigmaSquare) * 0.5};
    _sqrtEpsilonPerType = {std::sqrt(epsilon24 / 24.0)};
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

  const VectorDouble _cutoffSquaredVec{highway::Set(tag_double, 0.)};

  // epsilon, sigma and shift6 used in the AoS functor.
  double _sigmaSquareAoS{0.};
  double _epsilon24AoS{0.};
  double _shift6AoS{0.};
  std::vector<double> _halfSigmaPerType{};
  std::vector<double> _sqrtEpsilonPerType{};

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
