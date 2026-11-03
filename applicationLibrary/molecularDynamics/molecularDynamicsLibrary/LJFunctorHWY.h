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

 public:
  /**
   * Helper struct to manage temporary buffers for projected halfSigma and sqrtEpsilon values.
   *
   * @tparam StackCapacity Capacity of the stack-allocated buffer before falling back to thread-local heap memory.
   */
  template <size_t StackCapacity = 256>
  struct MixingBuffers {
    /**
     * Stack storage for projected halfSigma values.
     */
    alignas(64) std::array<double, StackCapacity> halfSigmaStack;
    /**
     * Stack storage for projected sqrtEpsilon values.
     */
    alignas(64) std::array<double, StackCapacity> sqrtEpsilonStack;

    /**
     * Fills the provided pointers with the projected halfSigma and sqrtEpsilon values for the given typeIDs.
     * If the number of particles is less than or equal to StackCapacity, the values are stored in the stack buffer.
     * Otherwise, they are stored in a thread-local heap buffer.
     *
     * @tparam BufferIndex An optional template parameter to differentiate between multiple MixingBuffers instances.
     * @param typeIDs Pointer to the array of type IDs for the particles.
     * @param soaSize The number of particles in the SoA.
     * @param halfSigmaPerType Vector containing halfSigma values for each particle type.
     * @param sqrtEpsilonPerType Vector containing sqrtEpsilon values for each particle type.
     * @param halfSigmaPtr Reference to a pointer that will be set to point to the projected halfSigma values.
     * @param sqrtEpsilonPtr Reference to a pointer that will be set to point to the projected sqrtEpsilon values.
     */
    template <size_t BufferIndex = 0>
    void getMixingPointers(const size_t *__restrict typeIDs, const size_t soaSize,
                           const std::vector<double> &halfSigmaPerType, const std::vector<double> &sqrtEpsilonPerType,
                           double *&halfSigmaPtr, double *&sqrtEpsilonPtr) {
      if (soaSize <= StackCapacity) [[likely]] {
        halfSigmaPtr = halfSigmaStack.data();
        sqrtEpsilonPtr = sqrtEpsilonStack.data();
      } else [[unlikely]] {
        thread_local std::vector<double, autopas::AlignedAllocator<double, 64>> halfSigmaHeap;
        thread_local std::vector<double, autopas::AlignedAllocator<double, 64>> sqrtEpsilonHeap;
        if (halfSigmaHeap.size() < soaSize) {
          halfSigmaHeap.resize(soaSize);
          sqrtEpsilonHeap.resize(soaSize);
        }
        halfSigmaPtr = halfSigmaHeap.data();
        sqrtEpsilonPtr = sqrtEpsilonHeap.data();
      }
#pragma omp simd
      for (size_t k = 0; k < soaSize; ++k) {
        const auto typeID = typeIDs[k];
        halfSigmaPtr[k] = halfSigmaPerType[typeID];
        sqrtEpsilonPtr[k] = sqrtEpsilonPerType[typeID];
      }
    }
  };

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
  explicit LJFunctorHWY(double cutoff,
                        const std::optional<std::reference_wrapper<ParticlePropertiesLibrary<double, size_t>>>
                            particlePropertiesLibrary = std::nullopt)
      : autopas::PairwiseFunctor<Particle_T, LJFunctorHWY>(cutoff),
        _cutoffSquared{cutoff * cutoff},
        _invCutoffSquared{1. / (cutoff * cutoff)} {
    if (calculateGlobals) {
      _threadData.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      AutoPasLog(DEBUG, "Using LJFunctorHWY with countFLOPs but FLOP counting is not implemented.");
    }

    if constexpr (useMixing) {
      if (not particlePropertiesLibrary.has_value()) {
        throw std::runtime_error("Mixing is enabled but no ParticlePropertiesLibrary was provided!");
      }
      initFromPPL(particlePropertiesLibrary->get());
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

 private:
  /**
   * Checks whether all particles within an SoA have identical type IDs.
   *
   * @param typeIDPtr Pointer to the particle type IDs array.
   * @param soaSize Number of particles in the SoA.
   * @return True if the cell is homogeneous (or if mixing is disabled), false otherwise.
   */
  static inline bool soaIsHomogeneous(const size_t *const __restrict typeIDPtr, const size_t soaSize) {
    if constexpr (useMixing) {
      const auto firstType = typeIDPtr[0];
      for (size_t k = 1; k < soaSize; ++k) {
        if (typeIDPtr[k] != firstType) {
          return false;
        }
      }
    }
    return true;
  }

 public:
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

    if (dr2 > _cutoffSquared) {
      return;
    }

    auto sigmaSquared = _sigmaSquared;
    auto epsilon24 = _epsilon24;
    auto shift6 = _shift6;
    if constexpr (useMixing) {
      const auto typeI = i.getTypeId();
      const auto typeJ = j.getTypeId();
      const auto sigma = _halfSigmaPerType[typeI] + _halfSigmaPerType[typeJ];
      sigmaSquared = sigma * sigma;
      epsilon24 = 24.0 * _sqrtEpsilonPerType[typeI] * _sqrtEpsilonPerType[typeJ];
      if constexpr (calculateGlobals and applyShift) {
        shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquared);
      }
    }

    const double invdr2 = 1. / dr2;
    double lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    const double lj12 = lj6 * lj6;
    const double lj12m6 = lj12 - lj6;
    const double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    const auto f = dr * fac;
    i.addF(f);
    if (newton3) {
      j.subF(f);
    }
    if (calculateGlobals) {
      const auto virial = dr * f;
      const double potentialEnergy6 = epsilon24 * lj12m6 + shift6;

      const int threadID = autopas::autopas_get_thread_num();

      if (i.isOwned()) {
        _threadData[threadID].potentialEnergySum += potentialEnergy6;
        _threadData[threadID].virialSum += virial;
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        // for non-newton3 the division is in the post-processing step.
        _threadData[threadID].potentialEnergySum += potentialEnergy6;
        _threadData[threadID].virialSum += virial;
      }
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always do a newton3 like traversal of the SoA using only the 1xVec Pattern.
   */
  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, const bool newton3) final {
    if (soa.size() == 0) return;

    // obtain iterators for the various values
    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    MixingBuffers mixingBuffers;
    double *halfSigmaPtr = nullptr;
    double *sqrtEpsilonPtr = nullptr;
    double homoSigmaSquared = _sigmaSquared;
    double homoEpsilon24 = _epsilon24;
    double homoShift6 = _shift6;

    if constexpr (useMixing) {
      const auto *const __restrict typeIDPtr = soa.template begin<Particle_T::AttributeNames::typeId>();
      if (soaIsHomogeneous(typeIDPtr, soa.size())) {
        const size_t type = soa.size() > 0 ? typeIDPtr[0] : 0;
        const double sigma = 2.0 * _halfSigmaPerType[type];
        homoSigmaSquared = sigma * sigma;
        homoEpsilon24 = 24.0 * _sqrtEpsilonPerType[type] * _sqrtEpsilonPerType[type];
        if constexpr (calculateGlobals and applyShift) {
          homoShift6 =
              ParticlePropertiesLibrary<double, size_t>::calcShift6(homoEpsilon24, homoSigmaSquared, _cutoffSquared);
        }
      } else {
        mixingBuffers.getMixingPointers(typeIDPtr, soa.size(), _halfSigmaPerType, _sqrtEpsilonPerType, halfSigmaPtr,
                                        sqrtEpsilonPtr);
      }
    }

    // initialize and declare vector variables
    auto virialSumX = highway::Zero(tag_double);
    auto virialSumY = highway::Zero(tag_double);
    auto virialSumZ = highway::Zero(tag_double);
    auto potentialEnergySum = highway::Zero(tag_double);

    for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(soa.size()) - 1; i >= 0; i -= 1) {
      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");

      handleILoopBody<true, true, false, VectorizationPattern::p1xVec>(
          i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, fxPtr, fyPtr, fzPtr,
          halfSigmaPtr, sqrtEpsilonPtr, halfSigmaPtr, sqrtEpsilonPtr, virialSumX, virialSumY, virialSumZ,
          potentialEnergySum, 0, 0, i, homoSigmaSquared, homoEpsilon24, homoShift6);
    }

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, potentialEnergySum);
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
   */
  inline void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                             bool newton3) final {
    switch (_vecPattern) {
      case VectorizationPattern::p1xVec: {
        if (newton3) {
          SoAFunctorPairImpl<true, false, VectorizationPattern::p1xVec>(soa1, soa2);
        } else {
          SoAFunctorPairImpl<false, false, VectorizationPattern::p1xVec>(soa1, soa2);
        }
        break;
      }
      case VectorizationPattern::p2xVecDiv2: {
        if (newton3) {
          SoAFunctorPairImpl<true, false, VectorizationPattern::p2xVecDiv2>(soa1, soa2);
        } else {
          SoAFunctorPairImpl<false, false, VectorizationPattern::p2xVecDiv2>(soa1, soa2);
        }
        break;
      }
      case VectorizationPattern::pVecDiv2x2: {
        if (newton3) {
          SoAFunctorPairImpl<true, false, VectorizationPattern::pVecDiv2x2>(soa1, soa2);
        } else {
          SoAFunctorPairImpl<false, false, VectorizationPattern::pVecDiv2x2>(soa1, soa2);
        }
        break;
      }
      case VectorizationPattern::pVecx1: {
        if (newton3) {
          SoAFunctorPairImpl<true, false, VectorizationPattern::pVecx1>(soa1, soa2);
        } else {
          SoAFunctorPairImpl<false, false, VectorizationPattern::pVecx1>(soa1, soa2);
        }
        break;
      }
      default:
        autopas::utils::ExceptionHandler::exception("Unknown VectorizationPattern!");
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPairSorted()
   */
  inline void SoAFunctorPairSorted(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                                   const autopas::SoASortingData &sortingData, bool newton3) final {
    if (soa1.size() == 0 or soa2.size() == 0) {
      return;
    }
    switch (_vecPattern) {
      case VectorizationPattern::p1xVec: {
        if (newton3) {
          SoAFunctorPairImpl<true, true, VectorizationPattern::p1xVec>(soa1, soa2, sortingData);
        } else {
          SoAFunctorPairImpl<false, true, VectorizationPattern::p1xVec>(soa1, soa2, sortingData);
        }
        break;
      }
      case VectorizationPattern::p2xVecDiv2: {
        if (newton3) {
          SoAFunctorPairImpl<true, true, VectorizationPattern::p2xVecDiv2>(soa1, soa2, sortingData);
        } else {
          SoAFunctorPairImpl<false, true, VectorizationPattern::p2xVecDiv2>(soa1, soa2, sortingData);
        }
        break;
      }
      case VectorizationPattern::pVecDiv2x2: {
        if (newton3) {
          SoAFunctorPairImpl<true, true, VectorizationPattern::pVecDiv2x2>(soa1, soa2, sortingData);
        } else {
          SoAFunctorPairImpl<false, true, VectorizationPattern::pVecDiv2x2>(soa1, soa2, sortingData);
        }
        break;
      }
      case VectorizationPattern::pVecx1: {
        if (newton3) {
          SoAFunctorPairImpl<true, true, VectorizationPattern::pVecx1>(soa1, soa2, sortingData);
        } else {
          SoAFunctorPairImpl<false, true, VectorizationPattern::pVecx1>(soa1, soa2, sortingData);
        }
        break;
      }
      default:
        autopas::utils::ExceptionHandler::exception("Unknown VectorizationPattern!");
    }
  }

 private:
  /**
   * Initializes the functor's internal constants from the provided ParticlePropertiesLibrary.
   * @param ppl The ParticlePropertiesLibrary containing particle properties.
   */
  void initFromPPL(const ParticlePropertiesLibrary<double, size_t> &ppl) {
    const size_t nTypes = ppl.getNumberRegisteredSiteTypes();
    _halfSigmaPerType.resize(nTypes);
    _sqrtEpsilonPerType.resize(nTypes);
    for (size_t t = 0; t < nTypes; ++t) {
      _halfSigmaPerType[t] = ppl.getSigma(t) * 0.5;
      _sqrtEpsilonPerType[t] = std::sqrt(ppl.getEpsilon(t));
    }
    if (nTypes > 0) {
      _sigmaSquared = 4.0 * _halfSigmaPerType[0] * _halfSigmaPerType[0];
      _epsilon24 = 24.0 * _sqrtEpsilonPerType[0] * _sqrtEpsilonPerType[0];
      if constexpr (applyShift) {
        _shift6 = ParticlePropertiesLibrary<>::calcShift6(_epsilon24, _sigmaSquared, _cutoffSquared);
      }
    }
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
   * Checks whether the inner loop over j should be continued. This depends on the respective VectorizationPattern.
   *
   * @tparam vecPattern
   * @param i
   * @param j
   * @return  true if the inner loop over j should be continued, else false.
   */
  template <VectorizationPattern vecPattern>
  static constexpr bool checkSecondLoopCondition(std::ptrdiff_t i, size_t j) {
    // Round i down to the nearest multiple of jStep (the j-lane width) to get the exclusive upper bound.
    const auto jStep = static_cast<std::ptrdiff_t>(jStepSize<vecPattern>());
    const std::ptrdiff_t limit = i - (i % jStep);
    return j < static_cast<size_t>(limit);
  }
  /**
   * Loads a single attribute array for the outer i-loop according to the specified VectorizationPattern.
   *
   * @tparam remainder Whether this is a remainder loop iteration with fewer lanes than the vector width.
   * @tparam reversed Whether iterating backwards over the outer loop.
   * @tparam vecPattern The vectorization pattern.
   * @tparam Tag The Highway tag type defining vector lane width and element type.
   * @tparam ValueType The element type being loaded.
   * @param tag Tag descriptor for vector operations.
   * @param dataPtr Pointer to the contiguous array of data to load.
   * @param i Particle index in the outer i loop.
   * @param restI Number of valid elements in case remainder is true.
   * @return SIMD vector filled according to the vectorization pattern.
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
   * Depending on the vectorization pattern, this function fills the registers for the loop over i. This includes the
   * positions and the ownership state.
   *
   * @tparam remainder Whether this is a remainder loop iteration with fewer lanes than the vector width.
   * @tparam reversed Whether iterating backwards over the outer loop.
   * @tparam vecPattern The vectorization pattern.
   * @param i Current loop index for i.
   * @param xPtr pointer to the x coordinate of the particle at position i.
   * @param yPtr pointer to the y coordinate of the particle at position i.
   * @param zPtr pointer to the z coordinate of the particle at position i.
   * @param ownedStatePtr pointer to the ownership state of the particle at position i.
   * @param x1 The register to be filled for x coordinates.
   * @param y1 The register to be filled for y coordinates.
   * @param z1 The register to be filled for z coordinates.
   * @param notDummyMaskI The mask to be filled for the 'not dummy' ownership state.
   * @param ownedMaskI The mask to be filled for the 'owned' ownership state.
   * @param restI The number of lanes filled in case of a remainder loop.
   */
  template <bool remainder, bool reversed, VectorizationPattern vecPattern>
  static void fillIRegisters(const size_t i, const double *const __restrict xPtr, const double *const __restrict yPtr,
                             const double *const __restrict zPtr,
                             const autopas::OwnershipState *const __restrict ownedStatePtr, VectorDouble &x1,
                             VectorDouble &y1, VectorDouble &z1, MaskDouble &notDummyMaskI,
                             [[maybe_unused]] MaskDouble &ownedMaskI, const size_t restI) {
    const auto ownedStateILong = loadIRegister<remainder, reversed, vecPattern>(
        tag_long, reinterpret_cast<const int64_t *>(ownedStatePtr), i, restI);
    x1 = loadIRegister<remainder, reversed, vecPattern>(tag_double, xPtr, i, restI);
    y1 = loadIRegister<remainder, reversed, vecPattern>(tag_double, yPtr, i, restI);
    z1 = loadIRegister<remainder, reversed, vecPattern>(tag_double, zPtr, i, restI);

    const MaskLong notDummyMaskILong = highway::Ne(ownedStateILong, highway::Zero(tag_long));
    notDummyMaskI = highway::RebindMask(tag_double, notDummyMaskILong);
    if constexpr (calculateGlobals) {
      const VectorLong ownedStateVal = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
      ownedMaskI = highway::RebindMask(tag_double, highway::Eq(ownedStateILong, ownedStateVal));
    }
  }

  /**
   * Loads the physical parameters (halfSigma and sqrtEpsilon) for particle(s) i according to the vectorization pattern.
   *
   * @tparam remainder Whether this is a remainder loop iteration with fewer lanes than the vector width.
   * @tparam reversed Whether iterating backwards over the outer loop.
   * @tparam vecPattern The vectorization pattern.
   * @param i Particle index in the outer i loop.
   * @param halfSigPtr Pointer to the halfSigma array.
   * @param sqrtEpsPtr Pointer to the sqrtEpsilon array.
   * @param halfSigI Output SIMD vector for halfSigma of particle(s) i.
   * @param sqrtEpsI Output SIMD vector for sqrtEpsilon of particle(s) i.
   * @param restI Number of valid elements in case remainder is true.
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
   *
   * @tparam remainder Whether this is a remainder loop iteration with fewer lanes than the vector width.
   * @tparam vecPattern The vectorization pattern.
   * @tparam Tag The Highway tag type defining vector lane width and element type.
   * @tparam ValueType The element type being loaded.
   * @param tag Tag descriptor for vector operations.
   * @param dataPtr Pointer to the contiguous array of data to load.
   * @param j Particle index in the inner j loop.
   * @param restJ Number of valid elements in case remainder is true.
   * @return SIMD vector filled according to the vectorization pattern.
   */
  template <bool remainder, VectorizationPattern vecPattern, typename Tag, typename ValueType>
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
   *
   * @tparam remainder Whether this is a remainder loop iteration with fewer lanes than the vector width.
   * @tparam vecPattern The vectorization pattern.
   * @param j Current loop index for j.
   * @param x2Ptr Pointer to the x coordinate of the particle at position j.
   * @param y2Ptr Pointer to the y coordinate of the particle at position j.
   * @param z2Ptr Pointer to the z coordinate of the particle at position j.
   * @param ownedStatePtr2 Pointer to the ownership state of the particle at position j.
   * @param x2 The register to be filled for x coordinates.
   * @param y2 The register to be filled for y coordinates.
   * @param z2 The register to be filled for z coordinates.
   * @param notDummyMaskJ The mask to be filled for the 'not dummy' ownership state.
   * @param ownedMaskJ The mask to be filled for the 'owned' ownership state.
   * @param rest The number of lanes filled in case of a remainder loop.
   */
  template <bool remainder, VectorizationPattern vecPattern>
  static void fillJRegisters(const size_t j, const double *const __restrict x2Ptr, const double *const __restrict y2Ptr,
                             const double *const __restrict z2Ptr, const int64_t *const __restrict ownedStatePtr2,
                             VectorDouble &x2, VectorDouble &y2, VectorDouble &z2, MaskDouble &notDummyMaskJ,
                             [[maybe_unused]] MaskDouble &ownedMaskJ, const unsigned int rest) {
    const auto ownedStateJLong = loadJRegister<remainder, vecPattern>(tag_long, ownedStatePtr2, j, rest);
    x2 = loadJRegister<remainder, vecPattern>(tag_double, x2Ptr, j, rest);
    y2 = loadJRegister<remainder, vecPattern>(tag_double, y2Ptr, j, rest);
    z2 = loadJRegister<remainder, vecPattern>(tag_double, z2Ptr, j, rest);

    const MaskLong notDummyMaskJLong = highway::Ne(ownedStateJLong, highway::Zero(tag_long));
    notDummyMaskJ = highway::RebindMask(tag_double, notDummyMaskJLong);
    if constexpr (calculateGlobals) {
      const VectorLong ownedStateVal = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
      ownedMaskJ = highway::RebindMask(tag_double, highway::Eq(ownedStateJLong, ownedStateVal));
    }
  }

  /**
   * Loads the physical parameters (halfSigma and sqrtEpsilon) for particle(s) j according to the vectorization pattern.
   *
   * @tparam remainder Whether this is a remainder loop iteration with fewer lanes than the vector width.
   * @tparam vecPattern The vectorization pattern.
   * @param j Particle index in the inner j loop.
   * @param halfSigPtr Pointer to the halfSigma array.
   * @param sqrtEpsPtr Pointer to the sqrtEpsilon array.
   * @param halfSigJ Output SIMD vector for halfSigma of particle(s) j.
   * @param sqrtEpsJ Output SIMD vector for sqrtEpsilon of particle(s) j.
   * @param restJ Number of valid elements in case remainder is true.
   */
  template <bool remainder, VectorizationPattern vecPattern>
  static void fillJPhysicsRegisters(const size_t j, const double *const __restrict halfSigPtr,
                                    const double *const __restrict sqrtEpsPtr, VectorDouble &halfSigJ,
                                    VectorDouble &sqrtEpsJ, const unsigned int restJ) {
    halfSigJ = loadJRegister<remainder, vecPattern>(tag_double, halfSigPtr, j, restJ);
    sqrtEpsJ = loadJRegister<remainder, vecPattern>(tag_double, sqrtEpsPtr, j, restJ);
  }

  /**
   * Helper function to update force arrays (either subtracting or adding delta forces).
   * Loads the existing forces, performs addition/subtraction, and writes them back.
   *
   * @tparam remainder Whether this is a remainder loop iteration with fewer lanes than the vector width.
   * @tparam subtract Whether to subtract (true, for Newton3 on particle j) or add (false, for particle i).
   * @tparam Tag The Highway tag type defining vector lane width and element type.
   * @param tag Tag descriptor for vector operations.
   * @param offset Starting index offset in the force arrays.
   * @param fxPtr Pointer to the x-component force array.
   * @param fyPtr Pointer to the y-component force array.
   * @param fzPtr Pointer to the z-component force array.
   * @param dfx Vector of force deltas in x direction.
   * @param dfy Vector of force deltas in y direction.
   * @param dfz Vector of force deltas in z direction.
   * @param lanes Number of valid lanes in case remainder is true.
   */
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
                                     const VectorDouble &fyAcc, const VectorDouble &fzAcc, const size_t restI) {
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
                             const VectorDouble &virialSumZ, const VectorDouble &potentialEnergySum) {
    const int threadID = autopas::autopas_get_thread_num();

    _threadData[threadID].virialSum[0] += highway::ReduceSum(tag_double, virialSumX);
    _threadData[threadID].virialSum[1] += highway::ReduceSum(tag_double, virialSumY);
    _threadData[threadID].virialSum[2] += highway::ReduceSum(tag_double, virialSumZ);
    _threadData[threadID].potentialEnergySum += highway::ReduceSum(tag_double, potentialEnergySum);
  }

  /**
   * Struct holding Lennard-Jones parameters loaded into SIMD vectors.
   */
  struct LJParameters {
    /**
     * Vector of sigma squared values.
     */
    VectorDouble sigmaSquared;
    /**
     * Vector of 24 * epsilon values.
     */
    VectorDouble epsilon24;
    /**
     * Vector of shift6 values.
     */
    VectorDouble shift6;
  };

  /**
   * Computes shift6 for SIMD vectors of sigmaSquared and epsilon24 using the inverse cutoff squared.
   *
   * @param sigmaSquared Vector of sigma squared values.
   * @param epsilon24 Vector of 24 * epsilon values.
   * @param invCutoffSquared Vector containing 1 / (cutoff^2).
   * @return Vector of shift6 values.
   */
  static inline VectorDouble computeShift6SIMD(const VectorDouble &sigmaSquared, const VectorDouble &epsilon24,
                                               const VectorDouble &invCutoffSquared) {
    const auto sig2Cut2 = highway::Mul(sigmaSquared, invCutoffSquared);
    const auto sig4Cut4 = highway::Mul(sig2Cut2, sig2Cut2);
    const auto sig6Cut6 = highway::Mul(sig4Cut4, sig2Cut2);
    const auto sig12Cut12 = highway::Mul(sig6Cut6, sig6Cut6);
    return highway::Mul(epsilon24, highway::Sub(sig6Cut6, sig12Cut12));
  }

  /**
   * Helper function that handles an iteration of the outer-loop (the loop that changes "i").
   * @tparam reversed Whether iterating backwards over the outer-loop.
   * @tparam newton3
   * @tparam remainderI True iff the "i" vector register is not fully filled, i.e. there are less than _vecLength
   * elements of the SoA to iterate over.
   * @tparam vecPattern
   * @param i index of SoA 1 from which the vector register is filled. E.g. with reversed=false, and the VecLength x 1
   * pattern used, the vector register will be filled with elements i, i+1, i+2, ..., i+_vecLength-1. With
   * reversed=true, this is instead i-_vecLength+1, i-_vecLength+2, ..., i-1, i.
   * @param xPtr1
   * @param yPtr1
   * @param zPtr1
   * @param ownedStatePtr1
   * @param xPtr2
   * @param yPtr2
   * @param zPtr2
   * @param ownedStatePtr2
   * @param fxPtr1
   * @param fyPtr1
   * @param fzPtr1
   * @param fxPtr2
   * @param fyPtr2
   * @param fzPtr2
   * @param halfSigma1Ptr
   * @param sqrtEpsilon1Ptr
   * @param halfSigma2Ptr
   * @param sqrtEpsilon2Ptr
   * @param virialSumX
   * @param virialSumY
   * @param virialSumZ
   * @param potentialEnergySum
   * @param restI In the remainder case, this is the number of elements left in the SoA.
   * @param jVecStart The start index for the inner "j" loop.
   * @param jVecEnd The end index for the inner "j" loop.
   * @param homoSigmaSquared
   * @param homoEpsilon24
   * @param homoShift6
   */
  template <bool reversed, bool newton3, bool remainderI, VectorizationPattern vecPattern>
  inline void handleILoopBody(
      const size_t i, const double *const __restrict xPtr1, const double *const __restrict yPtr1,
      const double *const __restrict zPtr1, const autopas::OwnershipState *const __restrict ownedStatePtr1,
      const double *const __restrict xPtr2, const double *const __restrict yPtr2, const double *const __restrict zPtr2,
      const autopas::OwnershipState *const __restrict ownedStatePtr2, double *const __restrict fxPtr1,
      double *const __restrict fyPtr1, double *const __restrict fzPtr1, double *const __restrict fxPtr2,
      double *const __restrict fyPtr2, double *const __restrict fzPtr2, const double *const __restrict halfSigma1Ptr,
      const double *const __restrict sqrtEpsilon1Ptr, const double *const __restrict halfSigma2Ptr,
      const double *const __restrict sqrtEpsilon2Ptr, VectorDouble &virialSumX, VectorDouble &virialSumY,
      VectorDouble &virialSumZ, VectorDouble &potentialEnergySum, const size_t restI, const size_t jVecStart,
      const size_t jVecEnd, const double homoSigmaSquared = 0.0, const double homoEpsilon24 = 0.0,
      const double homoShift6 = 0.0) const {
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
    const auto *const hs1Ptr = halfSigma1Ptr != nullptr ? &halfSigma1Ptr[i] : nullptr;
    const auto *const se1Ptr = sqrtEpsilon1Ptr != nullptr ? &sqrtEpsilon1Ptr[i] : nullptr;

    auto j = static_cast<std::ptrdiff_t>(jVecStart);
    for (; checkSecondLoopCondition<vecPattern>(jVecEnd, j);
         j += static_cast<std::ptrdiff_t>(jStepSize<vecPattern>())) {
      const auto *const hs2Ptr = halfSigma2Ptr != nullptr ? &halfSigma2Ptr[j] : nullptr;
      const auto *const se2Ptr = sqrtEpsilon2Ptr != nullptr ? &sqrtEpsilon2Ptr[j] : nullptr;
      SoAKernel<newton3, remainderI, false, reversed, vecPattern>(
          j, notDummyMaskI, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2,
          zPtr2, fxPtr2, fyPtr2, fzPtr2, hs1Ptr, se1Ptr, hs2Ptr, se2Ptr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
          virialSumZ, potentialEnergySum, restI, 0, homoSigmaSquared, homoEpsilon24, homoShift6);
    }

    const size_t restJ = jVecEnd & (jStepSize<vecPattern>() - 1);
    if (restJ > 0) {
      const auto *const hs2Ptr = halfSigma2Ptr != nullptr ? &halfSigma2Ptr[j] : nullptr;
      const auto *const se2Ptr = sqrtEpsilon2Ptr != nullptr ? &sqrtEpsilon2Ptr[j] : nullptr;
      SoAKernel<newton3, remainderI, true, reversed, vecPattern>(
          j, notDummyMaskI, ownedMaskI, reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2,
          zPtr2, fxPtr2, fyPtr2, fzPtr2, hs1Ptr, se1Ptr, hs2Ptr, se2Ptr, fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
          virialSumZ, potentialEnergySum, restI, restJ, homoSigmaSquared, homoEpsilon24, homoShift6);
    }

    reduceAccumulatedForce<reversed, remainderI, vecPattern>(i, fxPtr1, fyPtr1, fzPtr1, fxAcc, fyAcc, fzAcc, restI);
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
    if (soa1.size() == 0 || soa2.size() == 0) {
      return;
    }

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
    MixingBuffers mixingBuffers1;
    MixingBuffers mixingBuffers2;
    double *halfSigma1Ptr = nullptr;
    double *sqrtEpsilon1Ptr = nullptr;
    double *halfSigma2Ptr = nullptr;
    double *sqrtEpsilon2Ptr = nullptr;
    double homoSigmaSquared = _sigmaSquared;
    double homoEpsilon24 = _epsilon24;
    double homoShift6 = _shift6;

    if constexpr (useMixing) {
      const auto *const __restrict typeID1Ptr = soa1.template begin<Particle_T::AttributeNames::typeId>();
      const auto *const __restrict typeID2Ptr = soa2.template begin<Particle_T::AttributeNames::typeId>();
      if (soaIsHomogeneous(typeID1Ptr, n1) && soaIsHomogeneous(typeID2Ptr, n2)) {
        const size_t type1 = n1 > 0 ? typeID1Ptr[0] : 0;
        const size_t type2 = n2 > 0 ? typeID2Ptr[0] : 0;
        const double sigma = _halfSigmaPerType[type1] + _halfSigmaPerType[type2];
        homoSigmaSquared = sigma * sigma;
        homoEpsilon24 = 24.0 * _sqrtEpsilonPerType[type1] * _sqrtEpsilonPerType[type2];
        if constexpr (calculateGlobals and applyShift) {
          homoShift6 =
              ParticlePropertiesLibrary<double, size_t>::calcShift6(homoEpsilon24, homoSigmaSquared, _cutoffSquared);
        }
      } else {
        mixingBuffers1.template getMixingPointers<0>(typeID1Ptr, n1, _halfSigmaPerType, _sqrtEpsilonPerType,
                                                     halfSigma1Ptr, sqrtEpsilon1Ptr);
        mixingBuffers2.template getMixingPointers<1>(typeID2Ptr, n2, _halfSigmaPerType, _sqrtEpsilonPerType,
                                                     halfSigma2Ptr, sqrtEpsilon2Ptr);
      }
    }

    const std::ptrdiff_t startI =
        sorted && sortingData.has_value() ? static_cast<std::ptrdiff_t>(sortingData->get().startI) : 0;
    // endI bounds the outer loop from above, mirroring startI: i-particles from endI onwards cannot interact
    // with any j-particle, so the loop stops there instead of running to n1 and relying on the per-block
    // jVecStart >= jVecEnd skip below.
    const std::ptrdiff_t endI = sorted && sortingData.has_value() ? static_cast<std::ptrdiff_t>(sortingData->get().endI)
                                                                  : static_cast<std::ptrdiff_t>(n1);

    VectorDouble virialSumX = highway::Zero(tag_double);
    VectorDouble virialSumY = highway::Zero(tag_double);
    VectorDouble virialSumZ = highway::Zero(tag_double);
    VectorDouble potentialEnergySum = highway::Zero(tag_double);

    const size_t iStep = iStepSize<vecPattern>();
    const size_t jStep = jStepSize<vecPattern>();

    std::ptrdiff_t i = startI;
    for (; i + static_cast<std::ptrdiff_t>(iStep) <= endI; i += static_cast<std::ptrdiff_t>(iStep)) {
      size_t jVecEnd{};
      size_t jVecStart = 0;
      if constexpr (sorted) {
        // The get is always safe here since SoAFunctorPairSorted() will never call this with no sortingData.
        const auto &sd = sortingData->get();
        // maxIndex is monotonically non-decreasing, so the tightest valid bound for [i, i + iStep - 1] (the i values
        // handled in one iteration) is maxIndex of the last particle in the block. For p1xVec
        // (iStep=1) this collapses to maxIndex[i].
        jVecEnd = sd.maxIndex[i + iStep - 1];
        // minIndex is monotonically non-decreasing, so the tightest valid lower bound [i, i + iStep - 1] (the i values
        // handled in one iteration) is always minIndex[i], the minimum across the block.
        jVecStart = sd.minIndex[i];
        // If this check is true it means there are no particles in soa2 that can interact with particles in soa1
        // I.e. the hitrate in this case is 0%.
        if (jVecStart >= jVecEnd) {
          continue;
        }
        // Round down to the nearest full SIMD lane boundary so the j-loop starts aligned.
        jVecStart = jVecStart - (jVecStart % jStep);
        if constexpr (vecPattern == VectorizationPattern::p1xVec) {
          if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
            continue;
          }
        }
      } else {
        jVecEnd = n2;
      }
      handleILoopBody<false, newton3, false, vecPattern>(
          i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx1Ptr, fy1Ptr, fz1Ptr, fx2Ptr,
          fy2Ptr, fz2Ptr, halfSigma1Ptr, sqrtEpsilon1Ptr, halfSigma2Ptr, sqrtEpsilon2Ptr, virialSumX, virialSumY,
          virialSumZ, potentialEnergySum, 0, jVecStart, jVecEnd, homoSigmaSquared, homoEpsilon24, homoShift6);
    }
    if constexpr (vecPattern != VectorizationPattern::p1xVec) {
      // Rest I can't occur in 1xVec case. Bounded by endI, not n1: i-particles from endI onwards cannot
      // interact with any j-particle and are skipped entirely rather than processed as a no-op remainder.
      const std::ptrdiff_t restI = endI - i;
      if (restI > 0) {
        // Remainder block covers [i, i + restI - 1]. Same monotonicity argument as above.
        size_t jVecEnd = n2;
        size_t jVecStart = 0;
        if constexpr (sorted) {
          // The get is always safe here since SoAFunctorPairSorted() will never call this with no sortingData.
          const auto &sd = sortingData->get();
          jVecEnd = sd.maxIndex[i + restI - 1];
          jVecStart = sd.minIndex[i];
          if (jVecStart < jVecEnd) {
            // Round down to nearest full SIMD lane boundary.
            jVecStart = jVecStart - (jVecStart % jStep);
          }
        }
        // If this check is false it means there are no particles in soa2 that can interact with particles in soa1
        // I.e. the hitrate in this case is 0%.
        if (jVecStart < jVecEnd) {
          handleILoopBody<false, newton3, true, vecPattern>(
              i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx1Ptr, fy1Ptr, fz1Ptr,
              fx2Ptr, fy2Ptr, fz2Ptr, halfSigma1Ptr, sqrtEpsilon1Ptr, halfSigma2Ptr, sqrtEpsilon2Ptr, virialSumX,
              virialSumY, virialSumZ, potentialEnergySum, static_cast<size_t>(restI), jVecStart, jVecEnd,
              homoSigmaSquared, homoEpsilon24, homoShift6);
        }
      }
    }

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, potentialEnergySum);
    }
  }

  template <bool remainderI, bool remainderJ, bool reversed, VectorizationPattern vecPattern>
  inline void fillPhysicsRegisters(const double *const halfSigma1Ptr, const double *const sqrtEpsilon1Ptr,
                                   const double *const halfSigma2Ptr, const double *const sqrtEpsilon2Ptr,
                                   VectorDouble &epsilon24s, VectorDouble &sigmaSquareds, VectorDouble &shift6s,
                                   const unsigned int restI, const unsigned int restJ) const {
    // We overestimate the array size. This should be a tight/perfect upper bound on x86, and should never be large
    // enough on e.g. ARM/RISC-V to cause issues.
    HWY_ALIGN std::array<double, _maxVecLengthDouble> epsilons{};
    HWY_ALIGN std::array<double, _maxVecLengthDouble> sigmas{};

    if constexpr (vecPattern == VectorizationPattern::p1xVec) {
      const auto hs1 = *halfSigma1Ptr;
      const auto se1 = *sqrtEpsilon1Ptr;
      for (int j = 0; j < (remainderJ ? restJ : _vecLengthDouble); ++j) {
        const auto hs2 = *(halfSigma2Ptr + j);
        const auto se2 = *(sqrtEpsilon2Ptr + j);
        const auto sigma = hs1 + hs2;
        const auto sigmaSquared = sigma * sigma;
        const auto epsilon24 = 24.0 * se1 * se2;
        epsilons[j] = epsilon24;
        sigmas[j] = sigmaSquared;
      }
    } else if constexpr (vecPattern == VectorizationPattern::p2xVecDiv2) {
      for (int i = 0; i < (remainderI ? 1 : 2); ++i) {
        const auto hs1 = *(reversed ? halfSigma1Ptr - i : halfSigma1Ptr + i);
        const auto se1 = *(reversed ? sqrtEpsilon1Ptr - i : sqrtEpsilon1Ptr + i);
        for (int j = 0; j < (remainderJ ? restJ : _vecLengthDouble / 2); ++j) {
          const auto index = i * (_vecLengthDouble / 2) + j;
          const auto hs2 = *(halfSigma2Ptr + j);
          const auto se2 = *(sqrtEpsilon2Ptr + j);
          const auto sigma = hs1 + hs2;
          const auto sigmaSquared = sigma * sigma;
          const auto epsilon24 = 24.0 * se1 * se2;
          epsilons[index] = epsilon24;
          sigmas[index] = sigmaSquared;
        }
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecDiv2x2) {
      for (int i = 0; i < (remainderI ? restI : _vecLengthDouble / 2); ++i) {
        const auto hs1 = *(reversed ? halfSigma1Ptr - i : halfSigma1Ptr + i);
        const auto se1 = *(reversed ? sqrtEpsilon1Ptr - i : sqrtEpsilon1Ptr + i);
        for (int j = 0; j < (remainderJ ? 1 : 2); ++j) {
          const auto index = i + _vecLengthDouble / 2 * j;
          const auto hs2 = *(halfSigma2Ptr + j);
          const auto se2 = *(sqrtEpsilon2Ptr + j);
          const auto sigma = hs1 + hs2;
          const auto sigmaSquared = sigma * sigma;
          const auto epsilon24 = 24.0 * se1 * se2;

          epsilons[index] = epsilon24;
          sigmas[index] = sigmaSquared;
        }
      }
    } else if constexpr (vecPattern == VectorizationPattern::pVecx1) {
      for (int i = 0; i < (remainderI ? restI : _vecLengthDouble); ++i) {
        const auto hs1 = *(reversed ? halfSigma1Ptr - i : halfSigma1Ptr + i);
        const auto se1 = *(reversed ? sqrtEpsilon1Ptr - i : sqrtEpsilon1Ptr + i);
        const auto hs2 = *halfSigma2Ptr;
        const auto se2 = *sqrtEpsilon2Ptr;
        const auto sigma = hs1 + hs2;
        const auto sigmaSquared = sigma * sigma;
        const auto epsilon24 = 24.0 * se1 * se2;
        epsilons[i] = epsilon24;
        sigmas[i] = sigmaSquared;
      }
    }

    epsilon24s = highway::Load(tag_double, epsilons.data());
    sigmaSquareds = highway::Load(tag_double, sigmas.data());
    if constexpr (calculateGlobals and applyShift) {
      shift6s = computeShift6SIMD(sigmaSquareds, epsilon24s, highway::Set(tag_double, _invCutoffSquared));
    } else {
      shift6s = highway::Zero(tag_double);
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
   * @param j
   * @param notDummyMaskI
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
   * @param halfSigma1Ptr
   * @param sqrtEpsilon1Ptr
   * @param halfSigma2Ptr
   * @param sqrtEpsilon2Ptr
   * @param fxAcc
   * @param fyAcc
   * @param fzAcc
   * @param virialSumX
   * @param virialSumY
   * @param virialSumZ
   * @param potentialEnergySum
   * @param restI
   * @param restJ
   * @param homoSigmaSquared
   * @param homoEpsilon24
   * @param homoShift6
   */
  template <bool newton3, bool remainderI, bool remainderJ, bool reversed, VectorizationPattern vecPattern>
  inline void SoAKernel(const size_t j, const MaskDouble &notDummyMaskI, [[maybe_unused]] const MaskDouble &ownedMaskI,
                        const int64_t *const __restrict ownedStatePtr2, const VectorDouble &x1, const VectorDouble &y1,
                        const VectorDouble &z1, const double *const __restrict x2Ptr,
                        const double *const __restrict y2Ptr, const double *const __restrict z2Ptr,
                        double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                        double *const __restrict fz2Ptr, const double *const halfSigma1Ptr,
                        const double *const sqrtEpsilon1Ptr, const double *const halfSigma2Ptr,
                        const double *const sqrtEpsilon2Ptr, VectorDouble &fxAcc, VectorDouble &fyAcc,
                        VectorDouble &fzAcc, VectorDouble &virialSumX, VectorDouble &virialSumY,
                        VectorDouble &virialSumZ, VectorDouble &potentialEnergySum, const unsigned int restI,
                        const unsigned int restJ, const double homoSigmaSquared = 0.0, const double homoEpsilon24 = 0.0,
                        const double homoShift6 = 0.0) const {
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

    const auto drX2 = highway::Mul(drX, drX);
    const auto drY2 = highway::Mul(drY, drY);
    const auto drZ2 = highway::Mul(drZ, drZ);

    const auto dr2 = highway::Add(highway::Add(drX2, drY2), drZ2);

    VectorDouble cutoffSquared = highway::Set(tag_double, _cutoffSquared);

    const auto dummyMask = highway::And(notDummyMaskI, notDummyMaskJ);
    const auto cutoffDummyMask = highway::MaskedLe(dummyMask, dr2, cutoffSquared);

    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
      return;
    }

    VectorDouble epsilon24s = highway::Undefined(tag_double);
    VectorDouble sigmaSquareds = highway::Undefined(tag_double);
    VectorDouble shift6s = highway::Undefined(tag_double);

    if constexpr (useMixing) {
      if (halfSigma1Ptr != nullptr) {
        fillPhysicsRegisters<remainderI, remainderJ, reversed, vecPattern>(halfSigma1Ptr, sqrtEpsilon1Ptr,
                                                                           halfSigma2Ptr, sqrtEpsilon2Ptr, epsilon24s,
                                                                           sigmaSquareds, shift6s, restI, restJ);
      } else {
        epsilon24s = highway::Set(tag_double, homoEpsilon24);
        sigmaSquareds = highway::Set(tag_double, homoSigmaSquared);
        if constexpr (calculateGlobals and applyShift) {
          shift6s = highway::Set(tag_double, homoShift6);
        } else {
          shift6s = highway::Zero(tag_double);
        }
      }
    } else {
      epsilon24s = highway::Set(tag_double, _epsilon24);
      sigmaSquareds = highway::Set(tag_double, _sigmaSquared);
      if constexpr (calculateGlobals and applyShift) {
        shift6s = highway::Set(tag_double, _shift6);
      } else {
        shift6s = highway::Zero(tag_double);
      }
    }

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

      auto potentialEnergy = highway::MulAdd(epsilon24s, lj12m6, shift6s);
      auto potentialEnergyMasked = highway::IfThenElseZero(cutoffDummyMask, potentialEnergy);

      auto energyFactor = highway::MaskedSet(tag_double, ownedMaskI, 1.0);

      if constexpr (newton3) {
        energyFactor = highway::Add(energyFactor, highway::MaskedSet(tag_double, ownedMaskJ, 1.0));
      }

      potentialEnergySum = highway::MulAdd(energyFactor, potentialEnergyMasked, potentialEnergySum);
      virialSumX = highway::MulAdd(energyFactor, virialX, virialSumX);
      virialSumY = highway::MulAdd(energyFactor, virialY, virialSumY);
      virialSumZ = highway::MulAdd(energyFactor, virialZ, virialSumZ);
    }
  }

  template <bool newton3, bool remainder = false>
  inline void SoAKernelVerlet(const size_t j, const MaskDouble &notDummyMaskI,
                              [[maybe_unused]] const MaskDouble &ownedMaskI, const VectorDouble &x1,
                              const VectorDouble &y1, const VectorDouble &z1, const double *const __restrict xPtr,
                              const double *const __restrict yPtr, const double *const __restrict zPtr,
                              const int64_t *const __restrict ownedStatePtr, double *const __restrict fxPtr,
                              double *const __restrict fyPtr, double *const __restrict fzPtr,
                              const double *const halfSigma1Ptr, const double *const sqrtEpsilon1Ptr,
                              const size_t *const __restrict typeIDPtr, const size_t *const __restrict neighborList,
                              VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc, VectorDouble &virialSumX,
                              VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &potentialEnergySum,
                              [[maybe_unused]] const size_t rest = 0, const double homoSigmaSquared = 0.0,
                              const double homoEpsilon24 = 0.0, const double homoShift6 = 0.0) const {
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

    const auto drX2 = highway::Mul(drX, drX);
    const auto drY2 = highway::Mul(drY, drY);
    const auto drZ2 = highway::Mul(drZ, drZ);

    const auto dr2 = highway::Add(highway::Add(drX2, drY2), drZ2);

    VectorDouble cutoffSquared = highway::Set(tag_double, _cutoffSquared);

    const auto dummyMask = highway::And(notDummyMaskI, notDummyMaskJ);
    const auto cutoffDummyMask = highway::MaskedLe(dummyMask, dr2, cutoffSquared);

    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
      return;
    }

    VectorDouble epsilon24s = highway::Undefined(tag_double);
    VectorDouble sigmaSquareds = highway::Undefined(tag_double);
    VectorDouble shift6s = highway::Undefined(tag_double);

    if constexpr (useMixing) {
      if (halfSigma1Ptr != nullptr) {
        HWY_ALIGN std::array<double, _maxVecLengthDouble> epsilons{};
        HWY_ALIGN std::array<double, _maxVecLengthDouble> sigmas{};
        const size_t kMax = remainder ? rest : _vecLengthDouble;
        const auto hs1 = *halfSigma1Ptr;
        const auto se1 = *sqrtEpsilon1Ptr;
        for (size_t k = 0; k < kMax; ++k) {
          const auto neighborIdx = neighborList[j + k];
          const auto neighborType = typeIDPtr[neighborIdx];
          const auto hs2 = _halfSigmaPerType[neighborType];
          const auto se2 = _sqrtEpsilonPerType[neighborType];
          const auto sigma = hs1 + hs2;
          const auto sigmaSquared = sigma * sigma;
          const auto epsilon24 = 24.0 * se1 * se2;
          epsilons[k] = epsilon24;
          sigmas[k] = sigmaSquared;
        }
        epsilon24s = highway::Load(tag_double, epsilons.data());
        sigmaSquareds = highway::Load(tag_double, sigmas.data());
        if constexpr (calculateGlobals and applyShift) {
          shift6s = computeShift6SIMD(sigmaSquareds, epsilon24s, highway::Set(tag_double, _invCutoffSquared));
        } else {
          shift6s = highway::Zero(tag_double);
        }
      } else {
        epsilon24s = highway::Set(tag_double, homoEpsilon24);
        sigmaSquareds = highway::Set(tag_double, homoSigmaSquared);
        if constexpr (calculateGlobals and applyShift) {
          shift6s = highway::Set(tag_double, homoShift6);
        } else {
          shift6s = highway::Zero(tag_double);
        }
      }
    } else {
      epsilon24s = highway::Set(tag_double, _epsilon24);
      sigmaSquareds = highway::Set(tag_double, _sigmaSquared);
      if constexpr (calculateGlobals and applyShift) {
        shift6s = highway::Set(tag_double, _shift6);
      } else {
        shift6s = highway::Zero(tag_double);
      }
    }

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

      auto potentialEnergy = highway::MulAdd(epsilon24s, lj12m6, shift6s);
      auto potentialEnergyMasked = highway::IfThenElseZero(cutoffDummyMask, potentialEnergy);

      auto energyFactor = highway::MaskedSet(tag_double, ownedMaskI, 1.0);

      if constexpr (newton3) {
        const VectorLong ownedStateVal = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
        const MaskDouble ownedMaskJ = highway::RebindMask(tag_double, highway::Eq(ownedState2, ownedStateVal));
        energyFactor = highway::Add(energyFactor, highway::MaskedSet(tag_double, ownedMaskJ, 1.0));
      }

      potentialEnergySum = highway::MulAdd(energyFactor, potentialEnergyMasked, potentialEnergySum);
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

    double homoSigmaSquared = _sigmaSquared;
    double homoEpsilon24 = _epsilon24;
    double homoShift6 = _shift6;
    double halfSigmaI = 0.0;
    double sqrtEpsilonI = 0.0;
    const double *hs1Ptr = nullptr;
    const double *se1Ptr = nullptr;
    const size_t *typeIDPtr = nullptr;

    const bool isHomogeneous = not useMixing or _halfSigmaPerType.size() <= 1;

    if constexpr (useMixing) {
      typeIDPtr = soa.template begin<Particle_T::AttributeNames::typeId>();
      if (not isHomogeneous) {
        const size_t typeI = typeIDPtr[indexFirst];
        halfSigmaI = _halfSigmaPerType[typeI];
        sqrtEpsilonI = _sqrtEpsilonPerType[typeI];
        hs1Ptr = &halfSigmaI;
        se1Ptr = &sqrtEpsilonI;
      }
    }

    VectorDouble virialSumX = highway::Zero(tag_double);
    VectorDouble virialSumY = highway::Zero(tag_double);
    VectorDouble virialSumZ = highway::Zero(tag_double);
    VectorDouble potentialEnergySum = highway::Zero(tag_double);
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
      constexpr auto ownedVal = static_cast<double>(static_cast<int64_t>(autopas::OwnershipState::owned));
      ownedMaskI = highway::Eq(ownedStateI, highway::Set(tag_double, ownedVal));
    }

    size_t j = 0;
    const size_t neighborListSize = neighborList.size();
    const size_t vecEnd = (neighborListSize / _vecLengthDouble) * _vecLengthDouble;

    for (; j < vecEnd; j += _vecLengthDouble) {
      SoAKernelVerlet<newton3, false>(
          j, notDummyMaskI, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, reinterpret_cast<const int64_t *>(ownedStatePtr),
          fxPtr, fyPtr, fzPtr, hs1Ptr, se1Ptr, typeIDPtr, neighborList.data(), fxAcc, fyAcc, fzAcc, virialSumX,
          virialSumY, virialSumZ, potentialEnergySum, 0, homoSigmaSquared, homoEpsilon24, homoShift6);
    }

    const size_t rest = neighborListSize & (_vecLengthDouble - 1);

    if (rest > 0) {
      SoAKernelVerlet<newton3, true>(
          j, notDummyMaskI, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr, reinterpret_cast<const int64_t *>(ownedStatePtr),
          fxPtr, fyPtr, fzPtr, hs1Ptr, se1Ptr, typeIDPtr, neighborList.data(), fxAcc, fyAcc, fzAcc, virialSumX,
          virialSumY, virialSumZ, potentialEnergySum, rest, homoSigmaSquared, homoEpsilon24, homoShift6);
    }

    fxPtr[indexFirst] += highway::ReduceSum(tag_double, fxAcc);
    fyPtr[indexFirst] += highway::ReduceSum(tag_double, fyAcc);
    fzPtr[indexFirst] += highway::ReduceSum(tag_double, fzAcc);

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, potentialEnergySum);
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
    for (size_t i = 0; i < _threadData.size(); ++i) {
      _threadData[i].setZero();
    }
  }

  /**
   * Accumulates global values, e.g. potential energy and virial.
   * @param newton3
   */
  void endTraversal(const bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (_postProcessed) {
      throw autopas::utils::ExceptionHandler::AutoPasException(
          "Already postprocessed, endTraversal(bool newton3) was called twice without calling initTraversal().");
    }

    if (calculateGlobals) {
      for (size_t i = 0; i < _threadData.size(); ++i) {
        _potentialEnergySum += _threadData[i].potentialEnergySum;
        _virialSum += _threadData[i].virialSum;
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
   * @param sigmaSquared
   */
  void setParticleProperties(const double epsilon24, const double sigmaSquared) {
    _epsilon24 = epsilon24;
    _sigmaSquared = sigmaSquared;
    _halfSigmaPerType = {std::sqrt(sigmaSquared) * 0.5};
    _sqrtEpsilonPerType = {std::sqrt(epsilon24 / 24.0)};
    if constexpr (applyShift) {
      _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquared);
    } else {
      _shift6 = 0.;
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
  class ThreadData {
   public:
    ThreadData() : virialSum{0., 0., 0.}, potentialEnergySum{0.} {}

    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }

    std::array<double, 3> virialSum{};
    double potentialEnergySum{};

   private:
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };
  static_assert(sizeof(ThreadData) % 64 == 0, "ThreadData has wrong size");

  // cutoff squared used in the functor.
  const double _cutoffSquared{0.};
  // inverse cutoff squared used in the functor.
  const double _invCutoffSquared{0.};
  // epsilon, sigma and shift6 used in the functor.
  double _epsilon24{0.}, _sigmaSquared{0.}, _shift6{0.};
  std::vector<double> _halfSigmaPerType{};
  std::vector<double> _sqrtEpsilonPerType{};

  // accumulators for the globals (potential energy and virial).
  double _potentialEnergySum{0.};
  std::array<double, 3> _virialSum{0., 0., 0.};
  std::vector<ThreadData> _threadData{};

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
