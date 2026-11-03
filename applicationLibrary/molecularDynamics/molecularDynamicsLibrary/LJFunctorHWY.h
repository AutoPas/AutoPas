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
        _invCutoffSquared{1. / (cutoff * cutoff)},
        _cutoffSquaredVec{highway::Set(tag_double, cutoff * cutoff)} {
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
   * @copydoc autopas::Functor::setVecPattern()
   */
  void setVecPattern(const VectorizationPattern vecPattern) final { _vecPattern = vecPattern; }

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

    double halfSigI = 0., sqrtEpsI = 0.;
    size_t typeJ = 0;
    if constexpr (useMixing) {
      const auto typeI = i.getTypeId();
      typeJ = j.getTypeId();
      halfSigI = _halfSigmaPerType[typeI];
      sqrtEpsI = _sqrtEpsilonPerType[typeI];
    }

    const auto [fx, fy, fz, potentialEnergy6] =
        computeScalarLJInteraction(dr2, dr[0], dr[1], dr[2], halfSigI, sqrtEpsI, typeJ);
    const std::array<double, 3> f{fx, fy, fz};
    i.addF(f);
    if (newton3) {
      j.subF(f);
    }
    if constexpr (calculateGlobals) {
      const auto virial = dr * f;
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
  inline void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, [[maybe_unused]] const bool newton3) final {
    const auto soaSize = soa.size();
    if (soaSize <= 1) return;

    if (soaSize <= 7) {
      const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
      const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
      const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();
      const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();
      auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
      auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
      auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();
      const auto *const __restrict typeIDPtr =
          useMixing ? soa.template begin<Particle_T::AttributeNames::typeId>() : nullptr;

      scalarSoAKernel<true, true>(soaSize, soaSize, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, typeIDPtr,
                                  xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, typeIDPtr);
      return;
    }

    if constexpr (useMixing) {
      const auto *const __restrict typeIDPtr = soa.template begin<Particle_T::AttributeNames::typeId>();
      if (soaIsHomogeneous(typeIDPtr, soa.size())) {
        const size_t type = soa.size() > 0 ? typeIDPtr[0] : 0;
        auto paramPolicy = makeHomogeneousParameterPolicy(type, type);
        runSingleLoop(soa, paramPolicy);
      } else {
        MixingBuffers mixingBuffers;
        double *halfSigmaPtr = nullptr;
        double *sqrtEpsilonPtr = nullptr;
        mixingBuffers.getMixingPointers(typeIDPtr, soa.size(), _halfSigmaPerType, _sqrtEpsilonPerType, halfSigmaPtr,
                                        sqrtEpsilonPtr);
        HeterogeneousParameterPolicy paramPolicy{halfSigmaPtr, sqrtEpsilonPtr, halfSigmaPtr, sqrtEpsilonPtr,
                                                 _invCutoffSquared};
        runSingleLoop(soa, paramPolicy);
      }
    } else {
      auto paramPolicy = makeHomogeneousParameterPolicy();
      runSingleLoop(soa, paramPolicy);
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

  /**
   * Struct holding Lennard-Jones parameters loaded into SIMD vectors for an interaction.
   */
  struct LJPairParameters {
    VectorDouble sigmaSquared;
    VectorDouble epsilon24;
    VectorDouble shift6;
  };

  /**
   * Parameter policy for homogeneous cells (or interactions with mixing disabled), where all particles share the same
   * physical parameters.
   */
  struct HomogeneousParameterPolicy {
    const LJPairParameters params;

    HomogeneousParameterPolicy(const double sigmaSquaredConst, const double epsilon24Const, const double shift6Const)
        : params{highway::Set(tag_double, sigmaSquaredConst), highway::Set(tag_double, epsilon24Const),
                 highway::Set(tag_double, shift6Const)} {}

    template <bool remainderI, bool reversed, VectorizationPattern vecPattern>
    HWY_INLINE void prepareI([[maybe_unused]] const size_t i, [[maybe_unused]] const size_t restI) const {}

    template <bool remainderJ, VectorizationPattern vecPattern>
    HWY_INLINE const LJPairParameters &getJ([[maybe_unused]] const size_t j,
                                            [[maybe_unused]] const size_t restJ) const {
      return params;
    }

    template <bool remainder, VectorizationPattern vecPattern = VectorizationPattern::p1xVec>
    HWY_INLINE const LJPairParameters &getJ([[maybe_unused]] const VectorLong &indices,
                                            [[maybe_unused]] const MaskDouble &restMaskDouble) const {
      return params;
    }
  };

  /**
   * Parameter policy for heterogeneous cells where multiple particle types interact.
   * Handles per-particle loading and SIMD mixing for cell pairs, single cells, and Verlet lists.
   */
  struct HeterogeneousParameterPolicy {
    const double *const __restrict halfSigma1Ptr{nullptr};
    const double *const __restrict sqrtEpsilon1Ptr{nullptr};
    const double *const __restrict halfSigma2Ptr{nullptr};
    const double *const __restrict sqrtEpsilon2Ptr{nullptr};
    const int64_t *const __restrict typeID2Ptr{nullptr};
    const VectorDouble invCutoffSquaredVec;

    VectorDouble halfSigmaI{highway::Zero(tag_double)};
    VectorDouble sqrtEpsilon24I{highway::Zero(tag_double)};

    /**
     * Constructor for cell-pair and single-cell interactions using pre-unpacked mixing buffers.
     *
     * @param halfSigma1 Pointer to halfSigma buffer for particle(s) 1.
     * @param sqrtEpsilon1 Pointer to sqrtEpsilon buffer for particle(s) 1.
     * @param halfSigma2 Pointer to halfSigma buffer for particle(s) 2.
     * @param sqrtEpsilon2 Pointer to sqrtEpsilon buffer for particle(s) 2.
     * @param invCutoffSquared Scalar 1 / (cutoff^2).
     */
    HeterogeneousParameterPolicy(const double *const __restrict halfSigma1, const double *const __restrict sqrtEpsilon1,
                                 const double *const __restrict halfSigma2, const double *const __restrict sqrtEpsilon2,
                                 const double invCutoffSquared)
        : halfSigma1Ptr{halfSigma1},
          sqrtEpsilon1Ptr{sqrtEpsilon1},
          halfSigma2Ptr{halfSigma2},
          sqrtEpsilon2Ptr{sqrtEpsilon2},
          invCutoffSquaredVec{highway::Set(tag_double, invCutoffSquared)} {}

    /**
     * Constructor for Verlet list interactions using per-type lookup tables and SIMD gather.
     *
     * @param halfSigmaIScalar Scalar halfSigma of particle i.
     * @param sqrtEpsilonIScalar Scalar sqrtEpsilon of particle i.
     * @param typeIDPtr Pointer to particle type IDs in the SoA.
     * @param halfSigmaPerType Pointer to per-type halfSigma lookup array.
     * @param sqrtEpsilonPerType Pointer to per-type sqrtEpsilon lookup array.
     * @param invCutoffSquared Scalar 1 / (cutoff^2).
     */
    HeterogeneousParameterPolicy(const double halfSigmaIScalar, const double sqrtEpsilonIScalar,
                                 const int64_t *const __restrict typeIDPtr,
                                 const double *const __restrict halfSigmaPerType,
                                 const double *const __restrict sqrtEpsilonPerType, const double invCutoffSquared)
        : halfSigma2Ptr{halfSigmaPerType},
          sqrtEpsilon2Ptr{sqrtEpsilonPerType},
          typeID2Ptr{typeIDPtr},
          invCutoffSquaredVec{highway::Set(tag_double, invCutoffSquared)},
          halfSigmaI{highway::Set(tag_double, halfSigmaIScalar)},
          sqrtEpsilon24I{highway::Set(tag_double, 24.0 * sqrtEpsilonIScalar)} {}

    /**
     * Prepares physical parameters for particle(s) i.
     *
     * @tparam remainderI Whether this is a remainder iteration.
     * @tparam reversed Whether outer loop iterates backwards.
     * @tparam vecPattern The vectorization pattern.
     * @param i Loop index i.
     * @param restI Number of valid lanes in remainder case.
     */
    template <bool remainderI, bool reversed, VectorizationPattern vecPattern>
    HWY_INLINE void prepareI(const size_t i, const size_t restI) {
      fillIPhysicsRegisters<remainderI, reversed, vecPattern>(i, halfSigma1Ptr, sqrtEpsilon1Ptr, halfSigmaI,
                                                              sqrtEpsilon24I, restI);
      sqrtEpsilon24I = highway::Mul(highway::Set(tag_double, 24.0), sqrtEpsilon24I);
    }

    /**
     * Mixes parameters of particle(s) i with particle(s) j into an LJPairParameters struct.
     *
     * @param halfSigmaJ SIMD vector of halfSigma for particle(s) j.
     * @param sqrtEpsilonJ SIMD vector of sqrtEpsilon for particle(s) j.
     * @return Evaluated LJPairParameters for the interaction.
     */
    inline LJPairParameters mixParameters(const VectorDouble &halfSigmaJ, const VectorDouble &sqrtEpsilonJ) const {
      const auto sigma = highway::Add(halfSigmaI, halfSigmaJ);
      const auto sigmaSquared = highway::Mul(sigma, sigma);
      const auto epsilon24 = highway::Mul(sqrtEpsilon24I, sqrtEpsilonJ);
      VectorDouble shift6 = highway::Zero(tag_double);
      if constexpr (calculateGlobals and applyShift) {
        shift6 = computeShift6SIMD(sigmaSquared, epsilon24, invCutoffSquaredVec);
      }
      return {sigmaSquared, epsilon24, shift6};
    }

    /**
     * Retrieves the mixed Lennard-Jones parameters for particle(s) j in regular cell interactions.
     *
     * @tparam remainderJ Whether this is a remainder iteration.
     * @tparam vecPattern The vectorization pattern.
     * @param j Loop index j.
     * @param restJ Number of valid lanes in remainder case.
     * @return Evaluated LJPairParameters.
     */
    template <bool remainderJ, VectorizationPattern vecPattern>
    HWY_INLINE LJPairParameters getJ(const size_t j, const size_t restJ) const {
      VectorDouble halfSigmaJ;
      VectorDouble sqrtEpsilonJ;
      fillJPhysicsRegisters<remainderJ, vecPattern>(j, halfSigma2Ptr, sqrtEpsilon2Ptr, halfSigmaJ, sqrtEpsilonJ, restJ);
      return mixParameters(halfSigmaJ, sqrtEpsilonJ);
    }

    /**
     * Retrieves the mixed Lennard-Jones parameters for particle(s) j in Verlet list interactions via SIMD gather.
     *
     * @tparam remainder Whether this is a remainder iteration.
     * @tparam vecPattern The vectorization pattern.
     * @param indices SIMD vector of neighbor indices.
     * @param restMaskDouble Lane mask for remainder case.
     * @return Evaluated LJPairParameters.
     */
    template <bool remainder, VectorizationPattern vecPattern = VectorizationPattern::p1xVec>
    HWY_INLINE LJPairParameters getJ(const VectorLong &indices, const MaskDouble &restMaskDouble) const {
      VectorLong typeJ;
      VectorDouble halfSigmaJ;
      VectorDouble sqrtEpsilonJ;
      if constexpr (remainder) {
        const MaskLong restMaskLong = highway::RebindMask(tag_long, restMaskDouble);
        typeJ = highway::MaskedGatherIndex(restMaskLong, tag_long, typeID2Ptr, indices);
        halfSigmaJ = highway::MaskedGatherIndex(restMaskDouble, tag_double, halfSigma2Ptr, typeJ);
        sqrtEpsilonJ = highway::MaskedGatherIndex(restMaskDouble, tag_double, sqrtEpsilon2Ptr, typeJ);
      } else {
        typeJ = highway::GatherIndex(tag_long, typeID2Ptr, indices);
        halfSigmaJ = highway::GatherIndex(tag_double, halfSigma2Ptr, typeJ);
        sqrtEpsilonJ = highway::GatherIndex(tag_double, sqrtEpsilon2Ptr, typeJ);
      }
      return mixParameters(halfSigmaJ, sqrtEpsilonJ);
    }
  };

  /**
   * Helper method to construct a HomogeneousParameterPolicy for given particle type IDs.
   *
   * @param typeId1 Type ID of the first particle (used if mixing is enabled).
   * @param typeId2 Type ID of the second particle (used if mixing is enabled).
   * @return Initialized HomogeneousParameterPolicy.
   */
  inline HomogeneousParameterPolicy makeHomogeneousParameterPolicy([[maybe_unused]] const size_t typeId1 = 0,
                                                                   [[maybe_unused]] const size_t typeId2 = 0) const {
    double sigmaSquared = _sigmaSquared;
    double epsilon24 = _epsilon24;
    double shift6 = _shift6;
    if constexpr (useMixing) {
      const double sigma = _halfSigmaPerType[typeId1] + _halfSigmaPerType[typeId2];
      sigmaSquared = sigma * sigma;
      epsilon24 = 24.0 * _sqrtEpsilonPerType[typeId1] * _sqrtEpsilonPerType[typeId2];
      if constexpr (calculateGlobals and applyShift) {
        shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquared);
      }
    }
    return HomogeneousParameterPolicy{sigmaSquared, epsilon24, shift6};
  }

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
   * @param halfSigmaPtr Pointer to the halfSigma array.
   * @param sqrtEpsilonPtr Pointer to the sqrtEpsilon array.
   * @param halfSigmaI Output SIMD vector for halfSigma of particle(s) i.
   * @param sqrtEpsilonI Output SIMD vector for sqrtEpsilon of particle(s) i.
   * @param restI Number of valid elements in case remainder is true.
   */
  template <bool remainder, bool reversed, VectorizationPattern vecPattern>
  static void fillIPhysicsRegisters(const size_t i, const double *const __restrict halfSigmaPtr,
                                    const double *const __restrict sqrtEpsilonPtr, VectorDouble &halfSigmaI,
                                    VectorDouble &sqrtEpsilonI, const size_t restI) {
    halfSigmaI = loadIRegister<remainder, reversed, vecPattern>(tag_double, halfSigmaPtr, i, restI);
    sqrtEpsilonI = loadIRegister<remainder, reversed, vecPattern>(tag_double, sqrtEpsilonPtr, i, restI);
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
  HWY_INLINE static auto loadJRegister(Tag tag, const ValueType *const __restrict dataPtr, const size_t j,
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
   * Loads the physical parameters (halfSigma and sqrtEpsilon) for particle(s) j according to the vectorization pattern.
   *
   * @tparam remainder Whether this is a remainder loop iteration with fewer lanes than the vector width.
   * @tparam vecPattern The vectorization pattern.
   * @param j Particle index in the inner j loop.
   * @param halfSigmaPtr Pointer to the halfSigma array.
   * @param sqrtEpsilonPtr Pointer to the sqrtEpsilon array.
   * @param halfSigmaJ Output SIMD vector for halfSigma of particle(s) j.
   * @param sqrtEpsilonJ Output SIMD vector for sqrtEpsilon of particle(s) j.
   * @param restJ Number of valid elements in case remainder is true.
   */
  template <bool remainder, VectorizationPattern vecPattern>
  static void fillJPhysicsRegisters(const size_t j, const double *const __restrict halfSigmaPtr,
                                    const double *const __restrict sqrtEpsilonPtr, VectorDouble &halfSigmaJ,
                                    VectorDouble &sqrtEpsilonJ, const size_t restJ) {
    halfSigmaJ = loadJRegister<remainder, vecPattern>(tag_double, halfSigmaPtr, j, restJ);
    sqrtEpsilonJ = loadJRegister<remainder, vecPattern>(tag_double, sqrtEpsilonPtr, j, restJ);
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
  HWY_INLINE static void updateForceArrays(const Tag tag, const size_t offset, double *const __restrict fxPtr,
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

  /**
   * Recursively folds a vector in half down to 2 lanes by adding its lower and upper halves.
   *
   * @tparam D Tag type defining the vector properties.
   * @param d Tag descriptor.
   * @param v Input SIMD vector.
   * @return A 2-lane SIMD vector containing the folded sum.
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
   *
   * @tparam D Tag type defining the vector properties.
   * @param d Tag descriptor.
   * @param v Input SIMD vector.
   * @return The horizontal sum of all vector lanes.
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
   *
   * @tparam D Tag type defining the vector properties.
   * @param d Tag descriptor.
   * @param va First input SIMD vector.
   * @param vb Second input SIMD vector.
   * @param outA Output scalar accumulator for va.
   * @param outB Output scalar accumulator for vb.
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
   *
   * @tparam D Tag type defining the vector properties.
   * @param d Tag descriptor.
   * @param va First input SIMD vector.
   * @param vb Second input SIMD vector.
   * @param outA Output scalar accumulator for va.
   * @param outB Output scalar accumulator for vb.
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
   * Reduces 3 vectors (x, y, z) into scalar accumulators by interleaving (x, y) and halving z,
   * adding the results to output scalar accumulators.
   *
   * @tparam D Tag type defining the vector properties.
   * @param d Tag descriptor.
   * @param fx Vector of x force components.
   * @param fy Vector of y force components.
   * @param fz Vector of z force components.
   * @param outX Output scalar accumulator for fx.
   * @param outY Output scalar accumulator for fy.
   * @param outZ Output scalar accumulator for fz.
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
   *
   * @tparam D Tag type defining the vector properties.
   * @param d Tag descriptor.
   * @param fx Vector of x force components.
   * @param fy Vector of y force components.
   * @param fz Vector of z force components.
   * @param outX Output scalar accumulator for fx.
   * @param outY Output scalar accumulator for fy.
   * @param outZ Output scalar accumulator for fz.
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
   *
   * @tparam D Tag type defining the vector properties.
   * @param d Tag descriptor.
   * @param v0 First input SIMD vector.
   * @param v1 Second input SIMD vector.
   * @param v2 Third input SIMD vector.
   * @param v3 Fourth input SIMD vector.
   * @param out0 Output scalar accumulator for v0.
   * @param out1 Output scalar accumulator for v1.
   * @param out2 Output scalar accumulator for v2.
   * @param out3 Output scalar accumulator for v3.
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
   *
   * @tparam DHalf Half-width tag descriptor type.
   * @param dh Half-width tag descriptor.
   * @param lfx Lower half of fx.
   * @param lfy Lower half of fy.
   * @param lfz Lower half of fz.
   * @param ufx Upper half of fx.
   * @param ufy Upper half of fy.
   * @param ufz Upper half of fz.
   * @param fxPtr Pointer to the first of two contiguous force x elements.
   * @param fyPtr Pointer to the first of two contiguous force y elements.
   * @param fzPtr Pointer to the first of two contiguous force z elements.
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
   *
   * @tparam DHalf Half-width tag descriptor type.
   * @param dh Half-width tag descriptor.
   * @param lfx Lower half of fx.
   * @param lfy Lower half of fy.
   * @param lfz Lower half of fz.
   * @param ufx Upper half of fx.
   * @param ufy Upper half of fy.
   * @param ufz Upper half of fz.
   * @param fxPtr Pointer to the first of two contiguous force x elements.
   * @param fyPtr Pointer to the first of two contiguous force y elements.
   * @param fzPtr Pointer to the first of two contiguous force z elements.
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
  static void reduceAccumulatedForce(const size_t i, double *const __restrict fxPtr, double *const __restrict fyPtr,
                                     double *const __restrict fzPtr, const VectorDouble &fxAcc,
                                     const VectorDouble &fyAcc, const VectorDouble &fzAcc, const size_t restI) {
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

    reduceAdd4_Interleaved(tag_double, virialSumX, virialSumY, virialSumZ, potentialEnergySum,
                           _threadData[threadID].virialSum[0], _threadData[threadID].virialSum[1],
                           _threadData[threadID].virialSum[2], _threadData[threadID].potentialEnergySum);
  }

  /**
   * Struct holding the result of a scalar Lennard-Jones interaction calculation.
   */
  struct ScalarLJResult {
    double fx{0.};
    double fy{0.};
    double fz{0.};
    double potentialEnergy6{0.};
  };

  /**
   * Computes the scalar Lennard-Jones interaction between two particles.
   *
   * @param dr2 Distance squared between particles.
   * @param dx Distance vector x-component.
   * @param dy Distance vector y-component.
   * @param dz Distance vector z-component.
   * @param halfSigmaI Pre-looked-up halfSigma for particle i (used when useMixing=true).
   * @param sqrtEpsilonI Pre-looked-up sqrtEpsilon for particle i (used when useMixing=true).
   * @param typeJ Type ID of particle j (used when useMixing=true).
   * @return ScalarLJResult containing forces and potential energy.
   */
  HWY_INLINE ScalarLJResult computeScalarLJInteraction(const double dr2, const double dx, const double dy,
                                                       const double dz, [[maybe_unused]] const double halfSigmaI,
                                                       [[maybe_unused]] const double sqrtEpsilonI,
                                                       [[maybe_unused]] const size_t typeJ) const {
    double sigmaSquared = _sigmaSquared;
    double epsilon24 = _epsilon24;
    [[maybe_unused]] double shift6 = _shift6;
    if constexpr (useMixing) {
      const double sig = halfSigmaI + _halfSigmaPerType[typeJ];
      sigmaSquared = sig * sig;
      epsilon24 = 24.0 * sqrtEpsilonI * _sqrtEpsilonPerType[typeJ];
      if constexpr (calculateGlobals and applyShift) {
        shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(epsilon24, sigmaSquared, _cutoffSquared);
      }
    }
    const double invdr2 = 1.0 / dr2;
    const double lj2 = sigmaSquared * invdr2;
    const double lj6 = lj2 * lj2 * lj2;
    const double lj12 = lj6 * lj6;
    const double lj12m6 = lj12 - lj6;
    const double fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    ScalarLJResult result{dx * fac, dy * fac, dz * fac, 0.};
    if constexpr (calculateGlobals) {
      result.potentialEnergy6 = epsilon24 * lj12m6 + shift6;
    }
    return result;
  }

  /**
   * Scalar fallback kernel for small particle counts (N <= 7).
   *
   * @tparam isSingleCell Whether computing interactions within a single cell (i > j).
   * @tparam newton3 Whether Newton's third law is enabled.
   * @param soaSize1 Number of particles in first cell.
   * @param soaSize2 Number of particles in second cell.
   * @param x1Ptr Pointer to posX array of cell 1.
   * @param y1Ptr Pointer to posY array of cell 1.
   * @param z1Ptr Pointer to posZ array of cell 1.
   * @param ownedStatePtr1 Pointer to ownershipState array of cell 1.
   * @param fx1Ptr Pointer to forceX array of cell 1.
   * @param fy1Ptr Pointer to forceY array of cell 1.
   * @param fz1Ptr Pointer to forceZ array of cell 1.
   * @param typeID1Ptr Pointer to typeId array of cell 1.
   * @param x2Ptr Pointer to posX array of cell 2.
   * @param y2Ptr Pointer to posY array of cell 2.
   * @param z2Ptr Pointer to posZ array of cell 2.
   * @param ownedStatePtr2 Pointer to ownershipState array of cell 2.
   * @param fx2Ptr Pointer to forceX array of cell 2.
   * @param fy2Ptr Pointer to forceY array of cell 2.
   * @param fz2Ptr Pointer to forceZ array of cell 2.
   * @param typeID2Ptr Pointer to typeId array of cell 2.
   */
  template <bool isSingleCell, bool newton3>
  HWY_INLINE void scalarSoAKernel(const size_t soaSize1, const size_t soaSize2, const double *const x1Ptr,
                                  const double *const y1Ptr, const double *const z1Ptr,
                                  const autopas::OwnershipState *const ownedStatePtr1, double *const fx1Ptr,
                                  double *const fy1Ptr, double *const fz1Ptr, const size_t *const typeID1Ptr,
                                  const double *const x2Ptr, const double *const y2Ptr, const double *const z2Ptr,
                                  const autopas::OwnershipState *const ownedStatePtr2, double *const fx2Ptr,
                                  double *const fy2Ptr, double *const fz2Ptr, const size_t *const typeID2Ptr) {
    [[maybe_unused]] double virialX = 0., virialY = 0., virialZ = 0., potentialEnergy = 0.;

    // Local force accumulators kept in registers/stack to eliminate repeated memory read-modify-writes
    std::array<double, 8> fx1Acc{}, fy1Acc{}, fz1Acc{};
    std::array<double, 8> fx2Acc{}, fy2Acc{}, fz2Acc{};

    constexpr size_t startI = isSingleCell ? 1 : 0;
    for (size_t i = startI; i < soaSize1; ++i) {
      if (ownedStatePtr1[i] == autopas::OwnershipState::dummy) {
        continue;
      }
      const double xi = x1Ptr[i], yi = y1Ptr[i], zi = z1Ptr[i];
      double halfSigmaI = 0., sqrtEpsilonI = 0.;
      if constexpr (useMixing) {
        const size_t typeI = typeID1Ptr[i];
        halfSigmaI = _halfSigmaPerType[typeI];
        sqrtEpsilonI = _sqrtEpsilonPerType[typeI];
      }
      [[maybe_unused]] const bool isOwnedI = (ownedStatePtr1[i] == autopas::OwnershipState::owned);

      const size_t endJ = isSingleCell ? i : soaSize2;
      for (size_t j = 0; j < endJ; ++j) {
        if (ownedStatePtr2[j] == autopas::OwnershipState::dummy) {
          continue;
        }
        const double dx = xi - x2Ptr[j];
        const double dy = yi - y2Ptr[j];
        const double dz = zi - z2Ptr[j];
        const double dr2 = dx * dx + dy * dy + dz * dz;

        if (dr2 <= _cutoffSquared) {
          const auto [fx, fy, fz, pot6] =
              computeScalarLJInteraction(dr2, dx, dy, dz, halfSigmaI, sqrtEpsilonI, useMixing ? typeID2Ptr[j] : 0);
          fx1Acc[i] += fx;
          fy1Acc[i] += fy;
          fz1Acc[i] += fz;
          if constexpr (newton3) {
            if constexpr (isSingleCell) {
              fx1Acc[j] -= fx;
              fy1Acc[j] -= fy;
              fz1Acc[j] -= fz;
            } else {
              fx2Acc[j] -= fx;
              fy2Acc[j] -= fy;
              fz2Acc[j] -= fz;
            }
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
    }

    // Flush accumulated forces to memory in a single sequential pass
    for (size_t i = 0; i < soaSize1; ++i) {
      fx1Ptr[i] += fx1Acc[i];
      fy1Ptr[i] += fy1Acc[i];
      fz1Ptr[i] += fz1Acc[i];
    }
    if constexpr (newton3 and not isSingleCell) {
      for (size_t j = 0; j < soaSize2; ++j) {
        fx2Ptr[j] += fx2Acc[j];
        fy2Ptr[j] += fy2Acc[j];
        fz2Ptr[j] += fz2Acc[j];
      }
    }

    if constexpr (calculateGlobals) {
      const int threadID = autopas::autopas_get_thread_num();
      _threadData[threadID].potentialEnergySum += potentialEnergy;
      _threadData[threadID].virialSum[0] += virialX;
      _threadData[threadID].virialSum[1] += virialY;
      _threadData[threadID].virialSum[2] += virialZ;
    }
  }

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
   * Computes the Lennard-Jones pair interaction (force, energy, virial) for a SIMD register.
   *
   * @tparam newton3 Whether Newton's third law is enabled.
   * @param drX Difference vector in x.
   * @param drY Difference vector in y.
   * @param drZ Difference vector in z.
   * @param dr2 Squared distance vector.
   * @param cutoffDummyMask Mask indicating interacting particles within cutoff.
   * @param ownedMaskI Mask indicating owned particle(s) 1.
   * @param ownedMaskJ Mask indicating owned particle(s) 2.
   * @param params Interaction parameters (sigmaSquared, epsilon24, shift6).
   * @param fx Output force component in x direction.
   * @param fy Output force component in y direction.
   * @param fz Output force component in z direction.
   * @param fxAcc Force accumulator in x direction for particle(s) 1.
   * @param fyAcc Force accumulator in y direction for particle(s) 1.
   * @param fzAcc Force accumulator in z direction for particle(s) 1.
   * @param virialSumX Virial accumulator in x direction.
   * @param virialSumY Virial accumulator in y direction.
   * @param virialSumZ Virial accumulator in z direction.
   * @param potentialEnergySum Potential energy accumulator.
   */
  template <bool newton3>
  HWY_INLINE void computeLJInteractionSIMD(
      const VectorDouble &drX, const VectorDouble &drY, const VectorDouble &drZ, const VectorDouble &dr2,
      const MaskDouble &cutoffDummyMask, [[maybe_unused]] const MaskDouble &ownedMaskI,
      [[maybe_unused]] const MaskDouble &ownedMaskJ, const LJPairParameters &params, VectorDouble &fx, VectorDouble &fy,
      VectorDouble &fz, VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc,
      [[maybe_unused]] VectorDouble &virialSumX, [[maybe_unused]] VectorDouble &virialSumY,
      [[maybe_unused]] VectorDouble &virialSumZ, [[maybe_unused]] VectorDouble &potentialEnergySum) const {
    const auto one = highway::Set(tag_double, 1.0);
    const auto invDr2 = highway::Div(one, dr2);
    const auto lj2 = highway::Mul(params.sigmaSquared, invDr2);
    const auto eps24InvDr2 = highway::Mul(params.epsilon24, invDr2);

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
      fx = highway::Mul(drX, facMasked);
      fy = highway::Mul(drY, facMasked);
      fz = highway::Mul(drZ, facMasked);

      fxAcc = highway::Add(fxAcc, fx);
      fyAcc = highway::Add(fyAcc, fy);
      fzAcc = highway::Add(fzAcc, fz);
    }

    if constexpr (calculateGlobals) {
      const auto virialX = highway::Mul(fx, drX);
      const auto virialY = highway::Mul(fy, drY);
      const auto virialZ = highway::Mul(fz, drZ);

      const auto lj12m6 = highway::MulSub(lj6, lj6, lj6);
      const auto potentialEnergy = highway::MulAdd(params.epsilon24, lj12m6, params.shift6);
      const auto potentialEnergyMasked = highway::IfThenElseZero(cutoffDummyMask, potentialEnergy);

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

  /**
   * The kernel calculating the force between two particles (or particle vectors).
   *
   * @tparam newton3 Whether Newton's third law is enabled.
   * @tparam remainderJ Whether inner loop j is in remainder iteration.
   * @tparam vecPattern Vectorization pattern used.
   * @tparam ParameterPolicy Type of parameter policy providing pair parameters.
   * @param j Particle index in inner j loop.
   * @param notDummyMaskI Mask indicating valid particle(s) i.
   * @param ownedMaskI Mask indicating owned particle(s) i.
   * @param ownedStatePtr2 Pointer to ownership states in SoA 2.
   * @param x1 Position x of particle(s) i.
   * @param y1 Position y of particle(s) i.
   * @param z1 Position z of particle(s) i.
   * @param x2Ptr Pointer to x coordinates in SoA 2.
   * @param y2Ptr Pointer to y coordinates in SoA 2.
   * @param z2Ptr Pointer to z coordinates in SoA 2.
   * @param fx2Ptr Pointer to force x array in SoA 2.
   * @param fy2Ptr Pointer to force y array in SoA 2.
   * @param fz2Ptr Pointer to force z array in SoA 2.
   * @param paramPolicy Parameter policy providing Lennard-Jones parameters.
   * @param fxAcc Force accumulator in x direction for particle(s) i.
   * @param fyAcc Force accumulator in y direction for particle(s) i.
   * @param fzAcc Force accumulator in z direction for particle(s) i.
   * @param virialSumX Virial accumulator in x direction.
   * @param virialSumY Virial accumulator in y direction.
   * @param virialSumZ Virial accumulator in z direction.
   * @param potentialEnergySum Potential energy accumulator.
   * @param restJ Number of remaining elements in remainder iteration.
   */
  template <bool newton3, bool remainderJ, VectorizationPattern vecPattern, typename ParameterPolicy>
  HWY_INLINE void SoAKernel(const size_t j, const MaskDouble &notDummyMaskI,
                            [[maybe_unused]] const MaskDouble &ownedMaskI,
                            const int64_t *const __restrict ownedStatePtr2, const VectorDouble &x1,
                            const VectorDouble &y1, const VectorDouble &z1, const double *const __restrict x2Ptr,
                            const double *const __restrict y2Ptr, const double *const __restrict z2Ptr,
                            double *const __restrict fx2Ptr, double *const __restrict fy2Ptr,
                            double *const __restrict fz2Ptr, const ParameterPolicy &paramPolicy, VectorDouble &fxAcc,
                            VectorDouble &fyAcc, VectorDouble &fzAcc, VectorDouble &virialSumX,
                            VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &potentialEnergySum,
                            const size_t restJ) const {
    const auto x2 = loadJRegister<remainderJ, vecPattern>(tag_double, x2Ptr, j, restJ);
    const auto y2 = loadJRegister<remainderJ, vecPattern>(tag_double, y2Ptr, j, restJ);
    const auto z2 = loadJRegister<remainderJ, vecPattern>(tag_double, z2Ptr, j, restJ);

    const auto drX = highway::Sub(x1, x2);
    const auto drY = highway::Sub(y1, y2);
    const auto drZ = highway::Sub(z1, z2);
    const auto dr2 = highway::MulAdd(drZ, drZ, highway::MulAdd(drY, drY, highway::Mul(drX, drX)));

    const auto cutoffMaskI = highway::MaskedLe(notDummyMaskI, dr2, _cutoffSquaredVec);
    if (highway::AllFalse(tag_double, cutoffMaskI)) {
      return;
    }

    const auto ownedStateJLong = loadJRegister<remainderJ, vecPattern>(tag_long, ownedStatePtr2, j, restJ);
    const auto notDummyMaskJ = highway::RebindMask(tag_double, highway::Ne(ownedStateJLong, highway::Zero(tag_long)));
    const auto cutoffDummyMask = highway::And(cutoffMaskI, notDummyMaskJ);

    if (highway::AllFalse(tag_double, cutoffDummyMask)) [[unlikely]] {
      return;
    }

    [[maybe_unused]] MaskDouble ownedMaskJ;
    if constexpr (calculateGlobals and newton3) {
      const VectorLong ownedStateVal = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
      ownedMaskJ = highway::RebindMask(tag_double, highway::Eq(ownedStateJLong, ownedStateVal));
    }

    const auto &params = paramPolicy.template getJ<remainderJ, vecPattern>(j, restJ);

    VectorDouble fx;
    VectorDouble fy;
    VectorDouble fz;
    computeLJInteractionSIMD<newton3>(drX, drY, drZ, dr2, cutoffDummyMask, ownedMaskI, ownedMaskJ, params, fx, fy, fz,
                                      fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, potentialEnergySum);

    if constexpr (newton3) {
      handleNewton3Reduction<remainderJ, vecPattern>(fx, fy, fz, fx2Ptr, fy2Ptr, fz2Ptr, j, restJ);
    }
  }

  /**
   * Kernel for the Verlet lists.
   *
   * @tparam newton3 Whether Newton's third law is enabled.
   * @tparam remainder Whether this is a remainder iteration.
   * @tparam ParameterPolicy Type of parameter policy providing pair parameters.
   * @param j Loop index into the neighbor list.
   * @param ownedMaskI Mask indicating owned particle i.
   * @param x1 Position x of particle i.
   * @param y1 Position y of particle i.
   * @param z1 Position z of particle i.
   * @param xPtr Pointer to x coordinates in SoA.
   * @param yPtr Pointer to y coordinates in SoA.
   * @param zPtr Pointer to z coordinates in SoA.
   * @param ownedStatePtr Pointer to ownership states in SoA.
   * @param fxPtr Pointer to force x array in SoA.
   * @param fyPtr Pointer to force y array in SoA.
   * @param fzPtr Pointer to force z array in SoA.
   * @param paramPolicy Parameter policy providing Lennard-Jones parameters.
   * @param neighborList Pointer to neighbor particle indices.
   * @param fxAcc Force accumulator in x direction for particle i.
   * @param fyAcc Force accumulator in y direction for particle i.
   * @param fzAcc Force accumulator in z direction for particle i.
   * @param virialSumX Virial accumulator in x direction.
   * @param virialSumY Virial accumulator in y direction.
   * @param virialSumZ Virial accumulator in z direction.
   * @param potentialEnergySum Potential energy accumulator.
   * @param rest Number of remaining elements in remainder iteration.
   */
  template <bool newton3, bool remainder = false, typename ParameterPolicy>
  HWY_INLINE void SoAKernelVerlet(const size_t j, [[maybe_unused]] const MaskDouble &ownedMaskI, const VectorDouble &x1,
                                  const VectorDouble &y1, const VectorDouble &z1, const double *const __restrict xPtr,
                                  const double *const __restrict yPtr, const double *const __restrict zPtr,
                                  const int64_t *const __restrict ownedStatePtr, double *const __restrict fxPtr,
                                  double *const __restrict fyPtr, double *const __restrict fzPtr,
                                  const ParameterPolicy &paramPolicy, const size_t *const __restrict neighborList,
                                  VectorDouble &fxAcc, VectorDouble &fyAcc, VectorDouble &fzAcc,
                                  VectorDouble &virialSumX, VectorDouble &virialSumY, VectorDouble &virialSumZ,
                                  VectorDouble &potentialEnergySum, [[maybe_unused]] const size_t rest = 0) const {
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
    const auto dr2 = highway::MulAdd(drZ, drZ, highway::MulAdd(drY, drY, drX2));

    const auto cutoffDummyMask = highway::MaskedLe(notDummyMaskJ, dr2, _cutoffSquaredVec);

    if (highway::AllFalse(tag_double, cutoffDummyMask)) {
      return;
    }

    [[maybe_unused]] MaskDouble ownedMaskJ;
    if constexpr (calculateGlobals and newton3) {
      const VectorLong ownedStateVal = highway::Set(tag_long, static_cast<int64_t>(autopas::OwnershipState::owned));
      ownedMaskJ = highway::RebindMask(tag_double, highway::Eq(ownedState2, ownedStateVal));
    }

    const auto &params = paramPolicy.template getJ<remainder>(indices, restMaskDouble);

    VectorDouble fx;
    VectorDouble fy;
    VectorDouble fz;
    computeLJInteractionSIMD<newton3>(drX, drY, drZ, dr2, cutoffDummyMask, ownedMaskI, ownedMaskJ, params, fx, fy, fz,
                                      fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ, potentialEnergySum);

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
  }
  template <bool reversed, bool newton3, bool remainderI, VectorizationPattern vecPattern, typename ParameterPolicy>
  inline void handleILoopBody(
      const size_t i, const double *const __restrict xPtr1, const double *const __restrict yPtr1,
      const double *const __restrict zPtr1, const autopas::OwnershipState *const __restrict ownedStatePtr1,
      const double *const __restrict xPtr2, const double *const __restrict yPtr2, const double *const __restrict zPtr2,
      const autopas::OwnershipState *const __restrict ownedStatePtr2, double *const __restrict fxPtr1,
      double *const __restrict fyPtr1, double *const __restrict fzPtr1, double *const __restrict fxPtr2,
      double *const __restrict fyPtr2, double *const __restrict fzPtr2, ParameterPolicy &paramPolicy,
      VectorDouble &virialSumX, VectorDouble &virialSumY, VectorDouble &virialSumZ, VectorDouble &potentialEnergySum,
      const size_t restI, const size_t jVecStart, const size_t jVecEnd) const {
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
    paramPolicy.template prepareI<remainderI, reversed, vecPattern>(i, restI);

    constexpr size_t jStep = jStepSize<vecPattern>();
    const size_t restJ = jVecEnd & (jStep - 1);
    const size_t jVecEndFull = jVecEnd - restJ;
    size_t j = jVecStart;
    for (; j < jVecEndFull; j += jStep) {
      SoAKernel<newton3, false, vecPattern>(j, notDummyMaskI, ownedMaskI,
                                            reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2,
                                            zPtr2, fxPtr2, fyPtr2, fzPtr2, paramPolicy, fxAcc, fyAcc, fzAcc, virialSumX,
                                            virialSumY, virialSumZ, potentialEnergySum, 0);
    }
    if (restJ > 0) {
      SoAKernel<newton3, true, vecPattern>(j, notDummyMaskI, ownedMaskI,
                                           reinterpret_cast<const int64_t *>(ownedStatePtr2), x1, y1, z1, xPtr2, yPtr2,
                                           zPtr2, fxPtr2, fyPtr2, fzPtr2, paramPolicy, fxAcc, fyAcc, fzAcc, virialSumX,
                                           virialSumY, virialSumZ, potentialEnergySum, restJ);
    }

    reduceAccumulatedForce<reversed, remainderI, vecPattern>(i, fxPtr1, fyPtr1, fzPtr1, fxAcc, fyAcc, fzAcc, restI);
  }

  /**
   * Helper function running the single-cell reverse i-loop for a given ParameterPolicy.
   *
   * @tparam ParameterPolicy Type of parameter policy.
   * @param soa SoA view of the cell.
   * @param paramPolicy Parameter policy.
   */
  template <typename ParameterPolicy>
  inline void runSingleLoop(autopas::SoAView<SoAArraysType> soa, ParameterPolicy &paramPolicy) {
    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    auto virialSumX = highway::Zero(tag_double);
    auto virialSumY = highway::Zero(tag_double);
    auto virialSumZ = highway::Zero(tag_double);
    auto potentialEnergySum = highway::Zero(tag_double);

    for (size_t i = soa.size() - 1; i >= 1; i -= 1) {
      static_assert(std::is_same_v<std::underlying_type_t<autopas::OwnershipState>, int64_t>,
                    "OwnershipStates underlying type should be int64_t!");

      handleILoopBody<true, true, false, VectorizationPattern::p1xVec>(
          i, xPtr, yPtr, zPtr, ownedStatePtr, xPtr, yPtr, zPtr, ownedStatePtr, fxPtr, fyPtr, fzPtr, fxPtr, fyPtr, fzPtr,
          paramPolicy, virialSumX, virialSumY, virialSumZ, potentialEnergySum, 0, 0, i);
    }

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, potentialEnergySum);
    }
  }

  /**
   * Helper function running the cell-pair i-loop for a given ParameterPolicy.
   *
   * @tparam newton3 Whether Newton's third law is enabled.
   * @tparam sorted Whether sortingData is provided and used.
   * @tparam vecPattern Vectorization pattern used.
   * @tparam ParameterPolicy Type of parameter policy.
   * @param soa1 First SoA view.
   * @param soa2 Second SoA view.
   * @param sortingData Sorting data for index bounds pruning.
   * @param paramPolicy Parameter policy.
   */
  template <bool newton3, bool sorted, VectorizationPattern vecPattern, typename ParameterPolicy>
  inline void runPairLoop(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                          const std::optional<std::reference_wrapper<const autopas::SoASortingData>> &sortingData,
                          ParameterPolicy &paramPolicy) {
    const size_t n1 = soa1.size();
    const size_t n2 = soa2.size();

    const auto *const __restrict x1Ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y1Ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z1Ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle_T::AttributeNames::ownershipState>();

    const auto *const __restrict x2Ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y2Ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z2Ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();
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
    VectorDouble potentialEnergySum = highway::Zero(tag_double);

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
          fy2Ptr, fz2Ptr, paramPolicy, virialSumX, virialSumY, virialSumZ, potentialEnergySum, 0, jVecStart, jVecEnd);
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
          handleILoopBody<false, newton3, true, vecPattern>(
              i, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx1Ptr, fy1Ptr, fz1Ptr,
              fx2Ptr, fy2Ptr, fz2Ptr, paramPolicy, virialSumX, virialSumY, virialSumZ, potentialEnergySum,
              static_cast<size_t>(restI), jVecStart, jVecEnd);
        }
      }
    }

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, potentialEnergySum);
    }
  }

  /**
   * Helper function running the Verlet list loop for a given ParameterPolicy.
   *
   * @tparam newton3 Whether Newton's third law is enabled.
   * @tparam ParameterPolicy Type of parameter policy.
   * @param soa SoA view of the cell.
   * @param indexFirst Index of base particle.
   * @param neighborList Neighbor indices.
   * @param paramPolicy Parameter policy.
   */
  template <bool newton3, typename ParameterPolicy>
  inline void runVerletLoop(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            std::span<const size_t> neighborList, ParameterPolicy &paramPolicy) {
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();

    auto *const __restrict fxPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

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
    [[maybe_unused]] MaskDouble ownedMaskI;
    if constexpr (calculateGlobals) {
      constexpr auto ownedVal = static_cast<double>(static_cast<int64_t>(autopas::OwnershipState::owned));
      ownedMaskI = highway::Eq(ownedStateI, highway::Set(tag_double, ownedVal));
    }

    size_t j = 0;
    const size_t neighborListSize = neighborList.size();
    const size_t rest = neighborListSize & (_vecLengthDouble - 1);
    const size_t vecEnd = neighborListSize - rest;

    for (; j < vecEnd; j += _vecLengthDouble) {
      SoAKernelVerlet<newton3, false>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr,
                                      reinterpret_cast<const int64_t *>(ownedStatePtr), fxPtr, fyPtr, fzPtr,
                                      paramPolicy, neighborList.data(), fxAcc, fyAcc, fzAcc, virialSumX, virialSumY,
                                      virialSumZ, potentialEnergySum, 0);
    }

    if (rest > 0) {
      SoAKernelVerlet<newton3, true>(j, ownedMaskI, x1, y1, z1, xPtr, yPtr, zPtr,
                                     reinterpret_cast<const int64_t *>(ownedStatePtr), fxPtr, fyPtr, fzPtr, paramPolicy,
                                     neighborList.data(), fxAcc, fyAcc, fzAcc, virialSumX, virialSumY, virialSumZ,
                                     potentialEnergySum, rest);
    }

    reduceAdd3_InterleaveFirst(tag_double, fxAcc, fyAcc, fzAcc, fxPtr[indexFirst], fyPtr[indexFirst],
                               fzPtr[indexFirst]);

    if constexpr (calculateGlobals) {
      computeGlobals(virialSumX, virialSumY, virialSumZ, potentialEnergySum);
    }
  }

  /**
   * Dispatches SoAFunctorPair or SoAFunctorPairSorted to the appropriate VectorizationPattern implementation.
   *
   * @tparam sorted Whether sortingData is provided and utilized.
   * @param soa1 First SoA view.
   * @param soa2 Second SoA view.
   * @param newton3 Whether Newton's third law is enabled.
   * @param sortingData Optional sorting data for index bounds pruning.
   */
  template <bool sorted>
  inline void SoAFunctorPairDispatch(
      autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, const bool newton3,
      [[maybe_unused]] const std::optional<std::reference_wrapper<const autopas::SoASortingData>> &sortingData =
          std::nullopt) {
    const size_t soaSize1 = soa1.size();
    const size_t soaSize2 = soa2.size();

    if (soaSize1 <= 7 and soaSize2 <= 7) {
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
      const auto *const __restrict typeID1Ptr =
          useMixing ? soa1.template begin<Particle_T::AttributeNames::typeId>() : nullptr;

      const auto *const __restrict x2Ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
      const auto *const __restrict y2Ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
      const auto *const __restrict z2Ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();
      const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle_T::AttributeNames::ownershipState>();
      auto *const __restrict fx2Ptr = soa2.template begin<Particle_T::AttributeNames::forceX>();
      auto *const __restrict fy2Ptr = soa2.template begin<Particle_T::AttributeNames::forceY>();
      auto *const __restrict fz2Ptr = soa2.template begin<Particle_T::AttributeNames::forceZ>();
      const auto *const __restrict typeID2Ptr =
          useMixing ? soa2.template begin<Particle_T::AttributeNames::typeId>() : nullptr;

      if (newton3) {
        scalarSoAKernel<false, true>(soaSize1, soaSize2, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, fx1Ptr, fy1Ptr, fz1Ptr,
                                     typeID1Ptr, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr, fy2Ptr, fz2Ptr,
                                     typeID2Ptr);
      } else {
        scalarSoAKernel<false, false>(soaSize1, soaSize2, x1Ptr, y1Ptr, z1Ptr, ownedStatePtr1, fx1Ptr, fy1Ptr, fz1Ptr,
                                      typeID1Ptr, x2Ptr, y2Ptr, z2Ptr, ownedStatePtr2, fx2Ptr, fy2Ptr, fz2Ptr,
                                      typeID2Ptr);
      }
      return;
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
   * @tparam newton3 Whether Newton's third law is enabled.
   * @tparam sorted Whether to use pruned index data in sortingData.
   * @tparam vecPattern Vectorization pattern.
   * @param soa1 First SoA view.
   * @param soa2 Second SoA view.
   * @param sortingData Pruning data needed for sorting optimization.
   */
  template <bool newton3, bool sorted, VectorizationPattern vecPattern>
  inline void SoAFunctorPairImpl(
      autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
      const std::optional<std::reference_wrapper<const autopas::SoASortingData>> &sortingData = std::nullopt) {
    const size_t n1 = soa1.size();
    const size_t n2 = soa2.size();
    if (n1 == 0 or n2 == 0) {
      return;
    }

    if constexpr (useMixing) {
      const auto *const __restrict typeID1Ptr = soa1.template begin<Particle_T::AttributeNames::typeId>();
      const auto *const __restrict typeID2Ptr = soa2.template begin<Particle_T::AttributeNames::typeId>();
      if (soaIsHomogeneous(typeID1Ptr, n1) and soaIsHomogeneous(typeID2Ptr, n2)) {
        const size_t type1 = n1 > 0 ? typeID1Ptr[0] : 0;
        const size_t type2 = n2 > 0 ? typeID2Ptr[0] : 0;
        auto paramPolicy = makeHomogeneousParameterPolicy(type1, type2);
        runPairLoop<newton3, sorted, vecPattern>(soa1, soa2, sortingData, paramPolicy);
      } else {
        MixingBuffers mixingBuffers;
        double *halfSigma1Ptr = nullptr;
        double *sqrtEpsilon1Ptr = nullptr;
        double *halfSigma2Ptr = nullptr;
        double *sqrtEpsilon2Ptr = nullptr;
        mixingBuffers.template getMixingPointers<0>(typeID1Ptr, n1, _halfSigmaPerType, _sqrtEpsilonPerType,
                                                    halfSigma1Ptr, sqrtEpsilon1Ptr);
        mixingBuffers.template getMixingPointers<1>(typeID2Ptr, n2, _halfSigmaPerType, _sqrtEpsilonPerType,
                                                    halfSigma2Ptr, sqrtEpsilon2Ptr);
        HeterogeneousParameterPolicy paramPolicy{halfSigma1Ptr, sqrtEpsilon1Ptr, halfSigma2Ptr, sqrtEpsilon2Ptr,
                                                 _invCutoffSquared};
        runPairLoop<newton3, sorted, vecPattern>(soa1, soa2, sortingData, paramPolicy);
      }
    } else {
      auto paramPolicy = makeHomogeneousParameterPolicy();
      runPairLoop<newton3, sorted, vecPattern>(soa1, soa2, sortingData, paramPolicy);
    }
  }

  /**
   * Templatized implementation of SoAFunctorVerlet.
   *
   * @tparam newton3 Whether Newton's third law is enabled.
   * @param soa SoA view of the cell.
   * @param indexFirst Index of base particle.
   * @param neighborList Neighbor indices list.
   */
  template <bool newton3>
  inline void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                                   std::span<const size_t> neighborList) {
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();
    if (ownedStatePtr[indexFirst] == autopas::OwnershipState::dummy) {
      return;
    }

    const bool isHomogeneous = not useMixing or _halfSigmaPerType.size() <= 1;

    if constexpr (useMixing) {
      if (not isHomogeneous) {
        const auto *const typeIDPtr = soa.template begin<Particle_T::AttributeNames::typeId>();
        const size_t typeI = typeIDPtr[indexFirst];
        HeterogeneousParameterPolicy paramPolicy{
            _halfSigmaPerType[typeI], _sqrtEpsilonPerType[typeI], reinterpret_cast<const int64_t *>(typeIDPtr),
            _halfSigmaPerType.data(), _sqrtEpsilonPerType.data(), _invCutoffSquared};
        runVerletLoop<newton3>(soa, indexFirst, neighborList, paramPolicy);
      } else {
        auto paramPolicy = makeHomogeneousParameterPolicy();
        runVerletLoop<newton3>(soa, indexFirst, neighborList, paramPolicy);
      }
    } else {
      auto paramPolicy = makeHomogeneousParameterPolicy();
      runVerletLoop<newton3>(soa, indexFirst, neighborList, paramPolicy);
    }
  }
  // cutoff squared used in the functor.
  const double _cutoffSquared{0.};
  // inverse cutoff squared used in the functor.
  const double _invCutoffSquared{0.};
  // cutoff squared pre-broadcast into a vector.
  const VectorDouble _cutoffSquaredVec{highway::Set(tag_double, 0.)};
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
