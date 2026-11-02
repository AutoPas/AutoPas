/**
 * @file LJFunctor.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "ParticlePropertiesLibrary.h"
#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace mdLib {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 *
 * This functor is applicable with the 1xVectorLength pattern only.
 *
 * @tparam Particle_T The type of particle.
 * @tparam applyShift Switch for the lj potential to be truncated shifted.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle_T, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctor
    : public autopas::PairwiseFunctor<Particle_T, LJFunctor<Particle_T, applyShift, useMixing, useNewton3,
                                                            calculateGlobals, countFLOPs, relevantForTuning>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Precision of SoA entries.
   */
  using SoAFloatPrecision = typename Particle_T::ParticleSoAFloatPrecision;

  template <size_t StackCapacity = 256>
  struct MixingBuffers {
    alignas(64) std::array<SoAFloatPrecision, StackCapacity> halfSigmaStack;
    alignas(64) std::array<SoAFloatPrecision, StackCapacity> sqrtEpsilonStack;

    struct Pointers {
      const SoAFloatPrecision *__restrict halfSigmaPtr;
      const SoAFloatPrecision *__restrict sqrtEpsilonPtr;
    };

    /**
     * Projects halfSigma and sqrtEpsilon values for the given typeIDs into a 64-byte aligned buffer.
     * If the number of particles is less than or equal to StackCapacity, the values are stored in the stack buffer.
     * Otherwise, they are stored in a thread-local heap buffer.
     *
     * @tparam BufferIndex An optional template parameter to differentiate between multiple MixingBuffers instances.
     * @param typeIDs Pointer to the array of type IDs for the particles.
     * @param soaSize The number of particles in the SoA.
     * @param halfSigmaPerType Vector containing halfSigma values for each particle type.
     * @param sqrtEpsilonPerType Vector containing sqrtEpsilon values for each particle type.
     * @return Pointers containing the halfSigma and sqrtEpsilon __restrict pointers.
     */
    template <size_t BufferIndex = 0>
    Pointers getMixingPointers(const size_t *__restrict typeIDs, const size_t soaSize,
                               const std::vector<SoAFloatPrecision> &halfSigmaPerType,
                               const std::vector<SoAFloatPrecision> &sqrtEpsilonPerType) {
      SoAFloatPrecision *rawHalfSigmaPtr = nullptr;
      SoAFloatPrecision *rawSqrtEpsilonPtr = nullptr;
      if (soaSize <= StackCapacity) [[likely]] {
        rawHalfSigmaPtr = halfSigmaStack.data();
        rawSqrtEpsilonPtr = sqrtEpsilonStack.data();
      } else [[unlikely]] {
        thread_local std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision, 64>> halfSigmaHeap;
        thread_local std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision, 64>> sqrtEpsilonHeap;
        if (halfSigmaHeap.size() < soaSize) {
          halfSigmaHeap.resize(soaSize);
          sqrtEpsilonHeap.resize(soaSize);
        }
        rawHalfSigmaPtr = halfSigmaHeap.data();
        rawSqrtEpsilonPtr = sqrtEpsilonHeap.data();
      }
#pragma omp simd
      for (size_t k = 0; k < soaSize; ++k) {
        const auto typeID = typeIDs[k];
        rawHalfSigmaPtr[k] = halfSigmaPerType[typeID];
        rawSqrtEpsilonPtr[k] = sqrtEpsilonPerType[typeID];
      }
      return {rawHalfSigmaPtr, rawSqrtEpsilonPtr};
    }
  };

  struct HomogeneousParams {
    SoAFloatPrecision sigmaSquared;
    SoAFloatPrecision epsilon24;
    SoAFloatPrecision shift6;
  };

  struct SoAPointers {
    const SoAFloatPrecision *const __restrict x;
    const SoAFloatPrecision *const __restrict y;
    const SoAFloatPrecision *const __restrict z;
    const autopas::OwnershipState *const __restrict owned;
    SoAFloatPrecision *const __restrict fx;
    SoAFloatPrecision *const __restrict fy;
    SoAFloatPrecision *const __restrict fz;
    const size_t *const __restrict type;
    const SoAFloatPrecision *__restrict halfSigma;
    const SoAFloatPrecision *__restrict sqrtEpsilon;

    static SoAPointers extract(autopas::SoAView<SoAArraysType> soa,
                               const SoAFloatPrecision *const __restrict halfSigma = nullptr,
                               const SoAFloatPrecision *const __restrict sqrtEpsilon = nullptr) {
      return {soa.template begin<Particle_T::AttributeNames::posX>(),
              soa.template begin<Particle_T::AttributeNames::posY>(),
              soa.template begin<Particle_T::AttributeNames::posZ>(),
              soa.template begin<Particle_T::AttributeNames::ownershipState>(),
              soa.template begin<Particle_T::AttributeNames::forceX>(),
              soa.template begin<Particle_T::AttributeNames::forceY>(),
              soa.template begin<Particle_T::AttributeNames::forceZ>(),
              soa.template begin<Particle_T::AttributeNames::typeId>(),
              halfSigma,
              sqrtEpsilon};
    }
  };

  struct GlobalsAccumulator {
    SoAFloatPrecision potentialEnergySum{0.};
    SoAFloatPrecision virialX{0.};
    SoAFloatPrecision virialY{0.};
    SoAFloatPrecision virialZ{0.};
    size_t numDistCalls{0};
    size_t numKernelCallsN3{0};
    size_t numKernelCallsNoN3{0};
    size_t numGlobalCalcsN3{0};
    size_t numGlobalCalcsNoN3{0};
  };

  /**
   * Evaluates the Lennard-Jones interaction between two particles given their displacements and parameters.
   *
   * @param drx x displacement
   * @param dry y displacement
   * @param drz z displacement
   * @param dr2 distance squared
   * @param mask SIMD mask indicating if the interaction is valid
   * @param sigmaSquared sigma^2
   * @param epsilon24 24 * epsilon
   * @param fx output force x
   * @param fy output force y
   * @param fz output force z
   * @param lj12m6 output (lj12 - lj6) for potential energy calculation
   */
  [[gnu::always_inline]] static inline void computeLJForceCore(const SoAFloatPrecision drx, const SoAFloatPrecision dry,
                                                               const SoAFloatPrecision drz, const SoAFloatPrecision dr2,
                                                               const bool mask, const SoAFloatPrecision sigmaSquared,
                                                               const SoAFloatPrecision epsilon24, SoAFloatPrecision &fx,
                                                               SoAFloatPrecision &fy, SoAFloatPrecision &fz,
                                                               SoAFloatPrecision &lj12m6) {
    const SoAFloatPrecision invdr2 = 1. / dr2;
    const SoAFloatPrecision lj2 = sigmaSquared * invdr2;
    const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
    const SoAFloatPrecision lj12 = lj6 * lj6;
    lj12m6 = lj12 - lj6;
    const SoAFloatPrecision fac = mask ? epsilon24 * (lj12 + lj12m6) * invdr2 : 0.;
    fx = drx * fac;
    fy = dry * fac;
    fz = drz * fac;
  }

  /**
   * Retrieves homogeneous parameters for a cell starting with typeptr[0].
   */
  HomogeneousParams getHomogeneousParams(const size_t *const __restrict typeptr) const {
    SoAFloatPrecision const_sigmaSquared = _sigmaSquared;
    SoAFloatPrecision const_epsilon24 = _epsilon24;
    SoAFloatPrecision const_shift6 = _shift6;

    if constexpr (useMixing) {
      const size_t firstType = typeptr[0];
      const auto halfSig = _halfSigmaPerType[firstType];
      const auto sqrtEps = _sqrtEpsilonPerType[firstType];
      const_sigmaSquared = static_cast<SoAFloatPrecision>(4.0) * halfSig * halfSig;
      const_epsilon24 = static_cast<SoAFloatPrecision>(24.0) * sqrtEps * sqrtEps;
      if constexpr (calculateGlobals and applyShift) {
        const_shift6 = computeShift6(const_epsilon24, const_sigmaSquared, _invCutoffSquared);
      }
    }
    return {const_sigmaSquared, const_epsilon24, const_shift6};
  }

  /**
   * Flushes local thread accumulators to global thread-local arrays.
   */
  void flushThreadData(const GlobalsAccumulator &acc) {
    if constexpr (countFLOPs or calculateGlobals) {
      const auto threadID = autopas::autopas_get_thread_num();
      if constexpr (countFLOPs) {
        _aosThreadDataFLOPs[threadID].numDistCalls += acc.numDistCalls;
        _aosThreadDataFLOPs[threadID].numKernelCallsNoN3 += acc.numKernelCallsNoN3;
        _aosThreadDataFLOPs[threadID].numKernelCallsN3 += acc.numKernelCallsN3;
        _aosThreadDataFLOPs[threadID].numGlobalCalcsNoN3 += acc.numGlobalCalcsNoN3;
        _aosThreadDataFLOPs[threadID].numGlobalCalcsN3 += acc.numGlobalCalcsN3;
      }
      if constexpr (calculateGlobals) {
        _aosThreadDataGlobals[threadID].potentialEnergySum += acc.potentialEnergySum;
        _aosThreadDataGlobals[threadID].virialSum[0] += acc.virialX;
        _aosThreadDataGlobals[threadID].virialSum[1] += acc.virialY;
        _aosThreadDataGlobals[threadID].virialSum[2] += acc.virialZ;
      }
    }
  }

 public:
  /**
   * Deleted default constructor
   */
  LJFunctor() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit LJFunctor(double cutoff, void * /*dummy*/)
      : autopas::PairwiseFunctor<Particle_T, LJFunctor<Particle_T, applyShift, useMixing, useNewton3, calculateGlobals,
                                                       countFLOPs, relevantForTuning>>(cutoff),
        _cutoffSquared{cutoff * cutoff},
        _invCutoffSquared{1. / (cutoff * cutoff)},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _postProcessed{false} {
    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
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
  explicit LJFunctor(double cutoff) : LJFunctor(cutoff, nullptr) {
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
  explicit LJFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : LJFunctor(cutoff, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    initFromPPL(particlePropertiesLibrary);
  }

  std::string getName() final { return "LJFunctorAutoVec"; }

  bool isRelevantForTuning() final { return relevantForTuning; }

  bool allowsNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Only or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  bool allowsNonNewton3() final {
    return useNewton3 == autopas::FunctorN3Modes::Newton3Off or useNewton3 == autopas::FunctorN3Modes::Both;
  }

  void AoSFunctor(Particle_T &i, Particle_T &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy()) {
      return;
    }

    if constexpr (countFLOPs) {
      const auto threadID = autopas::autopas_get_thread_num();
      ++_aosThreadDataFLOPs[threadID].numDistCalls;
    }

    auto dr = i.getR() - j.getR();
    double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffSquared) {
      return;
    }

    auto sigmaSquared = _sigmaSquared;
    auto epsilon24 = _epsilon24;
    [[maybe_unused]] auto shift6 = _shift6;
    if constexpr (useMixing) {
      const auto typeI = i.getTypeId();
      const auto typeJ = j.getTypeId();
      const auto sigma = _halfSigmaPerType[typeI] + _halfSigmaPerType[typeJ];
      sigmaSquared = sigma * sigma;
      epsilon24 = 24.0 * _sqrtEpsilonPerType[typeI] * _sqrtEpsilonPerType[typeJ];
      if constexpr (calculateGlobals and applyShift) {
        shift6 = computeShift6(epsilon24, sigmaSquared, _invCutoffSquared);
      }
    }

    double invdr2 = 1. / dr2;
    double lj6 = sigmaSquared * invdr2;
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

    if constexpr (calculateGlobals or countFLOPs) {
      const auto threadID = autopas::autopas_get_thread_num();
      if constexpr (countFLOPs) {
        if (newton3) {
          ++_aosThreadDataFLOPs[threadID].numKernelCallsN3;
        } else {
          ++_aosThreadDataFLOPs[threadID].numKernelCallsNoN3;
        }
      }

      if constexpr (calculateGlobals) {
        // We always add the full contribution for each owned particle and divide the sums by 2 in endTraversal().
        // Potential energy has an additional factor of 6, which is also handled in endTraversal().

        auto virial = dr * f;
        double potentialEnergy6 = epsilon24 * lj12m6 + shift6;

        if (i.isOwned()) {
          _aosThreadDataGlobals[threadID].potentialEnergySum += potentialEnergy6;
          _aosThreadDataGlobals[threadID].virialSum += virial;
        }
        // for non-newton3 the second particle will be considered in a separate calculation
        if (newton3 and j.isOwned()) {
          _aosThreadDataGlobals[threadID].potentialEnergySum += potentialEnergy6;
          _aosThreadDataGlobals[threadID].virialSum += virial;
        }
        if constexpr (countFLOPs) {
          if (newton3) {
            ++_aosThreadDataFLOPs[threadID].numGlobalCalcsN3;
          } else {
            ++_aosThreadDataFLOPs[threadID].numGlobalCalcsNoN3;
          }
        }
      }
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversal of the soa.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    const size_t n = soa.size();
    if (n <= 1) return;

    auto ptrs = SoAPointers::extract(soa);
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    [[maybe_unused]] const SoAFloatPrecision invCutoffSquared = _invCutoffSquared;
    GlobalsAccumulator acc;

    // Small-N Direct Scalar Shortcut (N <= 8)
    // Avoids SIMD peel/tail loop and vector mask/reduction overhead on small cell particle counts
    if (n <= 8) {
      for (size_t i = 0; i < n - 1; ++i) {
        const auto ownedStateI = ptrs.owned[i];
        if (ownedStateI == autopas::OwnershipState::dummy) continue;

        const SoAFloatPrecision xi = ptrs.x[i], yi = ptrs.y[i], zi = ptrs.z[i];
        [[maybe_unused]] const SoAFloatPrecision ownedFactorI =
            (ownedStateI == autopas::OwnershipState::owned) ? 1. : 0.;

        [[maybe_unused]] SoAFloatPrecision halfSigI = 0., sqrtEpsI = 0.;
        if constexpr (useMixing) {
          const size_t ti = ptrs.type[i];
          halfSigI = _halfSigmaPerType[ti];
          sqrtEpsI = _sqrtEpsilonPerType[ti];
        }

        SoAFloatPrecision fxi = 0., fyi = 0., fzi = 0.;

        for (size_t j = i + 1; j < n; ++j) {
          const auto ownedStateJ = ptrs.owned[j];
          if (ownedStateJ == autopas::OwnershipState::dummy) continue;

          const auto dx = xi - ptrs.x[j], dy = yi - ptrs.y[j], dz = zi - ptrs.z[j];
          const auto dr2 = dx * dx + dy * dy + dz * dz;

          if (dr2 <= cutoffSquared) {
            SoAFloatPrecision sigSq = _sigmaSquared;
            SoAFloatPrecision eps24 = _epsilon24;
            [[maybe_unused]] SoAFloatPrecision shift6 = _shift6;

            if constexpr (useMixing) {
              const size_t tj = ptrs.type[j];
              const auto sig = halfSigI + _halfSigmaPerType[tj];
              sigSq = sig * sig;
              eps24 = static_cast<SoAFloatPrecision>(24.0) * sqrtEpsI * _sqrtEpsilonPerType[tj];
              if constexpr (calculateGlobals and applyShift) {
                shift6 = computeShift6(eps24, sigSq, invCutoffSquared);
              }
            }

            SoAFloatPrecision dfx, dfy, dfz, lj12m6;
            computeLJForceCore(dx, dy, dz, dr2, true, sigSq, eps24, dfx, dfy, dfz, lj12m6);

            fxi += dfx;
            fyi += dfy;
            fzi += dfz;
            ptrs.fx[j] -= dfx;
            ptrs.fy[j] -= dfy;
            ptrs.fz[j] -= dfz;

            if constexpr (countFLOPs) {
              ++acc.numDistCalls;
              ++acc.numKernelCallsN3;
            }

            if constexpr (calculateGlobals) {
              const auto pot6 = eps24 * lj12m6 + shift6;
              const SoAFloatPrecision energyFactor =
                  ownedFactorI + (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.);
              acc.potentialEnergySum += pot6 * energyFactor;
              acc.virialX += dfx * dx * energyFactor;
              acc.virialY += dfy * dy * energyFactor;
              acc.virialZ += dfz * dz * energyFactor;

              if constexpr (countFLOPs) {
                ++acc.numGlobalCalcsN3;
              }
            }
          } else if constexpr (countFLOPs) {
            ++acc.numDistCalls;
          }
        }
        ptrs.fx[i] += fxi;
        ptrs.fy[i] += fyi;
        ptrs.fz[i] += fzi;
      }

      flushThreadData(acc);
      return;
    }

    const bool isHomogeneous = soaIsHomogeneous(ptrs.type, n);

    if (not useMixing or isHomogeneous) {
      const auto params = getHomogeneousParams(ptrs.type);
      computeCellInteractions<true, true, true>(n, n, ptrs, ptrs, params, acc);
    } else if constexpr (useMixing) {
      MixingBuffers<512> mixingBuffers;
      const auto mixingPtrs = mixingBuffers.getMixingPointers(ptrs.type, n, _halfSigmaPerType, _sqrtEpsilonPerType);
      ptrs.halfSigma = mixingPtrs.halfSigmaPtr;
      ptrs.sqrtEpsilon = mixingPtrs.sqrtEpsilonPtr;
      computeCellInteractions<true, false, true>(n, n, ptrs, ptrs, {}, acc);
    }

    flushThreadData(acc);
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2,
                      const bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
    }
  }

 private:
  template <bool isSingleCell, bool isHomogeneous, bool newton3>
  void computeCellInteractions(const size_t n1, const size_t n2, const SoAPointers &cell1, const SoAPointers &cell2,
                               const HomogeneousParams &params, GlobalsAccumulator &acc) {
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    [[maybe_unused]] const SoAFloatPrecision invCutoffSquared = _invCutoffSquared;

    const auto *const __restrict x1ptr = cell1.x;
    const auto *const __restrict y1ptr = cell1.y;
    const auto *const __restrict z1ptr = cell1.z;
    const auto *const __restrict x2ptr = cell2.x;
    const auto *const __restrict y2ptr = cell2.y;
    const auto *const __restrict z2ptr = cell2.z;
    const auto *const __restrict ownedStatePtr1 = cell1.owned;
    const auto *const __restrict ownedStatePtr2 = cell2.owned;
    auto *const __restrict fx1ptr = cell1.fx;
    auto *const __restrict fy1ptr = cell1.fy;
    auto *const __restrict fz1ptr = cell1.fz;
    auto *const __restrict fx2ptr = cell2.fx;
    auto *const __restrict fy2ptr = cell2.fy;
    auto *const __restrict fz2ptr = cell2.fz;

    const auto *const __restrict halfSigma1Ptr = cell1.halfSigma;
    const auto *const __restrict sqrtEpsilon1Ptr = cell1.sqrtEpsilon;
    const auto *const __restrict halfSigma2Ptr = cell2.halfSigma;
    const auto *const __restrict sqrtEpsilon2Ptr = cell2.sqrtEpsilon;

    const SoAFloatPrecision const_sigmaSquared = params.sigmaSquared;
    const SoAFloatPrecision const_epsilon24 = params.epsilon24;
    const SoAFloatPrecision const_shift6 = params.shift6;

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    const size_t outerLimit = isSingleCell ? n1 - 1 : n1;

    for (unsigned int i = 0; i < outerLimit; ++i) {
      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision xi = x1ptr[i];
      const SoAFloatPrecision yi = y1ptr[i];
      const SoAFloatPrecision zi = z1ptr[i];
      [[maybe_unused]] const SoAFloatPrecision ownedFactorI = (ownedStateI == autopas::OwnershipState::owned) ? 1. : 0.;

      [[maybe_unused]] SoAFloatPrecision halfSigmaI = 0.;
      [[maybe_unused]] SoAFloatPrecision sqrtEpsilon24I = 0.;
      if constexpr (not isHomogeneous) {
        halfSigmaI = halfSigma1Ptr[i];
        sqrtEpsilon24I = 24.0 * sqrtEpsilon1Ptr[i];
      }

      SoAFloatPrecision fxacc = 0.;
      SoAFloatPrecision fyacc = 0.;
      SoAFloatPrecision fzacc = 0.;

      const size_t jStart = isSingleCell ? i + 1 : 0;

      if constexpr (not calculateGlobals and not countFLOPs) {
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
        for (unsigned int j = jStart; j < n2; ++j) {
          SoAFloatPrecision sigmaSquared = const_sigmaSquared;
          SoAFloatPrecision epsilon24 = const_epsilon24;
          if constexpr (not isHomogeneous) {
            const SoAFloatPrecision sig = halfSigmaI + halfSigma2Ptr[j];
            sigmaSquared = sig * sig;
            epsilon24 = sqrtEpsilon24I * sqrtEpsilon2Ptr[j];
          }

          const auto ownedStateJ = ownedStatePtr2[j];

          const SoAFloatPrecision drx = xi - x2ptr[j];
          const SoAFloatPrecision dry = yi - y2ptr[j];
          const SoAFloatPrecision drz = zi - z2ptr[j];

          const SoAFloatPrecision dr2 = drx * drx + dry * dry + drz * drz;

          const bool mask = (dr2 <= cutoffSquared) & (ownedStateJ != autopas::OwnershipState::dummy);

          SoAFloatPrecision fx, fy, fz, lj12m6;
          computeLJForceCore(drx, dry, drz, dr2, mask, sigmaSquared, epsilon24, fx, fy, fz, lj12m6);

          fxacc += fx;
          fyacc += fy;
          fzacc += fz;
          if constexpr (newton3) {
            fx2ptr[j] -= fx;
            fy2ptr[j] -= fy;
            fz2ptr[j] -= fz;
          }
        }
      } else {
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ, \
                               numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum,       \
                               numGlobalCalcsN3Sum, numGlobalCalcsNoN3Sum)
        for (unsigned int j = jStart; j < n2; ++j) {
          SoAFloatPrecision sigmaSquared = const_sigmaSquared;
          SoAFloatPrecision epsilon24 = const_epsilon24;
          SoAFloatPrecision shift6 = const_shift6;
          if constexpr (not isHomogeneous) {
            const SoAFloatPrecision sig = halfSigmaI + halfSigma2Ptr[j];
            sigmaSquared = sig * sig;
            epsilon24 = sqrtEpsilon24I * sqrtEpsilon2Ptr[j];
            if constexpr (calculateGlobals and applyShift) {
              shift6 = computeShift6(epsilon24, sigmaSquared, invCutoffSquared);
            }
          }

          const auto ownedStateJ = ownedStatePtr2[j];

          const SoAFloatPrecision drx = xi - x2ptr[j];
          const SoAFloatPrecision dry = yi - y2ptr[j];
          const SoAFloatPrecision drz = zi - z2ptr[j];

          const SoAFloatPrecision dr2 = drx * drx + dry * dry + drz * drz;

          const bool mask = (dr2 <= cutoffSquared) & (ownedStateJ != autopas::OwnershipState::dummy);

          SoAFloatPrecision fx, fy, fz, lj12m6;
          computeLJForceCore(drx, dry, drz, dr2, mask, sigmaSquared, epsilon24, fx, fy, fz, lj12m6);

          fxacc += fx;
          fyacc += fy;
          fzacc += fz;
          if constexpr (newton3) {
            fx2ptr[j] -= fx;
            fy2ptr[j] -= fy;
            fz2ptr[j] -= fz;
          }

          if constexpr (countFLOPs) {
            numDistanceCalculationSum += (ownedStateJ != autopas::OwnershipState::dummy) ? 1 : 0;
            if constexpr (newton3) {
              numKernelCallsN3Sum += mask ? 1 : 0;
            } else {
              numKernelCallsNoN3Sum += mask ? 1 : 0;
            }
          }

          if constexpr (calculateGlobals) {
            const SoAFloatPrecision virialx = drx * fx;
            const SoAFloatPrecision virialy = dry * fy;
            const SoAFloatPrecision virialz = drz * fz;
            const SoAFloatPrecision potentialEnergy6 = mask ? (epsilon24 * lj12m6 + shift6) : 0.;

            SoAFloatPrecision energyFactor = ownedFactorI;
            if constexpr (newton3) {
              energyFactor += (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.);
            }

            potentialEnergySum += potentialEnergy6 * energyFactor;
            virialSumX += virialx * energyFactor;
            virialSumY += virialy * energyFactor;
            virialSumZ += virialz * energyFactor;

            if constexpr (countFLOPs) {
              if constexpr (newton3) {
                numGlobalCalcsN3Sum += mask ? 1 : 0;
              } else {
                numGlobalCalcsNoN3Sum += mask ? 1 : 0;
              }
            }
          }
        }
      }
      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }

    if constexpr (calculateGlobals) {
      acc.potentialEnergySum += potentialEnergySum;
      acc.virialX += virialSumX;
      acc.virialY += virialSumY;
      acc.virialZ += virialSumZ;
    }
    if constexpr (countFLOPs) {
      acc.numDistCalls += numDistanceCalculationSum;
      acc.numKernelCallsN3 += numKernelCallsN3Sum;
      acc.numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      acc.numGlobalCalcsN3 += numGlobalCalcsN3Sum;
      acc.numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
    }
  }

  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    const size_t n1 = soa1.size();
    const size_t n2 = soa2.size();
    if (n1 == 0 || n2 == 0) return;

    auto ptrs1 = SoAPointers::extract(soa1);
    auto ptrs2 = SoAPointers::extract(soa2);
    GlobalsAccumulator acc;

    const bool isPairHomogeneous =
        not useMixing or
        (ptrs1.type[0] == ptrs2.type[0] and soaIsHomogeneous(ptrs1.type, n1) and soaIsHomogeneous(ptrs2.type, n2));

    if (not useMixing or isPairHomogeneous) {
      const auto params = getHomogeneousParams(ptrs1.type);
      computeCellInteractions<false, true, newton3>(n1, n2, ptrs1, ptrs2, params, acc);
    } else if constexpr (useMixing) {
      MixingBuffers<512> mixingBuffers1;
      MixingBuffers<512> mixingBuffers2;

      const auto mixingPtrs1 =
          mixingBuffers1.template getMixingPointers<0>(ptrs1.type, n1, _halfSigmaPerType, _sqrtEpsilonPerType);
      ptrs1.halfSigma = mixingPtrs1.halfSigmaPtr;
      ptrs1.sqrtEpsilon = mixingPtrs1.sqrtEpsilonPtr;

      const auto mixingPtrs2 =
          mixingBuffers2.template getMixingPointers<1>(ptrs2.type, n2, _halfSigmaPerType, _sqrtEpsilonPerType);
      ptrs2.halfSigma = mixingPtrs2.halfSigmaPtr;
      ptrs2.sqrtEpsilon = mixingPtrs2.sqrtEpsilonPtr;

      computeCellInteractions<false, false, newton3>(n1, n2, ptrs1, ptrs2, {}, acc);
    }

    flushThreadData(acc);
  }

 public:
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorVerlet()
   */
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        std::span<const size_t> neighborList, bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

  /**
   * Sets the particle properties constants for this functor.
   *
   * This is only necessary if no particlePropertiesLibrary is used.
   *
   * @param epsilon24
   * @param sigmaSquared
   */
  void setParticleProperties(SoAFloatPrecision epsilon24, SoAFloatPrecision sigmaSquared) {
    _epsilon24 = epsilon24;
    _sigmaSquared = sigmaSquared;
    if constexpr (applyShift) {
      _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(_epsilon24, _sigmaSquared, _cutoffSquared);
    } else {
      _shift6 = 0.;
    }
  }

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
    if constexpr (calculateGlobals) {
      for (auto &data : _aosThreadDataGlobals) {
        data.setZero();
      }
    }
    if constexpr (countFLOPs) {
      for (auto &data : _aosThreadDataFLOPs) {
        data.setZero();
      }
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
    if constexpr (calculateGlobals) {
      for (const auto &data : _aosThreadDataGlobals) {
        _potentialEnergySum += data.potentialEnergySum;
        _virialSum += data.virialSum;
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
   * Get the potential Energy.
   * @return the potential Energy
   */
  double getPotentialEnergy() {
    if constexpr (not calculateGlobals) {
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
   * Get the virial.
   * @return
   */
  double getVirial() {
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
   * Gets the number of useful FLOPs.
   *
   * For the distance squared calculation, this is:
   * - Displacement: 3
   * - DistanceSquared in each dimension: 3
   * - DistanceSquared Total: 2
   * - Total: 8
   *
   * For the force kernel, this is:
   * - inverse distance squared: 1 (assume division is 1 FLOP)
   * - lj2: 1
   * - lj6: 2
   * - lj12: 1
   * - lj12m6: 1
   * - scalar factor: 3
   * - force: 3
   * - accumulate force on mol i: 3
   * - accumulate force on mol j if n3: 3
   * - Total: 15 without n3, 18 with n3
   *
   * For the globals calculation, this is:
   * - virial: 3
   * - potential: 1, or 2 with shift
   * - accumulation: 4 without n3, 8 with n3
   * - Total: 8 or 9 without n3, 12 or 13 with n3
   *
   * Caveats:
   *
   * This function is supposed to return useful FLOPs, e.g. without counting masked vector instructions.
   * You could also argue that, strictly speaking, we redundantly calculate forces and globals twice in the newton3 case
   * on a owned/halo boundary. This function does not treat such "redundant" calculations as useless. Similarly, this
   * function does not treat halo-halo interactions as redundant useless calculations.
   *
   * @return number of FLOPs since initTraversal() is called.
   */
  [[nodiscard]] size_t getNumFLOPs() const override {
    if constexpr (countFLOPs) {
      const size_t numDistCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numDistCalls; });
      const size_t numKernelCallsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
      const size_t numKernelCallsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });
      const size_t numGlobalCalcsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numGlobalCalcsN3; });
      const size_t numGlobalCalcsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numGlobalCalcsNoN3; });

      constexpr size_t numFLOPsPerDistanceCall = 8;
      constexpr size_t numFLOPsPerN3KernelCall = 18;
      constexpr size_t numFLOPsPerNoN3KernelCall = 15;
      constexpr size_t numFLOPsPerN3GlobalCalc = applyShift ? 13 : 12;
      constexpr size_t numFLOPsPerNoN3GlobalCalc = applyShift ? 9 : 8;

      return numDistCallsAcc * numFLOPsPerDistanceCall + numKernelCallsN3Acc * numFLOPsPerN3KernelCall +
             numKernelCallsNoN3Acc * numFLOPsPerNoN3KernelCall + numGlobalCalcsN3Acc * numFLOPsPerN3GlobalCalc +
             numGlobalCalcsNoN3Acc * numFLOPsPerNoN3GlobalCalc;
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<size_t>::max();
    }
  }

  [[nodiscard]] double getHitRate() const override {
    if constexpr (countFLOPs) {
      const size_t numDistCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numDistCalls; });
      const size_t numKernelCallsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
      const size_t numKernelCallsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });

      return (static_cast<double>(numKernelCallsNoN3Acc) + static_cast<double>(numKernelCallsN3Acc)) /
             (static_cast<double>(numDistCallsAcc));
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

 private:
  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            std::span<const size_t> neighborList) {
    const auto p = SoAPointers::extract(soa, _halfSigmaPerType.data(), _sqrtEpsilonPerType.data());

    // checks whether particle i is owned.
    const auto ownedStateI = p.owned[indexFirst];
    if (ownedStateI == autopas::OwnershipState::dummy) {
      return;
    }

    const SoAFloatPrecision xi = p.x[indexFirst];
    const SoAFloatPrecision yi = p.y[indexFirst];
    const SoAFloatPrecision zi = p.z[indexFirst];
    [[maybe_unused]] const SoAFloatPrecision ownedFactorI = (ownedStateI == autopas::OwnershipState::owned) ? 1. : 0.;

    [[maybe_unused]] SoAFloatPrecision halfSigmaI = 0.;
    [[maybe_unused]] SoAFloatPrecision sqrtEpsilon24I = 0.;
    if constexpr (useMixing) {
      const size_t typeI = p.type[indexFirst];
      halfSigmaI = p.halfSigma[typeI];
      sqrtEpsilon24I = 24.0 * p.sqrtEpsilon[typeI];
    }

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    [[maybe_unused]] const SoAFloatPrecision invCutoffSquared = _invCutoffSquared;

    GlobalsAccumulator acc;
    SoAFloatPrecision fxacc = 0;
    SoAFloatPrecision fyacc = 0;
    SoAFloatPrecision fzacc = 0;

#ifdef __AVX512F__
    constexpr size_t vecsize = 16;
#else
    constexpr size_t vecsize = 12;
#endif
    size_t joff = 0;

    const size_t neighborListSize = neighborList.size();
    if (neighborListSize >= vecsize) {
      alignas(64) std::array<SoAFloatPrecision, vecsize> xArr, yArr, zArr, fxArr, fyArr, fzArr;
      [[maybe_unused]] alignas(64) std::array<SoAFloatPrecision, vecsize> halfSigmaArr, sqrtEpsilonArr;
      alignas(64) std::array<autopas::OwnershipState, vecsize> ownedStateArr{};

      for (; joff < neighborListSize - vecsize + 1; joff += vecsize) {
#pragma omp simd safelen(vecsize)
        for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
          const size_t j = neighborList[joff + tmpj];
          xArr[tmpj] = p.x[j];
          yArr[tmpj] = p.y[j];
          zArr[tmpj] = p.z[j];
          ownedStateArr[tmpj] = p.owned[j];
          if constexpr (useMixing) {
            const size_t typeJ = p.type[j];
            halfSigmaArr[tmpj] = p.halfSigma[typeJ];
            sqrtEpsilonArr[tmpj] = p.sqrtEpsilon[typeJ];
          }
        }

        if constexpr (not calculateGlobals and not countFLOPs) {
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc) safelen(vecsize)
          for (size_t j = 0; j < vecsize; j++) {
            SoAFloatPrecision sigmaSquared = _sigmaSquared;
            SoAFloatPrecision epsilon24 = _epsilon24;
            if constexpr (useMixing) {
              const SoAFloatPrecision sigma = halfSigmaI + halfSigmaArr[j];
              sigmaSquared = sigma * sigma;
              epsilon24 = sqrtEpsilon24I * sqrtEpsilonArr[j];
            }

            const auto ownedStateJ = ownedStateArr[j];
            const SoAFloatPrecision drx = xi - xArr[j];
            const SoAFloatPrecision dry = yi - yArr[j];
            const SoAFloatPrecision drz = zi - zArr[j];
            const SoAFloatPrecision dr2 = drx * drx + dry * dry + drz * drz;
            const bool mask = (dr2 <= cutoffSquared) & (ownedStateJ != autopas::OwnershipState::dummy);

            SoAFloatPrecision fx, fy, fz, lj12m6;
            computeLJForceCore(drx, dry, drz, dr2, mask, sigmaSquared, epsilon24, fx, fy, fz, lj12m6);

            fxacc += fx;
            fyacc += fy;
            fzacc += fz;
            if constexpr (newton3) {
              fxArr[j] = fx;
              fyArr[j] = fy;
              fzArr[j] = fz;
            }
          }
        } else {
          SoAFloatPrecision potentialEnergySum = 0.;
          SoAFloatPrecision virialSumX = 0.;
          SoAFloatPrecision virialSumY = 0.;
          SoAFloatPrecision virialSumZ = 0.;
          size_t numDistanceCalculationSum = 0;
          size_t numKernelCallsN3Sum = 0;
          size_t numKernelCallsNoN3Sum = 0;
          size_t numGlobalCalcsN3Sum = 0;
          size_t numGlobalCalcsNoN3Sum = 0;

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ, \
                               numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum,       \
                               numGlobalCalcsN3Sum, numGlobalCalcsNoN3Sum) safelen(vecsize)
          for (size_t j = 0; j < vecsize; j++) {
            SoAFloatPrecision sigmaSquared = _sigmaSquared;
            SoAFloatPrecision epsilon24 = _epsilon24;
            SoAFloatPrecision shift6 = _shift6;
            if constexpr (useMixing) {
              const SoAFloatPrecision sigma = halfSigmaI + halfSigmaArr[j];
              sigmaSquared = sigma * sigma;
              epsilon24 = sqrtEpsilon24I * sqrtEpsilonArr[j];
              if constexpr (calculateGlobals and applyShift) {
                shift6 = computeShift6(epsilon24, sigmaSquared, invCutoffSquared);
              }
            }

            const auto ownedStateJ = ownedStateArr[j];
            const SoAFloatPrecision drx = xi - xArr[j];
            const SoAFloatPrecision dry = yi - yArr[j];
            const SoAFloatPrecision drz = zi - zArr[j];
            const SoAFloatPrecision dr2 = drx * drx + dry * dry + drz * drz;
            const bool mask = (dr2 <= cutoffSquared) & (ownedStateJ != autopas::OwnershipState::dummy);

            SoAFloatPrecision fx, fy, fz, lj12m6;
            computeLJForceCore(drx, dry, drz, dr2, mask, sigmaSquared, epsilon24, fx, fy, fz, lj12m6);

            fxacc += fx;
            fyacc += fy;
            fzacc += fz;
            if constexpr (newton3) {
              fxArr[j] = fx;
              fyArr[j] = fy;
              fzArr[j] = fz;
            }

            if constexpr (countFLOPs) {
              numDistanceCalculationSum += (ownedStateJ != autopas::OwnershipState::dummy) ? 1 : 0;
              if constexpr (newton3) {
                numKernelCallsN3Sum += mask ? 1 : 0;
              } else {
                numKernelCallsNoN3Sum += mask ? 1 : 0;
              }
            }

            if constexpr (calculateGlobals) {
              const SoAFloatPrecision virialx = drx * fx;
              const SoAFloatPrecision virialy = dry * fy;
              const SoAFloatPrecision virialz = drz * fz;
              const SoAFloatPrecision potentialEnergy6 = mask ? (epsilon24 * lj12m6 + shift6) : 0.;

              SoAFloatPrecision energyFactor = ownedFactorI;
              if constexpr (newton3) {
                energyFactor += (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.);
              }

              potentialEnergySum += potentialEnergy6 * energyFactor;
              virialSumX += virialx * energyFactor;
              virialSumY += virialy * energyFactor;
              virialSumZ += virialz * energyFactor;

              if constexpr (countFLOPs) {
                if constexpr (newton3) {
                  numGlobalCalcsN3Sum += mask ? 1 : 0;
                } else {
                  numGlobalCalcsNoN3Sum += mask ? 1 : 0;
                }
              }
            }
          }

          acc.potentialEnergySum += potentialEnergySum;
          acc.virialX += virialSumX;
          acc.virialY += virialSumY;
          acc.virialZ += virialSumZ;
          acc.numDistCalls += numDistanceCalculationSum;
          acc.numKernelCallsN3 += numKernelCallsN3Sum;
          acc.numKernelCallsNoN3 += numKernelCallsNoN3Sum;
          acc.numGlobalCalcsN3 += numGlobalCalcsN3Sum;
          acc.numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
        }

        if constexpr (newton3) {
#pragma omp simd safelen(vecsize)
          for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
            const size_t j = neighborList[joff + tmpj];
            p.fx[j] -= fxArr[tmpj];
            p.fy[j] -= fyArr[tmpj];
            p.fz[j] -= fzArr[tmpj];
          }
        }
      }
    }

    // Remainder loop
    for (size_t jNeighIndex = joff; jNeighIndex < neighborListSize; ++jNeighIndex) {
      const size_t j = neighborList[jNeighIndex];
      if (indexFirst == j) continue;

      const auto ownedStateJ = p.owned[j];
      if (ownedStateJ == autopas::OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision drx = xi - p.x[j];
      const SoAFloatPrecision dry = yi - p.y[j];
      const SoAFloatPrecision drz = zi - p.z[j];
      const SoAFloatPrecision dr2 = drx * drx + dry * dry + drz * drz;

      if constexpr (countFLOPs) {
        acc.numDistCalls += 1;
      }

      if (dr2 > cutoffSquared) {
        continue;
      }

      SoAFloatPrecision sigmaSquared = _sigmaSquared;
      SoAFloatPrecision epsilon24 = _epsilon24;
      SoAFloatPrecision shift6 = _shift6;
      if constexpr (useMixing) {
        const size_t typeJ = p.type[j];
        const SoAFloatPrecision sigma = halfSigmaI + p.halfSigma[typeJ];
        sigmaSquared = sigma * sigma;
        epsilon24 = sqrtEpsilon24I * p.sqrtEpsilon[typeJ];
        if constexpr (calculateGlobals and applyShift) {
          shift6 = computeShift6(epsilon24, sigmaSquared, invCutoffSquared);
        }
      }

      SoAFloatPrecision fx, fy, fz, lj12m6;
      computeLJForceCore(drx, dry, drz, dr2, true, sigmaSquared, epsilon24, fx, fy, fz, lj12m6);

      fxacc += fx;
      fyacc += fy;
      fzacc += fz;
      if constexpr (newton3) {
        p.fx[j] -= fx;
        p.fy[j] -= fy;
        p.fz[j] -= fz;
      }

      if constexpr (countFLOPs) {
        if constexpr (newton3) {
          acc.numKernelCallsN3 += 1;
        } else {
          acc.numKernelCallsNoN3 += 1;
        }
      }

      if constexpr (calculateGlobals) {
        const SoAFloatPrecision virialx = drx * fx;
        const SoAFloatPrecision virialy = dry * fy;
        const SoAFloatPrecision virialz = drz * fz;
        const SoAFloatPrecision potentialEnergy6 = (epsilon24 * lj12m6 + shift6);

        const SoAFloatPrecision energyFactor =
            ownedFactorI + (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.);
        acc.potentialEnergySum += potentialEnergy6 * energyFactor;
        acc.virialX += virialx * energyFactor;
        acc.virialY += virialy * energyFactor;
        acc.virialZ += virialz * energyFactor;

        if constexpr (countFLOPs) {
          if constexpr (newton3) {
            acc.numGlobalCalcsN3 += 1;
          } else {
            acc.numGlobalCalcsNoN3 += 1;
          }
        }
      }
    }

    if (fxacc != 0 or fyacc != 0 or fzacc != 0) {
      p.fx[indexFirst] += fxacc;
      p.fy[indexFirst] += fyacc;
      p.fz[indexFirst] += fzacc;
    }

    flushThreadData(acc);
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
      _sigmaSquared = 4.0 * _halfSigmaPerType[0] * _halfSigmaPerType[0];
      _epsilon24 = 24.0 * _sqrtEpsilonPerType[0] * _sqrtEpsilonPerType[0];
      if constexpr (applyShift) {
        _shift6 = ppl.getMixingShift6(0, 0);
      }
    }
  }

  static inline SoAFloatPrecision computeShift6(const SoAFloatPrecision epsilon24, const SoAFloatPrecision sigmaSquared,
                                                const SoAFloatPrecision invCutoffSquared) {
    const SoAFloatPrecision sig2Cut2 = sigmaSquared * invCutoffSquared;
    const SoAFloatPrecision sig6Cut6 = sig2Cut2 * sig2Cut2 * sig2Cut2;
    const SoAFloatPrecision sig12Cut12 = sig6Cut6 * sig6Cut6;
    return epsilon24 * (sig6Cut6 - sig12Cut12);
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
   * This class stores internal data for global calculations for each thread. Make sure that this data has proper size,
   * i.e. k*64 Bytes!
   */
  class AoSThreadDataGlobals {
   public:
    AoSThreadDataGlobals() : virialSum{0., 0., 0.}, potentialEnergySum{0.}, __remainingTo64{} {}
    void setZero() {
      virialSum = {0., 0., 0.};
      potentialEnergySum = 0.;
    }

    // variables
    std::array<double, 3> virialSum;
    double potentialEnergySum;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };

  /**
   * This class stores internal data for FLOP counters for each thread. Make sure that this data has proper size, i.e.
   * k*64 Bytes!
   * The FLOP count and HitRate are not counted/calculated directly, but through helper counters (numKernelCallsNoN3,
   * numKernelCallsN3, numDistCalls, numGlobalCalcs) to reduce computational cost in the functors themselves and to
   * improve maintainability (e.g. if the cost of a kernel call changes).
   */
  class AoSThreadDataFLOPs {
   public:
    AoSThreadDataFLOPs() : __remainingTo64{} {}

    /**
     * Set all counters to zero.
     */
    void setZero() {
      numKernelCallsNoN3 = 0;
      numKernelCallsN3 = 0;
      numDistCalls = 0;
      numGlobalCalcsNoN3 = 0;
      numGlobalCalcsN3 = 0;
    }

    /**
     * Number of calls to Lennard-Jones Kernel with newton3 disabled.
     * Used for calculating number of FLOPs and hit rate.
     */
    size_t numKernelCallsNoN3 = 0;

    /**
     * Number of calls to Lennard-Jones Kernel with newton3 enabled.
     * Used for calculating number of FLOPs and hit rate.
     */
    size_t numKernelCallsN3 = 0;

    /**
     * Number of distance calculations.
     * Used for calculating number of FLOPs and hit rate.
     */
    size_t numDistCalls = 0;

    /**
     * Counter for the number of times the globals have been calculated. Excludes the special case that N3 is enabled
     * and we calculate globals for an owned-halo particle pair.
     */
    size_t numGlobalCalcsN3 = 0;

    /**
     * Counter for the number of times the globals have been calculated. Excludes the special case that N3 is enabled
     * and we calculate globals for an owned-halo particle pair.
     */
    size_t numGlobalCalcsNoN3 = 0;

   private:
    /**
     * dummy parameter to get the right size (64 bytes)
     */
    double __remainingTo64[(64 - 5 * sizeof(size_t)) / sizeof(size_t)];
  };

  // make sure of the size of AoSThreadDataGlobals and AoSThreadDataFLOPs
  static_assert(sizeof(AoSThreadDataGlobals) % 64 == 0, "AoSThreadDataGlobals has wrong size");
  static_assert(sizeof(AoSThreadDataFLOPs) % 64 == 0, "AoSThreadDataFLOPs has wrong size");

  const double _cutoffSquared;
  const double _invCutoffSquared;
  // not const because they might be reset through PPL
  double _epsilon24, _sigmaSquared, _shift6 = 0;

  std::vector<SoAFloatPrecision> _halfSigmaPerType{};
  std::vector<SoAFloatPrecision> _sqrtEpsilonPerType{};

  // sum of the potential energy, only calculated if calculateGlobals is true
  double _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<double, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals{};
  std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;
};
}  // namespace mdLib
