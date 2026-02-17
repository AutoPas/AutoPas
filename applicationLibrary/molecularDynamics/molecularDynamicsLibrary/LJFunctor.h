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
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

namespace mdLib {

/**
 * A functor to handle lennard-jones interactions between two particles (molecules).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
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
    _PPLibrary = &particlePropertiesLibrary;
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

    const auto threadnum = autopas::autopas_get_thread_num();

    if constexpr (countFLOPs) {
      ++_aosThreadDataFLOPs[threadnum].numDistCalls;
    }

    auto sigmaSquared = _sigmaSquared;
    auto epsilon24 = _epsilon24;
    auto shift6 = _shift6;
    if constexpr (useMixing) {
      sigmaSquared = _PPLibrary->getMixingSigmaSquared(i.getTypeId(), j.getTypeId());
      epsilon24 = _PPLibrary->getMixing24Epsilon(i.getTypeId(), j.getTypeId());
      if constexpr (applyShift) {
        shift6 = _PPLibrary->getMixingShift6(i.getTypeId(), j.getTypeId());
      }
    }
    auto dr = i.getR() - j.getR();
    double dr2 = autopas::utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffSquared) {
      return;
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

    if constexpr (countFLOPs) {
      if (newton3) {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsN3;
      } else {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsNoN3;
      }
    }

    if constexpr (calculateGlobals) {
      // We always add the full contribution for each owned particle and divide the sums by 2 in endTraversal().
      // Potential energy has an additional factor of 6, which is also handled in endTraversal().

      auto virial = dr * f;
      double potentialEnergy6 = epsilon24 * lj12m6 + shift6;

      if (i.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadDataGlobals[threadnum].virialSum += virial;
      }
      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergy6;
        _aosThreadDataGlobals[threadnum].virialSum += virial;
      }
      if constexpr (countFLOPs) {
        if (newton3) {
          ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsN3;
        } else {
          ++_aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3;
        }
      }
    }
  }

  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversal of the soa.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (soa.size() == 0) return;

    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle_T::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle_T::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle_T::AttributeNames::typeId>();
    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;

    SoAFloatPrecision potentialEnergySum = 0.;  // Note: This is not the potential energy but some fixed multiple of it.
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsSum = 0;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmaSquareds;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> epsilon24s;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> shift6s;
    if constexpr (useMixing) {
      // Preload all sigma and epsilons for next vectorized region.
      // Not preloading and directly using the values, will produce worse results.
      sigmaSquareds.resize(soa.size());
      epsilon24s.resize(soa.size());
      // if no mixing or mixing but no shift shift6 is constant therefore we do not need this vector.
      if constexpr (applyShift) {
        shift6s.resize(soa.size());
      }
    }

    const SoAFloatPrecision const_shift6 = _shift6;
    const SoAFloatPrecision const_sigmaSquared = _sigmaSquared;
    const SoAFloatPrecision const_epsilon24 = _epsilon24;

    for (unsigned int i = 0; i < soa.size(); ++i) {
      const auto ownedStateI = ownedStatePtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      SoAFloatPrecision fxacc = 0.;
      SoAFloatPrecision fyacc = 0.;
      SoAFloatPrecision fzacc = 0.;

      if constexpr (useMixing) {
        for (unsigned int j = 0; j < soa.size(); ++j) {
          auto mixingData = _PPLibrary->getLJMixingData(typeptr[i], typeptr[j]);
          sigmaSquareds[j] = mixingData.sigmaSquared;
          epsilon24s[j] = mixingData.epsilon24;
          if constexpr (applyShift) {
            shift6s[j] = mixingData.shift6;
          }
        }
      }

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ, numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, numGlobalCalcsSum)
      for (unsigned int j = i + 1; j < soa.size(); ++j) {
        SoAFloatPrecision shift6 = const_shift6;
        SoAFloatPrecision sigmaSquared = const_sigmaSquared;
        SoAFloatPrecision epsilon24 = const_epsilon24;
        if constexpr (useMixing) {
          sigmaSquared = sigmaSquareds[j];
          epsilon24 = epsilon24s[j];
          if constexpr (applyShift) {
            shift6 = shift6s[j];
          }
        }

        const auto ownedStateJ = ownedStatePtr[j];

        const SoAFloatPrecision drx = xptr[i] - xptr[j];
        const SoAFloatPrecision dry = yptr[i] - yptr[j];
        const SoAFloatPrecision drz = zptr[i] - zptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;

        // Mask away if distance is too large or any particle is a dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

        const SoAFloatPrecision invdr2 = 1. / dr2;
        const SoAFloatPrecision lj2 = sigmaSquared * invdr2;
        const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
        const SoAFloatPrecision lj12 = lj6 * lj6;
        const SoAFloatPrecision lj12m6 = lj12 - lj6;
        const SoAFloatPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * invdr2;

        const SoAFloatPrecision fx = drx * fac;
        const SoAFloatPrecision fy = dry * fac;
        const SoAFloatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        // newton 3
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;

        if constexpr (countFLOPs) {
          numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
          numKernelCallsN3Sum += mask;
        }

        if (calculateGlobals) {
          const SoAFloatPrecision virialx = drx * fx;
          const SoAFloatPrecision virialy = dry * fy;
          const SoAFloatPrecision virialz = drz * fz;
          const SoAFloatPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

          // We add 6 times the potential energy for each owned particle. The total sum is corrected in endTraversal().
          SoAFloatPrecision energyFactor = (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
                                           (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.);
          potentialEnergySum += potentialEnergy6 * energyFactor;

          virialSumX += virialx * energyFactor;
          virialSumY += virialy * energyFactor;
          virialSumZ += virialz * energyFactor;

          if constexpr (countFLOPs) {
            numGlobalCalcsSum += mask;
          }
        }
      }

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsSum;  // Always N3 in Single SoAFunctor
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
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
  /**
   * Implementation function of SoAFunctorPair(soa1, soa2, newton3)
   *
   * @tparam newton3
   * @param soa1
   * @param soa2
   */
  template <bool newton3>
  void SoAFunctorPairImpl(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2) {
    if (soa1.size() == 0 || soa2.size() == 0) return;

    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict x1ptr = soa1.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y1ptr = soa1.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z1ptr = soa1.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict x2ptr = soa2.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict y2ptr = soa2.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict z2ptr = soa2.template begin<Particle_T::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle_T::AttributeNames::ownershipState>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle_T::AttributeNames::ownershipState>();

    auto *const __restrict fx1ptr = soa1.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy1ptr = soa1.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz1ptr = soa1.template begin<Particle_T::AttributeNames::forceZ>();
    auto *const __restrict fx2ptr = soa2.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fy2ptr = soa2.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fz2ptr = soa2.template begin<Particle_T::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle_T::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle_T::AttributeNames::typeId>();

    // Checks whether the cells are halo cells.
    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    SoAFloatPrecision shift6 = _shift6;
    SoAFloatPrecision sigmaSquared = _sigmaSquared;
    SoAFloatPrecision epsilon24 = _epsilon24;

    // preload all sigma and epsilons for next vectorized region
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmaSquareds;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> epsilon24s;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> shift6s;
    if constexpr (useMixing) {
      sigmaSquareds.resize(soa2.size());
      epsilon24s.resize(soa2.size());
      // if no mixing or mixing but no shift shift6 is constant therefore we do not need this vector.
      if constexpr (applyShift) {
        shift6s.resize(soa2.size());
      }
    }

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      SoAFloatPrecision fxacc = 0;
      SoAFloatPrecision fyacc = 0;
      SoAFloatPrecision fzacc = 0;

      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      // preload all sigma and epsilons for next vectorized region
      if constexpr (useMixing) {
        for (unsigned int j = 0; j < soa2.size(); ++j) {
          sigmaSquareds[j] = _PPLibrary->getMixingSigmaSquared(typeptr1[i], typeptr2[j]);
          epsilon24s[j] = _PPLibrary->getMixing24Epsilon(typeptr1[i], typeptr2[j]);
          if constexpr (applyShift) {
            shift6s[j] = _PPLibrary->getMixingShift6(typeptr1[i], typeptr2[j]);
          }
        }
      }

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ, numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, numGlobalCalcsN3Sum, numGlobalCalcsNoN3Sum)
      for (unsigned int j = 0; j < soa2.size(); ++j) {
        if constexpr (useMixing) {
          sigmaSquared = sigmaSquareds[j];
          epsilon24 = epsilon24s[j];
          if constexpr (applyShift) {
            shift6 = shift6s[j];
          }
        }

        const auto ownedStateJ = ownedStatePtr2[j];

        const SoAFloatPrecision drx = x1ptr[i] - x2ptr[j];
        const SoAFloatPrecision dry = y1ptr[i] - y2ptr[j];
        const SoAFloatPrecision drz = z1ptr[i] - z2ptr[j];

        const SoAFloatPrecision drx2 = drx * drx;
        const SoAFloatPrecision dry2 = dry * dry;
        const SoAFloatPrecision drz2 = drz * drz;

        const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;

        // Mask away if distance is too large or any particle is a dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

        const SoAFloatPrecision invdr2 = 1. / dr2;
        const SoAFloatPrecision lj2 = sigmaSquared * invdr2;
        const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
        const SoAFloatPrecision lj12 = lj6 * lj6;
        const SoAFloatPrecision lj12m6 = lj12 - lj6;
        const SoAFloatPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * invdr2;

        const SoAFloatPrecision fx = drx * fac;
        const SoAFloatPrecision fy = dry * fac;
        const SoAFloatPrecision fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;
        if (newton3) {
          fx2ptr[j] -= fx;
          fy2ptr[j] -= fy;
          fz2ptr[j] -= fz;
        }

        if constexpr (countFLOPs) {
          numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
          if constexpr (newton3) {
            numKernelCallsN3Sum += mask;
          } else {
            numKernelCallsNoN3Sum += mask;
          }
        }

        if constexpr (calculateGlobals) {
          SoAFloatPrecision virialx = drx * fx;
          SoAFloatPrecision virialy = dry * fy;
          SoAFloatPrecision virialz = drz * fz;
          SoAFloatPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

          // We add 6 times the potential energy for each owned particle. The total sum is corrected in endTraversal().
          const SoAFloatPrecision energyFactor =
              (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
              (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.);
          potentialEnergySum += potentialEnergy6 * energyFactor;
          virialSumX += virialx * energyFactor;
          virialSumY += virialy * energyFactor;
          virialSumZ += virialz * energyFactor;

          if constexpr (countFLOPs) {
            if constexpr (newton3) {
              numGlobalCalcsN3Sum += mask;
            } else {
              numGlobalCalcsNoN3Sum += mask;
            }
          }
        }
      }
      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;
    }
    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

 public:
  // clang-format off
  /**
   * @copydoc autopas::PairwiseFunctor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<autopas::SoAIndexIntType, autopas::AlignedAllocator<autopas::SoAIndexIntType>> &neighborList,
                        bool newton3) final {
    if (soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletImpl<false>(soa, indexFirst, neighborList);
    }
  }

  void SoAFunctorVerletOptimized(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                      const std::vector<autopas::SoAIndexIntType, autopas::AlignedAllocator<autopas::SoAIndexIntType>> &neighborList,
                      bool newton3) final{
    if (soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletOptimizedImpl<true>(soa, indexFirst, neighborList);
    } else {
      SoAFunctorVerletOptimizedImpl<false>(soa, indexFirst, neighborList);
    }
  }

  void SoAFunctorVerletOptimizedCompactAoS(autopas::VerletListsLJCompactAoS<Particle_T> &compactAoS, const size_t indexFirst,
                    const std::vector<autopas::SoAIndexIntType, autopas::AlignedAllocator<autopas::SoAIndexIntType>> &neighborList,
                    bool newton3) final{
    if (compactAoS._soa.size() == 0 or neighborList.empty()) return;
    if (newton3) {
      SoAFunctorVerletOptimizedCompactAoSImpl<true>(compactAoS, indexFirst, neighborList);
    } else {
      SoAFunctorVerletOptimizedCompactAoSImpl<false>(compactAoS, indexFirst, neighborList);
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
    if (applyShift) {
      _shift6 = ParticlePropertiesLibrary<double, size_t>::calcShift6(_epsilon24, _sigmaSquared, _cutoffSquared);
    } else {
      _shift6 = 0.;
    }
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle_T::AttributeNames, 8>{Particle_T::AttributeNames::posX,
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
    if (calculateGlobals) {
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
   * Get the virial.
   * @return
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
                            const std::vector<autopas::SoAIndexIntType, autopas::AlignedAllocator<autopas::SoAIndexIntType>> &neighborList) {
    const auto *const __restrict xptr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle_T::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle_T::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa.template begin<Particle_T::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa.template begin<Particle_T::AttributeNames::typeId>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    SoAFloatPrecision shift6 = _shift6;
    SoAFloatPrecision sigmaSquared = _sigmaSquared;
    SoAFloatPrecision epsilon24 = _epsilon24;

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    // Counters for when countFLOPs is activated
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    SoAFloatPrecision fxacc = 0;
    SoAFloatPrecision fyacc = 0;
    SoAFloatPrecision fzacc = 0;
    const size_t neighborListSize = neighborList.size();
    const autopas::SoAIndexIntType *const __restrict neighborListPtr = neighborList.data();

    // checks whether particle i is owned.
    const auto ownedStateI = ownedStatePtr[indexFirst];
    if (ownedStateI == autopas::OwnershipState::dummy) {
      return;
    }

    const auto threadnum = autopas::autopas_get_thread_num();

    // this is a magic number, that should correspond to at least
    // vectorization width*N have testet multiple sizes:
    // 4: does not give a speedup, slower than original AoSFunctor
    // 8: small speedup compared to AoS
    // 12: highest speedup compared to Aos
    // 16: smaller speedup
    // in theory this is a variable, we could auto-tune over...
#ifdef __AVX512F__
    // use a multiple of 8 for avx
    constexpr size_t vecsize = 16;
#else
    // for everything else 12 is faster
    constexpr size_t vecsize = 12;
#endif
    size_t joff = 0;

    // if the size of the verlet list is larger than the given size vecsize,
    // we will use a vectorized version.
    if (neighborListSize >= vecsize) {
      alignas(64) std::array<SoAFloatPrecision, vecsize> xtmp, ytmp, ztmp, xArr, yArr, zArr, fxArr, fyArr, fzArr;
      alignas(64) std::array<autopas::OwnershipState, vecsize> ownedStateArr{};
      // broadcast of the position of particle i
      for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
        xtmp[tmpj] = xptr[indexFirst];
        ytmp[tmpj] = yptr[indexFirst];
        ztmp[tmpj] = zptr[indexFirst];
      }
      // loop over the verlet list from 0 to x*vecsize
      for (; joff < neighborListSize - vecsize + 1; joff += vecsize) {
        // in each iteration we calculate the interactions of particle i with
        // vecsize particles in the neighborlist of particle i starting at
        // particle joff

        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> sigmaSquareds;
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> epsilon24s;
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecsize> shift6s;
        if constexpr (useMixing) {
          for (size_t j = 0; j < vecsize; j++) {
            sigmaSquareds[j] =
                _PPLibrary->getMixingSigmaSquared(typeptr1[indexFirst], typeptr2[neighborListPtr[joff + j]]);
            epsilon24s[j] = _PPLibrary->getMixing24Epsilon(typeptr1[indexFirst], typeptr2[neighborListPtr[joff + j]]);
            if constexpr (applyShift) {
              shift6s[j] = _PPLibrary->getMixingShift6(typeptr1[indexFirst], typeptr2[neighborListPtr[joff + j]]);
            }
          }
        }

        // gather position of particle j
#pragma omp simd safelen(vecsize)
        for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
          xArr[tmpj] = xptr[neighborListPtr[joff + tmpj]];
          yArr[tmpj] = yptr[neighborListPtr[joff + tmpj]];
          zArr[tmpj] = zptr[neighborListPtr[joff + tmpj]];
          ownedStateArr[tmpj] = ownedStatePtr[neighborListPtr[joff + tmpj]];
        }
        // do omp simd with reduction of the interaction
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, potentialEnergySum, virialSumX, virialSumY, virialSumZ, numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, numGlobalCalcsN3Sum, numGlobalCalcsNoN3Sum) safelen(vecsize)
        for (size_t j = 0; j < vecsize; j++) {
          if constexpr (useMixing) {
            sigmaSquared = sigmaSquareds[j];
            epsilon24 = epsilon24s[j];
            if constexpr (applyShift) {
              shift6 = shift6s[j];
            }
          }
          // const size_t j = currentList[jNeighIndex];

          const auto ownedStateJ = ownedStateArr[j];

          const SoAFloatPrecision drx = xtmp[j] - xArr[j];
          const SoAFloatPrecision dry = ytmp[j] - yArr[j];
          const SoAFloatPrecision drz = ztmp[j] - zArr[j];

          const SoAFloatPrecision drx2 = drx * drx;
          const SoAFloatPrecision dry2 = dry * dry;
          const SoAFloatPrecision drz2 = drz * drz;

          const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;

          // Mask away if distance is too large or any particle is a dummy.
          // Particle ownedStateI was already checked previously.
          const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

          const SoAFloatPrecision invdr2 = 1. / dr2;
          const SoAFloatPrecision lj2 = sigmaSquared * invdr2;
          const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
          const SoAFloatPrecision lj12 = lj6 * lj6;
          const SoAFloatPrecision lj12m6 = lj12 - lj6;
          const SoAFloatPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * invdr2;

          const SoAFloatPrecision fx = drx * fac;
          const SoAFloatPrecision fy = dry * fac;
          const SoAFloatPrecision fz = drz * fac;

          fxacc += fx;
          fyacc += fy;
          fzacc += fz;
          if (newton3) {
            fxArr[j] = fx;
            fyArr[j] = fy;
            fzArr[j] = fz;
          }

          if constexpr (countFLOPs) {
            numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
            if constexpr (newton3) {
              numKernelCallsN3Sum += mask;
            } else {
              numKernelCallsNoN3Sum += mask;
            }
          }

          if (calculateGlobals) {
            SoAFloatPrecision virialx = drx * fx;
            SoAFloatPrecision virialy = dry * fy;
            SoAFloatPrecision virialz = drz * fz;
            SoAFloatPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

            // We add 6 times the potential energy for each owned particle. The total sum is corrected in
            // endTraversal().
            const SoAFloatPrecision energyFactor =
                (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
                (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.);
            potentialEnergySum += potentialEnergy6 * energyFactor;
            virialSumX += virialx * energyFactor;
            virialSumY += virialy * energyFactor;
            virialSumZ += virialz * energyFactor;

            if constexpr (countFLOPs) {
              if constexpr (newton3) {
                numGlobalCalcsN3Sum += mask;
              } else {
                numGlobalCalcsNoN3Sum += mask;
              }
            }
          }
        }
        // scatter the forces to where they belong, this is only needed for newton3
        if (newton3) {
#pragma omp simd safelen(vecsize)
          for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
            const size_t j = neighborListPtr[joff + tmpj];
            fxptr[j] -= fxArr[tmpj];
            fyptr[j] -= fyArr[tmpj];
            fzptr[j] -= fzArr[tmpj];
          }
        }
      }
    }
    // this loop goes over the remainder and uses no optimizations
    for (size_t jNeighIndex = joff; jNeighIndex < neighborListSize; ++jNeighIndex) {
      size_t j = neighborList[jNeighIndex];
      if (indexFirst == j) continue;
      if constexpr (useMixing) {
        sigmaSquared = _PPLibrary->getMixingSigmaSquared(typeptr1[indexFirst], typeptr2[j]);
        epsilon24 = _PPLibrary->getMixing24Epsilon(typeptr1[indexFirst], typeptr2[j]);
        if constexpr (applyShift) {
          shift6 = _PPLibrary->getMixingShift6(typeptr1[indexFirst], typeptr2[j]);
        }
      }

      const auto ownedStateJ = ownedStatePtr[j];
      if (ownedStateJ == autopas::OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision drx = xptr[indexFirst] - xptr[j];
      const SoAFloatPrecision dry = yptr[indexFirst] - yptr[j];
      const SoAFloatPrecision drz = zptr[indexFirst] - zptr[j];

      const SoAFloatPrecision drx2 = drx * drx;
      const SoAFloatPrecision dry2 = dry * dry;
      const SoAFloatPrecision drz2 = drz * drz;

      const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;

      if constexpr (countFLOPs) {
        numDistanceCalculationSum += 1;
      }

      if (dr2 > cutoffSquared) {
        continue;
      }

      const SoAFloatPrecision invdr2 = 1. / dr2;
      const SoAFloatPrecision lj2 = sigmaSquared * invdr2;
      const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
      const SoAFloatPrecision lj12 = lj6 * lj6;
      const SoAFloatPrecision lj12m6 = lj12 - lj6;
      const SoAFloatPrecision fac = epsilon24 * (lj12 + lj12m6) * invdr2;

      const SoAFloatPrecision fx = drx * fac;
      const SoAFloatPrecision fy = dry * fac;
      const SoAFloatPrecision fz = drz * fac;

      fxacc += fx;
      fyacc += fy;
      fzacc += fz;
      if (newton3) {
        fxptr[j] -= fx;
        fyptr[j] -= fy;
        fzptr[j] -= fz;
      }

      if constexpr (countFLOPs) {
        if constexpr (newton3) {
          numKernelCallsN3Sum += 1;
        } else {
          numKernelCallsNoN3Sum += 1;
        }
      }

      if (calculateGlobals) {
        SoAFloatPrecision virialx = drx * fx;
        SoAFloatPrecision virialy = dry * fy;
        SoAFloatPrecision virialz = drz * fz;
        SoAFloatPrecision potentialEnergy6 = (epsilon24 * lj12m6 + shift6);

        // We add 6 times the potential energy for each owned particle. The total sum is corrected in endTraversal().
        const SoAFloatPrecision energyFactor =
            (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
            (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.);
        potentialEnergySum += potentialEnergy6 * energyFactor;
        virialSumX += virialx * energyFactor;
        virialSumY += virialy * energyFactor;
        virialSumZ += virialz * energyFactor;

        if constexpr (countFLOPs) {
          if constexpr (newton3) {
            ++numGlobalCalcsN3Sum;
          } else {
            ++numGlobalCalcsNoN3Sum;
          }
        }
      }
    }

    if (fxacc != 0 or fyacc != 0 or fzacc != 0) {
      fxptr[indexFirst] += fxacc;
      fyptr[indexFirst] += fyacc;
      fzptr[indexFirst] += fzacc;
    }

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;
    }

    if (calculateGlobals) {
      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

#pragma code_align 32
  template <bool newton3>
  void SoAFunctorVerletOptimizedImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<autopas::SoAIndexIntType, autopas::AlignedAllocator<autopas::SoAIndexIntType>> &neighborList) {

    const auto *const __restrict xPtr = soa.template begin<Particle_T::AttributeNames::posX>();
    const auto *const __restrict yPtr = soa.template begin<Particle_T::AttributeNames::posY>();
    const auto *const __restrict zPtr = soa.template begin<Particle_T::AttributeNames::posZ>();

    auto *const __restrict forceXPtr = soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict forceYPtr = soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict forceZPtr = soa.template begin<Particle_T::AttributeNames::forceZ>();

    SoAFloatPrecision forceXAcc = 0;
    SoAFloatPrecision forceYAcc = 0;
    SoAFloatPrecision forceZAcc = 0;

    [[maybe_unused]] const auto *const __restrict typeIdPtr = soa.template begin<Particle_T::AttributeNames::typeId>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle_T::AttributeNames::ownershipState>();

    const size_t neighborListSize = neighborList.size();
    const autopas::SoAIndexIntType *const __restrict neighborListPtr = neighborList.data();

    const SoAFloatPrecision xI = xPtr[indexFirst];
    const SoAFloatPrecision yI = yPtr[indexFirst];
    const SoAFloatPrecision zI = zPtr[indexFirst];

    const size_t typeIdI = typeIdPtr[indexFirst];

    // checks whether particle i is owned.
    const auto ownedStateI = ownedStatePtr[indexFirst];
    if (ownedStateI == autopas::OwnershipState::dummy) {
      return;
    }

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    SoAFloatPrecision shift6 = _shift6;
    SoAFloatPrecision sigmaSquared = _sigmaSquared;
    SoAFloatPrecision epsilon24 = _epsilon24;

    using PackedLJMixingData = ParticlePropertiesLibrary<SoAFloatPrecision, unsigned long>::PackedLJMixingData;
    size_t numRegisteredSiteTypes = 0;
    const PackedLJMixingData * __restrict computedLJMixingDataRow = nullptr;

    if constexpr (useMixing) {
      const auto &computedLJMixingData = _PPLibrary->getComputedLJMixingData();
      numRegisteredSiteTypes = _PPLibrary->getNumberRegisteredSiteTypes();
      computedLJMixingDataRow = computedLJMixingData.data() + typeIdI * numRegisteredSiteTypes;
    }

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    // Counters for when countFLOPs is activated
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    const auto threadNum = autopas::autopas_get_thread_num();

    // this is a magic number, that should correspond to at least
    // vectorization width*N have testet multiple sizes:
    // 4: does not give a speedup, slower than original AoSFunctor
    // 8: small speedup compared to AoS
    // 12: highest speedup compared to Aos
    // 16: smaller speedup
    // in theory this is a variable, we could auto-tune over...
#ifdef __AVX512F__
    // use a multiple of 8 for avx
    constexpr size_t vecSize = 16;
#else
    // for everything else 12 is faster
    constexpr size_t vecSize = 48;
#endif

    size_t jOff = 0;

    // if the size of the verlet list is larger than the given size vecsize,
    // we will use a vectorized version.
    if (neighborListSize >= vecSize) {
      alignas(64) std::array<SoAFloatPrecision, vecSize> xArr, yArr, zArr, forceXArr, forceYArr, forceZArr, ownedStateArr;

      [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecSize> sigmaSquareds;
      [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecSize> epsilon24s;
      [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecSize> shift6s;
      // loop over the verlet list from 0 to x*vecsize
      // in each iteration we calculate the interactions of particle i with
      // vecsize particles in the neighborlist of particle i starting at
      // particle joff
      #pragma code_align 64
      for (; jOff < neighborListSize - vecSize + 1; jOff += vecSize) {
        const auto *const __restrict neighborListPtrWithOffset = neighborListPtr + jOff;
        // gather position of particle j
        #pragma omp simd
        for (size_t tmpJ = 0; tmpJ < vecSize; tmpJ++) {
          const autopas::SoAIndexIntType indexInSoAJ = neighborListPtrWithOffset[tmpJ];
          xArr[tmpJ] = xPtr[indexInSoAJ];
          yArr[tmpJ] = yPtr[indexInSoAJ];
          zArr[tmpJ] = zPtr[indexInSoAJ];
          ownedStateArr[tmpJ] = ownedStatePtr[indexInSoAJ] == autopas::OwnershipState::dummy ? SoAFloatPrecision{0} : SoAFloatPrecision{1.};

          if constexpr (useMixing) {
            const auto typeIdJ = typeIdPtr[indexInSoAJ];
            const PackedLJMixingData &mixingDataJ = computedLJMixingDataRow[typeIdJ];
            sigmaSquareds[tmpJ] = mixingDataJ.sigmaSquared;
            epsilon24s[tmpJ] = mixingDataJ.epsilon24;
            if constexpr (applyShift) {
              shift6s[tmpJ] = mixingDataJ.shift6;
            }
          }
        }

        if constexpr (!calculateGlobals && !countFLOPs) {
        #pragma omp simd reduction(+ : forceXAcc, forceYAcc, forceZAcc) safelen(vecSize)
        #pragma code_align 64
          for (size_t j = 0; j < vecSize; j++) {
            if constexpr (useMixing) {
              sigmaSquared = sigmaSquareds[j];
              epsilon24 = epsilon24s[j];
              if constexpr (applyShift) {
                shift6 = shift6s[j];
              }
            }

            const SoAFloatPrecision drX = xI - xArr[j];
            const SoAFloatPrecision drY = yI - yArr[j];
            const SoAFloatPrecision drZ = zI - zArr[j];

            const SoAFloatPrecision drX2 = drX * drX;
            const SoAFloatPrecision drY2 = drY * drY;
            const SoAFloatPrecision drZ2 = drZ * drZ;

            const SoAFloatPrecision dr2 = drX2 + drY2 + drZ2;

            // Mask away if distance is too large or any particle is a dummy. ownedStateI was already checked
            // previously.
            const SoAFloatPrecision belowCutOff = dr2 <= cutoffSquared ? SoAFloatPrecision{1.} : SoAFloatPrecision{0};
            const SoAFloatPrecision mask = belowCutOff * ownedStateArr[j];

            const SoAFloatPrecision inverseDr2 = SoAFloatPrecision{1} / dr2;
            const SoAFloatPrecision lj2 = sigmaSquared * inverseDr2;
            const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
            const SoAFloatPrecision lj12 = lj6 * lj6;
            const SoAFloatPrecision lj12m6 = lj12 - lj6;
            const SoAFloatPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * inverseDr2;

            const SoAFloatPrecision forceXIJ = drX * fac;
            const SoAFloatPrecision forceYIJ = drY * fac;
            const SoAFloatPrecision forceZIJ = drZ * fac;

            forceXAcc += forceXIJ;
            forceYAcc += forceYIJ;
            forceZAcc += forceZIJ;

            if constexpr (newton3) {
              forceXArr[j] = forceXIJ;
              forceYArr[j] = forceYIJ;
              forceZArr[j] = forceZIJ;
            }
          }
        } else {
        #pragma omp simd reduction(+ : forceXAcc, forceYAcc, forceZAcc, potentialEnergySum, virialSumX, virialSumY,\
                               virialSumZ, numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum,\
                               numGlobalCalcsN3Sum, numGlobalCalcsNoN3Sum) safelen(vecSize)
#pragma code_align 64
          for (size_t j = 0; j < vecSize; j++) {
            if constexpr (useMixing) {
              sigmaSquared = sigmaSquareds[j];
              epsilon24 = epsilon24s[j];
              if constexpr (applyShift) {
                shift6 = shift6s[j];
              }
            }

            const SoAFloatPrecision drX = xI - xArr[j];
            const SoAFloatPrecision drY = yI - yArr[j];
            const SoAFloatPrecision drZ = zI - zArr[j];

            const SoAFloatPrecision drX2 = drX * drX;
            const SoAFloatPrecision drY2 = drY * drY;
            const SoAFloatPrecision drZ2 = drZ * drZ;

            const SoAFloatPrecision dr2 = drX2 + drY2 + drZ2;

            // Mask away if distance is too large or any particle is a dummy. ownedStateI was already checked
            // previously.
            const SoAFloatPrecision belowCutOff = dr2 <= cutoffSquared ? SoAFloatPrecision{1.} : SoAFloatPrecision{0};
            const SoAFloatPrecision mask = belowCutOff * ownedStateArr[j];

            const SoAFloatPrecision inverseDr2 = SoAFloatPrecision{1} / dr2;
            const SoAFloatPrecision lj2 = sigmaSquared * inverseDr2;
            const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
            const SoAFloatPrecision lj12 = lj6 * lj6;
            const SoAFloatPrecision lj12m6 = lj12 - lj6;
            const SoAFloatPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * inverseDr2;

            const SoAFloatPrecision forceXIJ = drX * fac;
            const SoAFloatPrecision forceYIJ = drY * fac;
            const SoAFloatPrecision forceZIJ = drZ * fac;

            forceXAcc += forceXIJ;
            forceYAcc += forceYIJ;
            forceZAcc += forceZIJ;

            if constexpr (newton3) {
              forceXArr[j] = forceXIJ;
              forceYArr[j] = forceYIJ;
              forceZArr[j] = forceZIJ;
            }

            if constexpr (countFLOPs) {
              numDistanceCalculationSum += ownedStateArr[j] != 0 ? 1 : 0;
              if constexpr (newton3) {
                numKernelCallsN3Sum += mask;
              } else {
                numKernelCallsNoN3Sum += mask;
              }
            }

            if constexpr (calculateGlobals) {
              SoAFloatPrecision virialx = drX * forceXIJ;
              SoAFloatPrecision virialy = drY * forceYIJ;
              SoAFloatPrecision virialz = drZ * forceZIJ;
              SoAFloatPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

              // We add 6 times the potential energy for each owned particle. Correction of total sum in endTraversal().
              const SoAFloatPrecision energyFactor =
                  (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
                  (newton3 ? (ownedStateArr[j] == 1. ? 1. : 0.) : 0.); // This is incorrect rn as 1 is set for both halo and owned but newton3 is not used here anyways... for I it is correct
              potentialEnergySum += potentialEnergy6 * energyFactor;
              virialSumX += virialx * energyFactor;
              virialSumY += virialy * energyFactor;
              virialSumZ += virialz * energyFactor;

              if constexpr (countFLOPs) {
                if constexpr (newton3) {
                  numGlobalCalcsN3Sum += mask;
                } else {
                  numGlobalCalcsNoN3Sum += mask;
                }
              }
            }
          }
        }

        // scatter the forces to where they belong, this is only needed for newton3
        if constexpr (newton3) {
        #pragma omp simd safelen(vecSize)
          for (size_t tmpj = 0; tmpj < vecSize; tmpj++) {
            const size_t j = neighborListPtr[jOff + tmpj];
            forceXPtr[j] -= forceXArr[tmpj];
            forceYPtr[j] -= forceYArr[tmpj];
            forceZPtr[j] -= forceZArr[tmpj];
          }
        }
      }
    }

    // this loop goes over the remainder and uses no optimizations
    for (size_t jNeighIndex = jOff; jNeighIndex < neighborListSize; ++jNeighIndex) {
      size_t j = neighborList[jNeighIndex];
      if (j == indexFirst) {
        continue;
      }

      const auto typeIdJ = typeIdPtr[j];

      if constexpr (useMixing) {
        sigmaSquared = _PPLibrary->getMixingSigmaSquared(typeIdI, typeIdJ);
        epsilon24 = _PPLibrary->getMixing24Epsilon(typeIdI, typeIdJ);
        if constexpr (applyShift) {
          shift6 = _PPLibrary->getMixingShift6(typeIdI, typeIdJ);
        }
      }

      const auto ownedStateJ = ownedStatePtr[j];
      if (ownedStateJ == autopas::OwnershipState::dummy) {
        continue;
      }

      const SoAFloatPrecision drX = xI - xPtr[j];
      const SoAFloatPrecision drY = yI - yPtr[j];
      const SoAFloatPrecision drZ = zI - zPtr[j];

      const SoAFloatPrecision drX2 = drX * drX;
      const SoAFloatPrecision drY2 = drY * drY;
      const SoAFloatPrecision drZ2 = drZ * drZ;

      const SoAFloatPrecision dr2 = drX2 + drY2 + drZ2;

      if constexpr (countFLOPs) {
        numDistanceCalculationSum += 1;
      }

      if (dr2 > cutoffSquared) {
        continue;
      }

      const SoAFloatPrecision inverseDr2 = 1. / dr2;
      const SoAFloatPrecision lj2 = sigmaSquared * inverseDr2;
      const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
      const SoAFloatPrecision lj12 = lj6 * lj6;
      const SoAFloatPrecision lj12m6 = lj12 - lj6;
      const SoAFloatPrecision fac = epsilon24 * (lj12 + lj12m6) * inverseDr2;

      const SoAFloatPrecision forceXIJ = drX * fac;
      const SoAFloatPrecision forceYIJ = drY * fac;
      const SoAFloatPrecision forceZIJ = drZ * fac;

      forceXAcc += forceXIJ;
      forceYAcc += forceYIJ;
      forceZAcc += forceZIJ;

      if constexpr (newton3) {
        forceXPtr[j] -= forceXIJ;
        forceYPtr[j] -= forceYIJ;
        forceZPtr[j] -= forceZIJ;
      }

      if constexpr (countFLOPs) {
        if constexpr (newton3) {
          numKernelCallsN3Sum += 1;
        } else {
          numKernelCallsNoN3Sum += 1;
        }
      }

      if constexpr (calculateGlobals) {
        SoAFloatPrecision virialx = drX * forceXIJ;
        SoAFloatPrecision virialy = drY * forceYIJ;
        SoAFloatPrecision virialz = drZ * forceZIJ;
        SoAFloatPrecision potentialEnergy6 = (epsilon24 * lj12m6 + shift6);

        // We add 6 times the potential energy for each owned particle. The total sum is corrected in endTraversal().
        const SoAFloatPrecision energyFactor =
            (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
            (newton3 ? (ownedStateJ == autopas::OwnershipState::owned ? 1. : 0.) : 0.); //once again incorrect check for J
        potentialEnergySum += potentialEnergy6 * energyFactor;
        virialSumX += virialx * energyFactor;
        virialSumY += virialy * energyFactor;
        virialSumZ += virialz * energyFactor;

        if constexpr (countFLOPs) {
          if constexpr (newton3) {
            ++numGlobalCalcsN3Sum;
          } else {
            ++numGlobalCalcsNoN3Sum;
          }
        }
      }
    }

    if (forceXAcc != 0 or forceYAcc != 0 or forceZAcc != 0) {
      forceXPtr[indexFirst] += forceXAcc;
      forceYPtr[indexFirst] += forceYAcc;
      forceZPtr[indexFirst] += forceZAcc;
    }

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadNum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadNum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadNum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadNum].numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
      _aosThreadDataFLOPs[threadNum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;
    }

    if constexpr (calculateGlobals) {
      _aosThreadDataGlobals[threadNum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadNum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadNum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadNum].virialSum[2] += virialSumZ;
    }
  }


  #pragma code_align (64)
  template <bool newton3>
  void SoAFunctorVerletOptimizedCompactAoSImpl(autopas::VerletListsLJCompactAoS<Particle_T> &compactAoS, const size_t indexFirst,
                            const std::vector<autopas::SoAIndexIntType, autopas::AlignedAllocator<autopas::SoAIndexIntType>> &neighborList) {

    auto *const __restrict forceXPtr = compactAoS._soa.template begin<Particle_T::AttributeNames::forceX>();
    auto *const __restrict forceYPtr = compactAoS._soa.template begin<Particle_T::AttributeNames::forceY>();
    auto *const __restrict forceZPtr = compactAoS._soa.template begin<Particle_T::AttributeNames::forceZ>();

    SoAFloatPrecision forceXAcc = 0;
    SoAFloatPrecision forceYAcc = 0;
    SoAFloatPrecision forceZAcc = 0;

    const auto *const __restrict compactDataPtr = compactAoS._compactParticles.data();

    const size_t neighborListSize = neighborList.size();
    const autopas::SoAIndexIntType *const __restrict neighborListPtr = neighborList.data();

    const auto particleI = compactDataPtr[indexFirst];


    const auto typeIdI = particleI.typeId;
    const auto ownedStateI = particleI.ownershipState;

    if (ownedStateI == 0) {
      return;
    }

    const SoAFloatPrecision xI = particleI.posX;
    const SoAFloatPrecision yI = particleI.posY;
    const SoAFloatPrecision zI = particleI.posZ;

    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    SoAFloatPrecision shift6 = _shift6;
    SoAFloatPrecision sigmaSquared = _sigmaSquared;
    SoAFloatPrecision epsilon24 = _epsilon24;

    using PackedLJMixingData = ParticlePropertiesLibrary<SoAFloatPrecision, unsigned long>::PackedLJMixingData;
    size_t numRegisteredSiteTypes = 0;
    const PackedLJMixingData * __restrict computedLJMixingDataRow = nullptr;

    if constexpr (useMixing) {
      const auto &computedLJMixingData = _PPLibrary->getComputedLJMixingData();
      numRegisteredSiteTypes = _PPLibrary->getNumberRegisteredSiteTypes();
      computedLJMixingDataRow = computedLJMixingData.data() + typeIdI * numRegisteredSiteTypes;
    }

    SoAFloatPrecision potentialEnergySum = 0.;
    SoAFloatPrecision virialSumX = 0.;
    SoAFloatPrecision virialSumY = 0.;
    SoAFloatPrecision virialSumZ = 0.;

    // Counters for when countFLOPs is activated
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    const auto threadNum = autopas::autopas_get_thread_num();

    // this is a magic number, that should correspond to at least
    // vectorization width*N have testet multiple sizes:
    // 4: does not give a speedup, slower than original AoSFunctor
    // 8: small speedup compared to AoS
    // 12: highest speedup compared to Aos
    // 16: smaller speedup
    // in theory this is a variable, we could auto-tune over...
#ifdef __AVX512F__
    // use a multiple of 8 for avx
    constexpr size_t vecSize = 16;
#else
    // for everything else 12 is faster
    constexpr size_t vecSize = 48;
#endif

    size_t jOff = 0;

    // if the size of the verlet list is larger than the given size vecsize,
    // we will use a vectorized version.
    if (neighborListSize >= vecSize) {
      alignas(64) std::array<SoAFloatPrecision, vecSize> xArr, yArr, zArr, forceXArr, forceYArr, forceZArr, ownedStateArr;

      [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecSize> sigmaSquareds;
      [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecSize> epsilon24s;
      [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<SoAFloatPrecision, vecSize> shift6s;
      // loop over the verlet list from 0 to x*vecsize
      // in each iteration we calculate the interactions of particle i with
      // vecsize particles in the neighborlist of particle i starting at
      // particle joff
      for (; jOff < neighborListSize - vecSize + 1; jOff += vecSize) {
        const auto *neighborListPtrWithOffset = neighborListPtr + jOff;
        // gather position of particle j
        #pragma omp simd safelen(vecSize)
        for (size_t tmpJ = 0; tmpJ < vecSize; tmpJ++) {
          const autopas::SoAIndexIntType indexInSoAJ = neighborListPtrWithOffset[tmpJ];
          const auto particleJ = compactDataPtr[indexInSoAJ];
          /*const auto ownedStateAndTypeIdJ = ownedStateAndTypeIdPtr[indexInSoAJ];*/
          /*xArr[tmpJ] = xPtr[indexInSoAJ];
          yArr[tmpJ] = yPtr[indexInSoAJ];
          zArr[tmpJ] = zPtr[indexInSoAJ];*/
          xArr[tmpJ] = particleJ.posX;
          yArr[tmpJ] = particleJ.posY;
          zArr[tmpJ] = particleJ.posZ;
          ownedStateArr[tmpJ] = particleJ.ownershipState == 0 ? SoAFloatPrecision{0} : SoAFloatPrecision{1.};

          if constexpr (useMixing) {
            const size_t typeIdJ = particleJ.typeId;
            const PackedLJMixingData &mixingDataJ = computedLJMixingDataRow[typeIdJ];
            sigmaSquareds[tmpJ] = mixingDataJ.sigmaSquared;
            epsilon24s[tmpJ] = mixingDataJ.epsilon24;
            if constexpr (applyShift) {
              shift6s[tmpJ] = mixingDataJ.shift6;
            }
          }
        }

        if constexpr (!calculateGlobals && !countFLOPs) {
        #pragma omp simd reduction(+ : forceXAcc, forceYAcc, forceZAcc) safelen(vecSize)
        #pragma code_align 64
          for (size_t j = 0; j < vecSize; j++) {
            if constexpr (useMixing) {
              sigmaSquared = sigmaSquareds[j];
              epsilon24 = epsilon24s[j];
              if constexpr (applyShift) {
                shift6 = shift6s[j];
              }
            }

            const SoAFloatPrecision drX = xI - xArr[j];
            const SoAFloatPrecision drY = yI - yArr[j];
            const SoAFloatPrecision drZ = zI - zArr[j];

            const SoAFloatPrecision drX2 = drX * drX;
            const SoAFloatPrecision drY2 = drY * drY;
            const SoAFloatPrecision drZ2 = drZ * drZ;

            const SoAFloatPrecision dr2 = drX2 + drY2 + drZ2;

            // Mask away if distance is too large or any particle is a dummy. ownedStateI was already checked
            // previously.
            const SoAFloatPrecision belowCutOff = dr2 <= cutoffSquared ? SoAFloatPrecision{1.} : SoAFloatPrecision{0};
            const SoAFloatPrecision mask = belowCutOff * ownedStateArr[j];

            const SoAFloatPrecision inverseDr2 = 1. / dr2;
            const SoAFloatPrecision lj2 = sigmaSquared * inverseDr2;
            const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
            const SoAFloatPrecision lj12 = lj6 * lj6;
            const SoAFloatPrecision lj12m6 = lj12 - lj6;
            const SoAFloatPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * inverseDr2;

            const SoAFloatPrecision forceXIJ = drX * fac;
            const SoAFloatPrecision forceYIJ = drY * fac;
            const SoAFloatPrecision forceZIJ = drZ * fac;

            forceXAcc += forceXIJ;
            forceYAcc += forceYIJ;
            forceZAcc += forceZIJ;

            if constexpr (newton3) {
              forceXArr[j] = forceXIJ;
              forceYArr[j] = forceYIJ;
              forceZArr[j] = forceZIJ;
            }
          }
        } else {
        #pragma omp simd reduction(+ : forceXAcc, forceYAcc, forceZAcc, potentialEnergySum, virialSumX, virialSumY,\
                               virialSumZ, numDistanceCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum,\
                               numGlobalCalcsN3Sum, numGlobalCalcsNoN3Sum) safelen(vecSize)
          for (size_t j = 0; j < vecSize; j++) {
            if constexpr (useMixing) {
              sigmaSquared = sigmaSquareds[j];
              epsilon24 = epsilon24s[j];
              if constexpr (applyShift) {
                shift6 = shift6s[j];
              }
            }

            const auto ownedStateJ = ownedStateArr[j];

            const SoAFloatPrecision drX = xI - xArr[j];
            const SoAFloatPrecision drY = yI - yArr[j];
            const SoAFloatPrecision drZ = zI - zArr[j];

            const SoAFloatPrecision drX2 = drX * drX;
            const SoAFloatPrecision drY2 = drY * drY;
            const SoAFloatPrecision drZ2 = drZ * drZ;

            const SoAFloatPrecision dr2 = drX2 + drY2 + drZ2;

            // Mask away if distance is too large or any particle is a dummy. ownedStateI was already checked
            // previously.
            const bool mask = dr2 <= cutoffSquared and ownedStateJ != 0;

            const SoAFloatPrecision inverseDr2 = 1. / dr2;
            const SoAFloatPrecision lj2 = sigmaSquared * inverseDr2;
            const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
            const SoAFloatPrecision lj12 = lj6 * lj6;
            const SoAFloatPrecision lj12m6 = lj12 - lj6;
            const SoAFloatPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * inverseDr2;

            const SoAFloatPrecision forceXIJ = drX * fac;
            const SoAFloatPrecision forceYIJ = drY * fac;
            const SoAFloatPrecision forceZIJ = drZ * fac;

            forceXAcc += forceXIJ;
            forceYAcc += forceYIJ;
            forceZAcc += forceZIJ;

            if constexpr (newton3) {
              forceXArr[j] = forceXIJ;
              forceYArr[j] = forceYIJ;
              forceZArr[j] = forceZIJ;
            }

            if constexpr (countFLOPs) {
              numDistanceCalculationSum += ownedStateJ != 0 ? 1 : 0;
              if constexpr (newton3) {
                numKernelCallsN3Sum += mask;
              } else {
                numKernelCallsNoN3Sum += mask;
              }
            }

            if constexpr (calculateGlobals) {
              SoAFloatPrecision virialx = drX * forceXIJ;
              SoAFloatPrecision virialy = drY * forceYIJ;
              SoAFloatPrecision virialz = drZ * forceZIJ;
              SoAFloatPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

              // We add 6 times the potential energy for each owned particle. Correction of total sum in endTraversal().
              const SoAFloatPrecision energyFactor =
                  (ownedStateI == 1 ? 1. : 0.) +
                  (newton3 ? (ownedStateJ == 1 ? 1. : 0.) : 0.);
              potentialEnergySum += potentialEnergy6 * energyFactor;
              virialSumX += virialx * energyFactor;
              virialSumY += virialy * energyFactor;
              virialSumZ += virialz * energyFactor;

              if constexpr (countFLOPs) {
                if constexpr (newton3) {
                  numGlobalCalcsN3Sum += mask;
                } else {
                  numGlobalCalcsNoN3Sum += mask;
                }
              }
            }
          }
        }

        // scatter the forces to where they belong, this is only needed for newton3
        if constexpr (newton3) {
        #pragma omp simd safelen(vecSize)
          for (size_t tmpj = 0; tmpj < vecSize; tmpj++) {
            const size_t j = neighborListPtr[jOff + tmpj];
            forceXPtr[j] -= forceXArr[tmpj];
            forceYPtr[j] -= forceYArr[tmpj];
            forceZPtr[j] -= forceZArr[tmpj];
          }
        }
      }
    }

    // this loop goes over the remainder and uses no optimizations
    for (size_t jNeighIndex = jOff; jNeighIndex < neighborListSize; ++jNeighIndex) {
      size_t j = neighborList[jNeighIndex];
      if (j == indexFirst) {
        continue;
      }

      if constexpr (useMixing) {
        sigmaSquared = _PPLibrary->getMixingSigmaSquared(compactDataPtr[indexFirst].typeId, compactDataPtr[j].typeId);
        epsilon24 = _PPLibrary->getMixing24Epsilon(compactDataPtr[indexFirst].typeId, compactDataPtr[j].typeId);
        if constexpr (applyShift) {
          shift6 = _PPLibrary->getMixingShift6(compactDataPtr[indexFirst].typeId, compactDataPtr[j].typeId);
        }
      }

      const auto ownedStateJ = compactDataPtr[j].ownershipState;
      if (ownedStateJ == 0) {
        continue;
      }

      const SoAFloatPrecision drX = xI - compactDataPtr[j].posX;
      const SoAFloatPrecision drY = yI - compactDataPtr[j].posY;
      const SoAFloatPrecision drZ = zI - compactDataPtr[j].posZ;

      const SoAFloatPrecision drX2 = drX * drX;
      const SoAFloatPrecision drY2 = drY * drY;
      const SoAFloatPrecision drZ2 = drZ * drZ;

      const SoAFloatPrecision dr2 = drX2 + drY2 + drZ2;

      if constexpr (countFLOPs) {
        numDistanceCalculationSum += 1;
      }

      if (dr2 > cutoffSquared) {
        continue;
      }

      const SoAFloatPrecision inverseDr2 = 1. / dr2;
      const SoAFloatPrecision lj2 = sigmaSquared * inverseDr2;
      const SoAFloatPrecision lj6 = lj2 * lj2 * lj2;
      const SoAFloatPrecision lj12 = lj6 * lj6;
      const SoAFloatPrecision lj12m6 = lj12 - lj6;
      const SoAFloatPrecision fac = epsilon24 * (lj12 + lj12m6) * inverseDr2;

      const SoAFloatPrecision forceXIJ = drX * fac;
      const SoAFloatPrecision forceYIJ = drY * fac;
      const SoAFloatPrecision forceZIJ = drZ * fac;

      forceXAcc += forceXIJ;
      forceYAcc += forceYIJ;
      forceZAcc += forceZIJ;

      if constexpr (newton3) {
        forceXPtr[j] -= forceXIJ;
        forceYPtr[j] -= forceYIJ;
        forceZPtr[j] -= forceZIJ;
      }

      if constexpr (countFLOPs) {
        if constexpr (newton3) {
          numKernelCallsN3Sum += 1;
        } else {
          numKernelCallsNoN3Sum += 1;
        }
      }

      if constexpr (calculateGlobals) {
        SoAFloatPrecision virialx = drX * forceXIJ;
        SoAFloatPrecision virialy = drY * forceYIJ;
        SoAFloatPrecision virialz = drZ * forceZIJ;
        SoAFloatPrecision potentialEnergy6 = (epsilon24 * lj12m6 + shift6);

        // We add 6 times the potential energy for each owned particle. The total sum is corrected in endTraversal().
        const SoAFloatPrecision energyFactor =
            (ownedStateI == 1 ? 1. : 0.) +
            (newton3 ? (ownedStateJ == 1 ? 1. : 0.) : 0.);
        potentialEnergySum += potentialEnergy6 * energyFactor;
        virialSumX += virialx * energyFactor;
        virialSumY += virialy * energyFactor;
        virialSumZ += virialz * energyFactor;

        if constexpr (countFLOPs) {
          if constexpr (newton3) {
            ++numGlobalCalcsN3Sum;
          } else {
            ++numGlobalCalcsNoN3Sum;
          }
        }
      }
    }

    if (forceXAcc != 0 or forceYAcc != 0 or forceZAcc != 0) {
      forceXPtr[indexFirst] += forceXAcc;
      forceYPtr[indexFirst] += forceYAcc;
      forceZPtr[indexFirst] += forceZAcc;
    }

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadNum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadNum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadNum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadNum].numGlobalCalcsNoN3 += numGlobalCalcsNoN3Sum;
      _aosThreadDataFLOPs[threadNum].numGlobalCalcsN3 += numGlobalCalcsN3Sum;
    }

    if (calculateGlobals) {
      _aosThreadDataGlobals[threadNum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadNum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadNum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadNum].virialSum[2] += virialSumZ;
    }
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
  // not const because they might be reset through PPL
  double _epsilon24, _sigmaSquared, _shift6 = 0;

  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

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
