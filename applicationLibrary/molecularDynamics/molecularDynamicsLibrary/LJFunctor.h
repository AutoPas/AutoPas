/**
 * @file LJFunctor.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include "ParticlePropertiesLibrary.h"
#include "autopas/pairwiseFunctors/Functor.h"
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
 * @tparam Particle The type of particle.
 * @tparam applyShift Switch for the lj potential to be truncated shifted.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool applyShift = false, bool useMixing = false,
          autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both, bool calculateGlobals = false,
          bool countFLOPs = false, bool relevantForTuning = true>
class LJFunctor : public autopas::Functor<Particle, LJFunctor<Particle, applyShift, useMixing, useNewton3,
                                                              calculateGlobals, countFLOPs, relevantForTuning>> {
  /**
   * FloatType used for calculations
   */
  using CalcPrecision = typename Particle::ParticleCalcPrecision;

  /**
   * FloatType used for accumulations or more relevant calculations
   */
  using AccuPrecision = typename Particle::ParticleAccuPrecision;

  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Precision of SoA entries.
   * TODO MP Remove
   */
  using SoAFloatPrecision = typename AccuPrecision;

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
  explicit LJFunctor(CalcPrecision cutoff, void * /*dummy*/)
      : autopas::Functor<Particle, LJFunctor<Particle, applyShift, useMixing, useNewton3, calculateGlobals, countFLOPs,
                                             relevantForTuning>>(cutoff),
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
  explicit LJFunctor(CalcPrecision cutoff) : LJFunctor(cutoff, nullptr) {
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
  explicit LJFunctor(CalcPrecision cutoff, ParticlePropertiesLibrary<CalcPrecision, size_t> &particlePropertiesLibrary)
      : LJFunctor(cutoff, nullptr) {
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

  void AoSFunctor(Particle &i, Particle &j, bool newton3) final {
    using namespace autopas::utils::ArrayMath::literals;

    if (i.isDummy() or j.isDummy()) {
      return;
    }

    const auto threadnum = autopas::autopas_get_thread_num();

    if constexpr (countFLOPs) {
      ++_aosThreadDataFLOPs[threadnum].numDistCalls;
    }

    CalcPrecision sigmaSquared = _sigmaSquared;
    CalcPrecision epsilon24 = _epsilon24;
    CalcPrecision shift6 = _shift6;
    if constexpr (useMixing) {
      sigmaSquared = _PPLibrary->getMixingSigmaSquared(i.getTypeId(), j.getTypeId());
      epsilon24 = _PPLibrary->getMixing24Epsilon(i.getTypeId(), j.getTypeId());
      if constexpr (applyShift) {
        shift6 = _PPLibrary->getMixingShift6(i.getTypeId(), j.getTypeId());
      }
    }
    std::array<CalcPrecision, 3> dr = i.getR() - j.getR();
    CalcPrecision dr2 = autopas::utils::ArrayMath::dot(dr, dr);

    if (dr2 > _cutoffSquared) {
      return;
    }

    CalcPrecision invdr2 = static_cast<CalcPrecision>(1.) / dr2;
    CalcPrecision lj6 = sigmaSquared * invdr2;
    lj6 = lj6 * lj6 * lj6;
    CalcPrecision lj12 = lj6 * lj6;
    CalcPrecision lj12m6 = lj12 - lj6;
    CalcPrecision fac = epsilon24 * (lj12 + lj12m6) * invdr2;
    std::array<CalcPrecision, 3> f = dr * fac;
    std::array<CalcPrecision, 3> convertedF;
    if constexpr (std::is_same_v<CalcPrecision, AccuPrecision>) {
      convertedF = f;
    } else {
      convertedF = std::copy(f.begin(), f.end(), convertedF.begin());
    }
    i.addF(convertedF);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(convertedF);
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
      AccuPrecision potentialEnergy6 = static_cast<AccuPrecision>(epsilon24) * static_cast<AccuPrecision>(lj12m6) +
                                       static_cast<AccuPrecision>(shift6);

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
   * @copydoc autopas::Functor::SoAFunctorSingle()
   * This functor will always use a newton3 like traversal of the soa.
   * However, it still needs to know about newton3 to correctly add up the global values.
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
    if (soa.size() == 0) return;

    const auto threadnum = autopas::autopas_get_thread_num();

    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    AccuPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    AccuPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    AccuPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();
    // the local redeclaration of the following values helps the SoAFloatPrecision-generation of various compilers.
    const CalcPrecision cutoffSquared = _cutoffSquared;

    AccuPrecision potentialEnergySum = 0.;  // Note: This is not the potential energy but some fixed multiple of it.
    AccuPrecision virialSumX = 0.;
    AccuPrecision virialSumY = 0.;
    AccuPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsSum = 0;

    std::vector<CalcPrecision, autopas::AlignedAllocator<CalcPrecision>> sigmaSquareds;
    std::vector<CalcPrecision, autopas::AlignedAllocator<CalcPrecision>> epsilon24s;
    std::vector<CalcPrecision, autopas::AlignedAllocator<CalcPrecision>> shift6s;
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

    const CalcPrecision const_shift6 = _shift6;
    const CalcPrecision const_sigmaSquared = _sigmaSquared;
    const CalcPrecision const_epsilon24 = _epsilon24;

    for (unsigned int i = 0; i < soa.size(); ++i) {
      const auto ownedStateI = ownedStatePtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      AccuPrecision fxacc = 0.;
      AccuPrecision fyacc = 0.;
      AccuPrecision fzacc = 0.;

      if constexpr (useMixing) {
        for (unsigned int j = 0; j < soa.size(); ++j) {
          auto mixingData = _PPLibrary->getMixingData(typeptr[i], typeptr[j]);
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
        CalcPrecision shift6 = const_shift6;
        CalcPrecision sigmaSquared = const_sigmaSquared;
        CalcPrecision epsilon24 = const_epsilon24;
        if constexpr (useMixing) {
          sigmaSquared = sigmaSquareds[j];
          epsilon24 = epsilon24s[j];
          if constexpr (applyShift) {
            shift6 = shift6s[j];
          }
        }

        const auto ownedStateJ = ownedStatePtr[j];

        const CalcPrecision drx = xptr[i] - xptr[j];
        const CalcPrecision dry = yptr[i] - yptr[j];
        const CalcPrecision drz = zptr[i] - zptr[j];

        const CalcPrecision drx2 = drx * drx;
        const CalcPrecision dry2 = dry * dry;
        const CalcPrecision drz2 = drz * drz;

        const CalcPrecision dr2 = drx2 + dry2 + drz2;

        // Mask away if distance is too large or any particle is a dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

        const CalcPrecision invdr2 = static_cast<CalcPrecision>(1.) / dr2;
        const CalcPrecision lj2 = sigmaSquared * invdr2;
        const CalcPrecision lj6 = lj2 * lj2 * lj2;
        const CalcPrecision lj12 = lj6 * lj6;
        const CalcPrecision lj12m6 = lj12 - lj6;
        const CalcPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * invdr2;

        const CalcPrecision fx = drx * fac;
        const CalcPrecision fy = dry * fac;
        const CalcPrecision fz = drz * fac;

        // implicit cast to (potentially) different precision
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
          const CalcPrecision virialx = drx * fx;
          const CalcPrecision virialy = dry * fy;
          const CalcPrecision virialz = drz * fz;
          const CalcPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

          // We add 6 times the potential energy for each owned particle. The total sum is corrected in endTraversal().
          CalcPrecision energyFactor = (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
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
      const int threadnum = autopas::autopas_get_thread_num();

      _aosThreadDataGlobals[threadnum].potentialEnergySum += potentialEnergySum;
      _aosThreadDataGlobals[threadnum].virialSum[0] += virialSumX;
      _aosThreadDataGlobals[threadnum].virialSum[1] += virialSumY;
      _aosThreadDataGlobals[threadnum].virialSum[2] += virialSumZ;
    }
  }

  /**
   * @copydoc autopas::Functor::SoAFunctorPair()
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
    [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle::AttributeNames::typeId>();

    // Checks whether the cells are halo cells.
    AccuPrecision potentialEnergySum = 0.;
    AccuPrecision virialSumX = 0.;
    AccuPrecision virialSumY = 0.;
    AccuPrecision virialSumZ = 0.;

    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    const CalcPrecision cutoffSquared = _cutoffSquared;
    CalcPrecision shift6 = _shift6;
    CalcPrecision sigmaSquared = _sigmaSquared;
    CalcPrecision epsilon24 = _epsilon24;

    // preload all sigma and epsilons for next vectorized region
    std::vector<CalcPrecision, autopas::AlignedAllocator<CalcPrecision>> sigmaSquareds;
    std::vector<CalcPrecision, autopas::AlignedAllocator<CalcPrecision>> epsilon24s;
    std::vector<CalcPrecision, autopas::AlignedAllocator<CalcPrecision>> shift6s;
    if constexpr (useMixing) {
      sigmaSquareds.resize(soa2.size());
      epsilon24s.resize(soa2.size());
      // if no mixing or mixing but no shift shift6 is constant therefore we do not need this vector.
      if constexpr (applyShift) {
        shift6s.resize(soa2.size());
      }
    }

    for (unsigned int i = 0; i < soa1.size(); ++i) {
      AccuPrecision fxacc = 0;
      AccuPrecision fyacc = 0;
      AccuPrecision fzacc = 0;

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

        const CalcPrecision drx = x1ptr[i] - x2ptr[j];
        const CalcPrecision dry = y1ptr[i] - y2ptr[j];
        const CalcPrecision drz = z1ptr[i] - z2ptr[j];

        const CalcPrecision drx2 = drx * drx;
        const CalcPrecision dry2 = dry * dry;
        const CalcPrecision drz2 = drz * drz;

        const CalcPrecision dr2 = drx2 + dry2 + drz2;

        // Mask away if distance is too large or any particle is a dummy.
        // Particle ownedStateI was already checked previously.
        const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

        const CalcPrecision invdr2 = 1. / dr2;
        const CalcPrecision lj2 = sigmaSquared * invdr2;
        const CalcPrecision lj6 = lj2 * lj2 * lj2;
        const CalcPrecision lj12 = lj6 * lj6;
        const CalcPrecision lj12m6 = lj12 - lj6;
        const CalcPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * invdr2;

        const CalcPrecision fx = drx * fac;
        const CalcPrecision fy = dry * fac;
        const CalcPrecision fz = drz * fac;

        // implicit cast to (potentially) different precision
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
          CalcPrecision virialx = drx * fx;
          CalcPrecision virialy = dry * fy;
          CalcPrecision virialz = drz * fz;
          CalcPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

          // We add 6 times the potential energy for each owned particle. The total sum is corrected in endTraversal().
          const CalcPrecision energyFactor = (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
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
   * @copydoc autopas::Functor::SoAFunctorVerlet()
   * @note If you want to parallelize this by openmp, please ensure that there
   * are no dependencies, i.e. introduce colors!
   */
  // clang-format on
  void SoAFunctorVerlet(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                        const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList,
                        bool newton3) final {
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

      AutoPasLog(TRACE, "Final potential energy {}", _potentialEnergySum);
      AutoPasLog(TRACE, "Final virial           {}", _virialSum[0] + _virialSum[1] + _virialSum[2]);
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
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();
    [[maybe_unused]] auto *const __restrict typeptr1 = soa.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa.template begin<Particle::AttributeNames::typeId>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    const CalcPrecision cutoffSquared = _cutoffSquared;
    CalcPrecision shift6 = _shift6;
    CalcPrecision sigmaSquared = _sigmaSquared;
    CalcPrecision epsilon24 = _epsilon24;

    AccuPrecision potentialEnergySum = 0.;
    AccuPrecision virialSumX = 0.;
    AccuPrecision virialSumY = 0.;
    AccuPrecision virialSumZ = 0.;

    // Counters for when countFLOPs is activated
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numGlobalCalcsN3Sum = 0;
    size_t numGlobalCalcsNoN3Sum = 0;

    AccuPrecision fxacc = 0;
    AccuPrecision fyacc = 0;
    AccuPrecision fzacc = 0;
    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict neighborListPtr = neighborList.data();

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
      alignas(64) std::array<CalcPrecision, vecsize> xtmp, ytmp, ztmp, xArr, yArr, zArr, fxArr, fyArr, fzArr;
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

        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<CalcPrecision, vecsize> sigmaSquareds;
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<CalcPrecision, vecsize> epsilon24s;
        [[maybe_unused]] alignas(autopas::DEFAULT_CACHE_LINE_SIZE) std::array<CalcPrecision, vecsize> shift6s;
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

          const CalcPrecision drx = xtmp[j] - xArr[j];
          const CalcPrecision dry = ytmp[j] - yArr[j];
          const CalcPrecision drz = ztmp[j] - zArr[j];

          const CalcPrecision drx2 = drx * drx;
          const CalcPrecision dry2 = dry * dry;
          const CalcPrecision drz2 = drz * drz;

          const CalcPrecision dr2 = drx2 + dry2 + drz2;

          // Mask away if distance is too large or any particle is a dummy.
          // Particle ownedStateI was already checked previously.
          const bool mask = dr2 <= cutoffSquared and ownedStateJ != autopas::OwnershipState::dummy;

          const CalcPrecision invdr2 = 1. / dr2;
          const CalcPrecision lj2 = sigmaSquared * invdr2;
          const CalcPrecision lj6 = lj2 * lj2 * lj2;
          const CalcPrecision lj12 = lj6 * lj6;
          const CalcPrecision lj12m6 = lj12 - lj6;
          const CalcPrecision fac = mask * epsilon24 * (lj12 + lj12m6) * invdr2;

          const CalcPrecision fx = drx * fac;
          const CalcPrecision fy = dry * fac;
          const CalcPrecision fz = drz * fac;

          // implicit cast to (potentially) different precision
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
            CalcPrecision virialx = drx * fx;
            CalcPrecision virialy = dry * fy;
            CalcPrecision virialz = drz * fz;
            CalcPrecision potentialEnergy6 = mask * (epsilon24 * lj12m6 + shift6);

            // We add 6 times the potential energy for each owned particle. The total sum is corrected in
            // endTraversal().
            const CalcPrecision energyFactor =
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

      const CalcPrecision drx = xptr[indexFirst] - xptr[j];
      const CalcPrecision dry = yptr[indexFirst] - yptr[j];
      const CalcPrecision drz = zptr[indexFirst] - zptr[j];

      const CalcPrecision drx2 = drx * drx;
      const CalcPrecision dry2 = dry * dry;
      const CalcPrecision drz2 = drz * drz;

      const CalcPrecision dr2 = drx2 + dry2 + drz2;

      if constexpr (countFLOPs) {
        numDistanceCalculationSum += 1;
      }

      if (dr2 > cutoffSquared) {
        continue;
      }

      const CalcPrecision invdr2 = 1. / dr2;
      const CalcPrecision lj2 = sigmaSquared * invdr2;
      const CalcPrecision lj6 = lj2 * lj2 * lj2;
      const CalcPrecision lj12 = lj6 * lj6;
      const CalcPrecision lj12m6 = lj12 - lj6;
      const CalcPrecision fac = epsilon24 * (lj12 + lj12m6) * invdr2;

      const CalcPrecision fx = drx * fac;
      const CalcPrecision fy = dry * fac;
      const CalcPrecision fz = drz * fac;

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
        CalcPrecision virialx = drx * fx;
        CalcPrecision virialy = dry * fy;
        CalcPrecision virialz = drz * fz;
        CalcPrecision potentialEnergy6 = (epsilon24 * lj12m6 + shift6);

        // We add 6 times the potential energy for each owned particle. The total sum is corrected in endTraversal().
        const CalcPrecision energyFactor = (ownedStateI == autopas::OwnershipState::owned ? 1. : 0.) +
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
    std::array<AccuPrecision, 3> virialSum;
    AccuPrecision potentialEnergySum;

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 4 * sizeof(AccuPrecision)) / sizeof(double)];
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

  const CalcPrecision _cutoffSquared;
  // not const because they might be reset through PPL
  CalcPrecision _epsilon24, _sigmaSquared, _shift6 = 0;

  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

  // sum of the potential energy, only calculated if calculateGlobals is true
  AccuPrecision _potentialEnergySum;

  // sum of the virial, only calculated if calculateGlobals is true
  std::array<AccuPrecision, 3> _virialSum;

  // thread buffer for aos
  std::vector<AoSThreadDataGlobals> _aosThreadDataGlobals{};
  std::vector<AoSThreadDataFLOPs> _aosThreadDataFLOPs{};

  // defines whether or whether not the global values are already preprocessed
  bool _postProcessed;
};
}  // namespace mdLib
