/**
 * @file DEMFunctor.h
 * @author J. Kim
 * @date 04/11/24
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
 * A functor to handle the interactions between particles using Discrete Element Method (DEM).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle The type of particle.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon and _sigma need to be set and the constructor with PPL can be omitted.
 * @tparam useNewton3 Switch for the functor to support newton3 on, off or both. See FunctorN3Modes for possible values.
 * @tparam calculateGlobals Defines whether the global values are to be calculated (energy, virial).
 * @tparam countFLOPs counts FLOPs and hitrate
 * @tparam relevantForTuning Whether or not the auto-tuner should consider this functor.
 */
template <class Particle, bool useMixing = false, autopas::FunctorN3Modes useNewton3 = autopas::FunctorN3Modes::Both,
          bool calculateGloabls = false, bool countFLOPs = false, bool relevantForTuning = true>
class DEMFunctor
    : public autopas::Functor<
          Particle, DEMFunctor<Particle, useMixing, useNewton3, calculateGloabls, countFLOPs, relevantForTuning>> {
  /**
   * Structure of the SoAs defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Precision of SoA entries.
   */
  using SoAFloatPrecision = typename Particle::ParticleSoAFloatPrecision;

 public:
  /**
   * Deleted default constructor
   */
  DEMFunctor() = delete;

 private:
  /**
   * Internal, actual constructor.
   * @param cutoff
   * @param springStiffness
   * @param viscousDamping
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit DEMFunctor(double cutoff, double elasticStiffness, double adhesiveStiffness, double frictionStiffness,
                      double normalViscosity, double frictionViscosity, double staticFrictionCoeff,
                      double dynamicFrictionCoeff, void * /*dummy*/)
      : autopas::Functor<Particle,
                         DEMFunctor<Particle, useMixing, useNewton3, calculateGloabls, countFLOPs, relevantForTuning>>(
            cutoff),
        _radius{cutoff / 5.},  // default initialization
        _cutoffSquared{cutoff * cutoff},
        _elasticStiffness{elasticStiffness},
        _adhesiveStiffness{adhesiveStiffness},
        _frictionStiffness{frictionStiffness},
        _normalViscosity{normalViscosity},
        _frictionViscosity{frictionViscosity},
        _staticFrictionCoeff{staticFrictionCoeff},
        _dynamicFrictionCoeff{dynamicFrictionCoeff},
        _potentialEnergySum{0.},
        _virialSum{0., 0., 0.},
        _postProcessed{false} {
    if constexpr (calculateGloabls) {
      _aosThreadDataGlobals.resize(autopas::autopas_get_max_threads());
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs.resize(autopas::autopas_get_max_threads());
    }
  }

 public:
  /**
   * Minimal Constructor for Functor with mixing disabled. Other Parameters are set with default values.
   * Values taken from https://www2.msm.ctw.utwente.nl/sluding/PAPERS/Alert_Luding2008.pdf page 15.
   *
   * @note Only to be used with mixing == false;
   */
  explicit DEMFunctor(double cutoff) : DEMFunctor(cutoff, 5., 2.5, 1., 5e-5, 1e-5, 1., 1., nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Minimal Constructor for Functor with mixing disabled. Other Parameters are set with default values.
   * Values taken from https://www2.msm.ctw.utwente.nl/sluding/PAPERS/Alert_Luding2008.pdf page 15.
   *
   * @note Only to be used with mixing == true.
   */
  explicit DEMFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : DEMFunctor(cutoff, 5., 2.5, 1., 5e-5, 1e-5, 1., 1., nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  /**
   * Constructor for Functor with mixing disabled. When using this functor it is necessary to call
   * setParticleProperties() to set internal constants because it does not use a particle properties library.
   *
   * @note Only to be used with mixing == false.
   *
   * @param cutoff
   */
  explicit DEMFunctor(double cutoff, double elasticStiffness, double adhesiveStiffness, double frictionStiffness,
                      double normalViscosity, double frictionViscosity, double staticFrictionCoeff,
                      double dynamicFrictionCoeff)
      : DEMFunctor(cutoff, elasticStiffness, adhesiveStiffness, frictionStiffness, normalViscosity, frictionViscosity,
                   staticFrictionCoeff, dynamicFrictionCoeff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like radius.
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit DEMFunctor(double cutoff, double elasticStiffness, double adhesiveStiffness, double frictionStiffness,
                      double normalViscosity, double frictionViscosity, double staticFrictionCoeff,
                      double dynamicFrictionCoeff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : DEMFunctor(cutoff, elasticStiffness, adhesiveStiffness, frictionStiffness, normalViscosity, frictionViscosity,
                   staticFrictionCoeff, dynamicFrictionCoeff, nullptr) {
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

    // Compute necessary values and check for distance and contact
    const double radius = _radius;
    const std::array<double, 3> displacement = i.getR() - j.getR();
    const double distSquared = autopas::utils::ArrayMath::dot(displacement, displacement);
    if (distSquared > _cutoffSquared) return;

    const double dist = autopas::utils::ArrayMath::L2Norm(displacement);
    const std::array<double, 3> normalUnit = displacement / dist;
    double overlap = 2. * radius - dist;
    if (useMixing) {
     const double radius_i = _PPLibrary->getRadius(i.getTypeId());
     const double radius_j = _PPLibrary->getRadius(j.getTypeId());
    overlap = radius_i + radius_j - dist;
    }
    if (overlap <= 0) return;

    const std::array<double, 3> relVel = i.getV() - j.getV();
    const double normalRelVelMag = autopas::utils::ArrayMath::dot(normalUnit, relVel);
    const std::array<double, 3> normalRelVel = autopas::utils::ArrayMath::mulScalar(normalUnit, normalRelVelMag);
    const std::array<double, 3> tanRelVel = relVel - normalRelVel;  // TODO: add angle velocities

    // Compute normal force
    const double normalContactFMag = _elasticStiffness * overlap - _normalViscosity * normalRelVelMag;
    const std::array<double, 3> normalContactF = autopas::utils::ArrayMath::mulScalar(normalUnit, normalContactFMag);

    const std::array<double, 3> normalLongRangeF = {0., 0., 0.};
    const std::array<double, 3> normalF = normalContactF + normalLongRangeF;

    // Compute tangential force
    const double coulombLimit = _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);
    std::array<double, 3> tanF =
        autopas::utils::ArrayMath::mulScalar(tanRelVel, -1. * _frictionViscosity);  // TODO: add tangential spring
    const double tanFMag = autopas::utils::ArrayMath::L2Norm(tanF);
    if (tanFMag > coulombLimit) {
      const double scale = _dynamicFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap) / tanFMag;
      tanF = autopas::utils::ArrayMath::mulScalar(tanF, scale);
    }

    // Compute total force
    auto totalF = normalF + tanF;

    // Apply forces
    i.addF(totalF);
    if (newton3) {
      // only if we use newton 3 here, we want to
      j.subF(totalF);
    }

    if constexpr (countFLOPs) {
      if (newton3) {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsN3;
      } else {
        ++_aosThreadDataFLOPs[threadnum].numKernelCallsNoN3;
      }
    }

    if constexpr (calculateGloabls) {
      // TODO

      // for non-newton3 the second particle will be considered in a separate calculation
      if (newton3 and j.isOwned()) {
        // TODO
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
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3)
      final {  // TODo: calculate the global values (potenetialEnergy, virial), add cases for useMixing
    if (soa.size() == 0) return;

    const auto threadnum = autopas::autopas_get_thread_num();

    // Initialize pointers to SoA Data
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    const auto *const __restrict vxptr = soa.template begin<Particle::AttributeNames::velocityX>();
    const auto *const __restrict vyptr = soa.template begin<Particle::AttributeNames::velocityY>();
    const auto *const __restrict vzptr = soa.template begin<Particle::AttributeNames::velocityZ>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    // Initialize local variables for constants
    const SoAFloatPrecision cutoffSquared = _cutoffSquared;
    const SoAFloatPrecision radius = _radius;

    // Loop over Particles
    for (unsigned int i = 0; i < soa.size(); i++) {
      const auto ownedStateI = ownedStatePtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      SoAFloatPrecision fxacc = 0.;
      SoAFloatPrecision fyacc = 0.;
      SoAFloatPrecision fzacc = 0.;

      // TODO: useMixing-case

      // TODO: parallelize with omp
      for (unsigned int j = i + 1; j < soa.size(); ++j) {
        // TODO: useMixing-case

        // Compute necessary values for computations of forces
        const auto ownedStateJ = ownedStatePtr[j];

        const SoAFloatPrecision drx = xptr[i] - xptr[j];
        const SoAFloatPrecision dry = yptr[i] - yptr[j];
        const SoAFloatPrecision drz = zptr[i] - zptr[j];

        const SoAFloatPrecision distSquared = drx * drx + dry * dry + drz * drz;
        const SoAFloatPrecision dist = std::sqrt(distSquared);
        SoAFloatPrecision overlap = 2. * radius - dist;
        if (useMixing) {
            const SoAFloatPrecision radius_i = _PPLibrary->getRadius(typeptr[i]);
            const SoAFloatPrecision radius_j = _PPLibrary->getRadius(typeptr[j]);
            overlap = radius_i + radius_j - dist;
        }

        const SoAFloatPrecision invDist = 1.0 / dist;
        const SoAFloatPrecision normalUnitX = drx * invDist;
        const SoAFloatPrecision normalUnitY = dry * invDist;
        const SoAFloatPrecision normalUnitZ = drz * invDist;

        const SoAFloatPrecision relVelX = vxptr[i] - vxptr[j];
        const SoAFloatPrecision relVelY = vyptr[i] - vyptr[j];
        const SoAFloatPrecision relVelZ = vzptr[i] - vzptr[j];

        const SoAFloatPrecision relVelDotNormalUnit =
            relVelX * normalUnitX + relVelY * normalUnitY + relVelZ * normalUnitZ;

        const SoAFloatPrecision tanRelVelX = relVelX - relVelDotNormalUnit * normalUnitX;
        const SoAFloatPrecision tanRelVelY = relVelY - relVelDotNormalUnit * normalUnitY;
        const SoAFloatPrecision tanRelVelZ = relVelZ - relVelDotNormalUnit * normalUnitZ;

        // Mask away if distance is too large or overlap is non-positive or any particle is a dummy
        const bool mask =
            distSquared <= cutoffSquared and overlap > 0 and ownedStateJ != autopas::OwnershipState::dummy;

        // Compute normal force
        const SoAFloatPrecision normalFX =
            _elasticStiffness * overlap * normalUnitX - _normalViscosity * relVelDotNormalUnit * normalUnitX;
        const SoAFloatPrecision normalFY =
            _elasticStiffness * overlap * normalUnitY - _normalViscosity * relVelDotNormalUnit * normalUnitY;
        const SoAFloatPrecision normalFZ =
            _elasticStiffness * overlap * normalUnitZ - _normalViscosity * relVelDotNormalUnit * normalUnitZ;

        // Compute tangential force
        SoAFloatPrecision tanFX = -_frictionViscosity * tanRelVelX;  // TODO: add tangential spring
        SoAFloatPrecision tanFY = -_frictionViscosity * tanRelVelY;
        SoAFloatPrecision tanFZ = -_frictionViscosity * tanRelVelZ;

        const SoAFloatPrecision tanFMag = std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ);
        const SoAFloatPrecision normalFMag = std::sqrt(normalFX * normalFX + normalFY * normalFY + normalFZ * normalFZ);
        const SoAFloatPrecision coulombLimit = _staticFrictionCoeff * (normalFMag + _adhesiveStiffness * overlap);

        if (tanFMag > coulombLimit) {
          const SoAFloatPrecision scale = _dynamicFrictionCoeff * (normalFMag + _adhesiveStiffness * overlap) / tanFMag;
          tanFX *= scale;
          tanFY *= scale;
          tanFZ *= scale;
        }

        // Compute total force
        fxacc += mask * (normalFX + tanFX);
        fyacc += mask * (normalFY + tanFY);
        fzacc += mask * (normalFZ + tanFZ);

        // Apply total force
        if (newton3) {
          // only if we use newton 3 here, we want to
          fxptr[j] -= mask * (normalFX + tanFX);
          fyptr[j] -= mask * (normalFY + tanFY);
          fzptr[j] -= mask * (normalFZ + tanFZ);
        }

        fxptr[i] += fxacc;
        fyptr[i] += fyacc;
        fzptr[i] += fzacc;
      }
    }
  }

  /**
   * @copydoc autopas::Functor::SoAFunctorPair()
   */
  void SoAFunctorPair(autopas::SoAView<SoAArraysType> soa1, autopas::SoAView<SoAArraysType> soa2, bool newton3) final {
    if (newton3) {
      SoAFunctorPairImpl<true>(soa1, soa2);
    } else {
      SoAFunctorPairImpl<false>(soa1, soa2);
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
  void setParticleProperties(SoAFloatPrecision radius) {  // TODO: add tangential spring etc.
    _radius = radius;
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {  // TODO
    return std::array<typename Particle::AttributeNames, 9>{
        Particle::AttributeNames::id,     Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,   Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {  // TODO
    return std::array<typename Particle::AttributeNames, 6>{
        Particle::AttributeNames::id,   Particle::AttributeNames::posX,   Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ, Particle::AttributeNames::typeId, Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {  // TODO
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
    if constexpr (calculateGloabls) {
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
    if (calculateGloabls) {
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
    if (not calculateGloabls) {
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
    if (not calculateGloabls) {
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

  [[nodiscard]] size_t getNumFLOPs() const override {
    // TODO
    return 0;
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

     // Initialize pointers to SoA Data for soa1
     const auto *const __restrict xptr1 = soa1.template begin<Particle::AttributeNames::posX>();
     const auto *const __restrict yptr1 = soa1.template begin<Particle::AttributeNames::posY>();
     const auto *const __restrict zptr1 = soa1.template begin<Particle::AttributeNames::posZ>();
     const auto *const __restrict ownedStatePtr1 = soa1.template begin<Particle::AttributeNames::ownershipState>();

     const auto *const __restrict vxptr1 = soa1.template begin<Particle::AttributeNames::velocityX>();
     const auto *const __restrict vyptr1 = soa1.template begin<Particle::AttributeNames::velocityY>();
     const auto *const __restrict vzptr1 = soa1.template begin<Particle::AttributeNames::velocityZ>();

     SoAFloatPrecision *const __restrict fxptr1 = soa1.template begin<Particle::AttributeNames::forceX>();
     SoAFloatPrecision *const __restrict fyptr1 = soa1.template begin<Particle::AttributeNames::forceY>();
     SoAFloatPrecision *const __restrict fzptr1 = soa1.template begin<Particle::AttributeNames::forceZ>();

     // Initialize pointers to SoA Data for soa2
     const auto *const __restrict xptr2 = soa2.template begin<Particle::AttributeNames::posX>();
     const auto *const __restrict yptr2 = soa2.template begin<Particle::AttributeNames::posY>();
     const auto *const __restrict zptr2 = soa2.template begin<Particle::AttributeNames::posZ>();
     const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

     const auto *const __restrict vxptr2 = soa2.template begin<Particle::AttributeNames::velocityX>();
     const auto *const __restrict vyptr2 = soa2.template begin<Particle::AttributeNames::velocityY>();
     const auto *const __restrict vzptr2 = soa2.template begin<Particle::AttributeNames::velocityZ>();

     SoAFloatPrecision *const __restrict fxptr2 = soa2.template begin<Particle::AttributeNames::forceX>();
     SoAFloatPrecision *const __restrict fyptr2 = soa2.template begin<Particle::AttributeNames::forceY>();
     SoAFloatPrecision *const __restrict fzptr2 = soa2.template begin<Particle::AttributeNames::forceZ>();

     [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle::AttributeNames::typeId>();
     [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle::AttributeNames::typeId>();

     // Initialize local variables for constants
     const SoAFloatPrecision cutoffSquared = _cutoffSquared;
     const SoAFloatPrecision radius = _radius;

     // Loop over Particles in soa1
     for (unsigned int i = 0; i < soa1.size(); ++i) {
       const auto ownedStateI = ownedStatePtr1[i];
       if (ownedStateI == autopas::OwnershipState::dummy) continue;

       SoAFloatPrecision fxacc = 0.;
       SoAFloatPrecision fyacc = 0.;
       SoAFloatPrecision fzacc = 0.;

       // Loop over Particles in soa2
       // TODO: parallelize with omp
       for (unsigned int j = 0; j < soa2.size(); ++j) {
         const auto ownedStateJ = ownedStatePtr2[j];

         const SoAFloatPrecision drx = xptr1[i] - xptr2[j];
         const SoAFloatPrecision dry = yptr1[i] - yptr2[j];
         const SoAFloatPrecision drz = zptr1[i] - zptr2[j];

         const SoAFloatPrecision distSquared = drx * drx + dry * dry + drz * drz;
         const SoAFloatPrecision dist = std::sqrt(distSquared);
         SoAFloatPrecision overlap = 2. * radius - dist;
         if (useMixing) {
           const SoAFloatPrecision radius_i = _PPLibrary->getRadius(typeptr1[i]);
           const SoAFloatPrecision radius_j = _PPLibrary->getRadius(typeptr2[j]);
           overlap = radius_i + radius_j - dist;
         }

         const SoAFloatPrecision invDist = 1.0 / dist;
         const SoAFloatPrecision normalUnitX = drx * invDist;
         const SoAFloatPrecision normalUnitY = dry * invDist;
         const SoAFloatPrecision normalUnitZ = drz * invDist;

         const SoAFloatPrecision relVelX = vxptr1[i] - vxptr2[j];
         const SoAFloatPrecision relVelY = vyptr1[i] - vyptr2[j];
         const SoAFloatPrecision relVelZ = vzptr1[i] - vzptr2[j];

         const SoAFloatPrecision relVelDotNormalUnit = relVelX * normalUnitX + relVelY * normalUnitY + relVelZ * normalUnitZ;

         const SoAFloatPrecision tanRelVelX = relVelX - relVelDotNormalUnit * normalUnitX;
         const SoAFloatPrecision tanRelVelY = relVelY - relVelDotNormalUnit * normalUnitY;
         const SoAFloatPrecision tanRelVelZ = relVelZ - relVelDotNormalUnit * normalUnitZ;

         // Compute normal force
         const SoAFloatPrecision normalFX = _elasticStiffness * overlap * normalUnitX - _normalViscosity * relVelDotNormalUnit * normalUnitX;
         const SoAFloatPrecision normalFY = _elasticStiffness * overlap * normalUnitY - _normalViscosity * relVelDotNormalUnit * normalUnitY;
         const SoAFloatPrecision normalFZ = _elasticStiffness * overlap * normalUnitZ - _normalViscosity * relVelDotNormalUnit * normalUnitZ;

         // Compute tangential force
         SoAFloatPrecision tanFX = -_frictionViscosity * tanRelVelX;  // TODO: add tangential spring
         SoAFloatPrecision tanFY = -_frictionViscosity * tanRelVelY;
         SoAFloatPrecision tanFZ = -_frictionViscosity * tanRelVelZ;

         const SoAFloatPrecision tanFMag = std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ);
         const SoAFloatPrecision normalFMag = std::sqrt(normalFX * normalFX + normalFY * normalFY + normalFZ * normalFZ);
         const SoAFloatPrecision coulombLimit = _staticFrictionCoeff * (normalFMag + _adhesiveStiffness * overlap);

         if (tanFMag > coulombLimit) {
           const SoAFloatPrecision scale = _dynamicFrictionCoeff * (normalFMag + _adhesiveStiffness * overlap) / tanFMag;
           tanFX *= scale;
           tanFY *= scale;
           tanFZ *= scale;
         }

         // Compute total force
         fxacc += normalFX + tanFX;
         fyacc += normalFY + tanFY;
         fzacc += normalFZ + tanFZ;

         // Apply total force
         if (newton3) {
           fxptr2[j] -= normalFX + tanFX;
           fyptr2[j] -= normalFY + tanFY;
           fzptr2[j] -= normalFZ + tanFZ;
         }
       }

       fxptr1[i] += fxacc;
       fyptr1[i] += fyacc;
       fzptr1[i] += fzacc;
     }

   }

  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    autopas::utils::ExceptionHandler::exception("DEMFunctor::SoAFunctorVerletImpl() is not implemented.");
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
  const double _elasticStiffness;
  const double _adhesiveStiffness;
  const double _frictionStiffness;
  const double _normalViscosity;
  const double _frictionViscosity;
  const double _staticFrictionCoeff;
  const double _dynamicFrictionCoeff;
  // not const because they might be reset through PPL
  double _epsilon, _sigma, _radius = 0;

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