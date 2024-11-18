/**
 * @file DEMFunctor.h
 * @author J. Kim
 * @date 04/11/24
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"

namespace demLib {

/**
 * A functor to handle the interactions between particles using Discrete Element Method (DEM).
 * This functor assumes that duplicated calculations are always happening, which is characteristic for a Full-Shell
 * scheme.
 * @tparam Particle The type of particle.
 * @tparam useMixing Switch for the functor to be used with multiple particle types.
 * If set to false, _epsilon6 and _sigma need to be set and the constructor with PPL can be omitted.
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
   * @param elasticStiffness
   * @param adhesiveStiffness
   * @param frictionStiffness
   * @param normalViscosity
   * @param frictionViscosity
   * @param rollingViscosity
   * @param staticFrictionCoeff
   * @param dynamicFrictionCoeff
   * @param rollingFrictionCoeff
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit DEMFunctor(double cutoff, double elasticStiffness, double adhesiveStiffness, double frictionStiffness,
                      double normalViscosity, double frictionViscosity, double rollingViscosity,
                      double staticFrictionCoeff, double dynamicFrictionCoeff, double rollingFrictionCoeff,
                      void * /*dummy*/)
      : autopas::Functor<Particle,
                         DEMFunctor<Particle, useMixing, useNewton3, calculateGloabls, countFLOPs, relevantForTuning>>(
            cutoff),
        _radius{cutoff / 5.},  // default initialization
        _cutoff{cutoff},
        _elasticStiffness{elasticStiffness},
        _adhesiveStiffness{adhesiveStiffness},
        _frictionStiffness{frictionStiffness},
        _normalViscosity{normalViscosity},
        _frictionViscosity{frictionViscosity},
        _rollingViscosity{rollingViscosity},
        _staticFrictionCoeff{staticFrictionCoeff},
        _dynamicFrictionCoeff{dynamicFrictionCoeff},
        _rollingFrictionCoeff{rollingFrictionCoeff},
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
   *
   * @param cutoff
   */
  explicit DEMFunctor(double cutoff) : DEMFunctor(cutoff, 5., 2.5, 1., 1e-2, 1e-1, 1e-3, 50., 25., 15., nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Minimal Constructor for Functor with mixing disabled. Other Parameters are set with default values.
   * Values taken from https://www2.msm.ctw.utwente.nl/sluding/PAPERS/Alert_Luding2008.pdf page 15.
   *
   * @note Only to be used with mixing == true.
   *
   * @param cutoff
   * @param particlePropertiesLibrary
   */
  explicit DEMFunctor(double cutoff, ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : DEMFunctor(cutoff, 5., 2.5, 1., 1e-2, 1e-1, 1e-3, 50., 25., 15., nullptr) {
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
   * @param elasticStiffness
   * @param adhesiveStiffness
   * @param frictionStiffness
   * @param normalViscosity
   * @param frictionViscosity
   * @param rollingViscosity
   * @param staticFrictionCoeff
   * @param dynamicFrictionCoeff
   * @param rollingFrictionCoeff
   */
  explicit DEMFunctor(double cutoff, double elasticStiffness, double adhesiveStiffness, double frictionStiffness,
                      double normalViscosity, double frictionViscosity, double rollingViscosity,
                      double staticFrictionCoeff, double dynamicFrictionCoeff, double rollingFrictionCoeff)
      : DEMFunctor(cutoff, elasticStiffness, adhesiveStiffness, frictionStiffness, normalViscosity, frictionViscosity,
                   rollingViscosity, staticFrictionCoeff, dynamicFrictionCoeff, rollingFrictionCoeff, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like radius.
   * @param cutoff
   * @param elasticStiffness
   * @param adhesiveStiffness
   * @param frictionStiffness
   * @param normalViscosity
   * @param frictionViscosity
   * @param rollingViscosity
   * @param staticFrictionCoeff
   * @param dynamicFrictionCoeff
   * @param rollingFrictionCoeff
   * @param particlePropertiesLibrary
   */
  explicit DEMFunctor(double cutoff, double elasticStiffness, double adhesiveStiffness, double frictionStiffness,
                      double normalViscosity, double frictionViscosity, double rollingViscosity,
                      double staticFrictionCoeff, double dynamicFrictionCoeff, double rollingFrictionCoeff,
                      ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : DEMFunctor(cutoff, elasticStiffness, adhesiveStiffness, frictionStiffness, normalViscosity, frictionViscosity,
                   rollingViscosity, staticFrictionCoeff, dynamicFrictionCoeff, rollingFrictionCoeff, nullptr) {
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
    const double cutoff = _cutoff;
    const std::array<double, 3> displacement = i.getR() - j.getR();
    const double dist = autopas::utils::ArrayMath::L2Norm(displacement);

    if (dist > cutoff) return;  // cutoff check

    // Retrieve particle-specific parameters, calculate generally necessary values
    const auto [sigma, epsilon6, radiusI, radiusJ] = computeMaterialProperties(i, j);
    const std::array<double, 3> normalUnit = displacement / dist;
    const double overlap = radiusI + radiusJ - dist;
    const double radiusIReduced = radiusI - overlap / 2.;
    const double radiusJReduced = radiusJ - overlap / 2.;
    const std::array<double, 3> relVel = i.getV() - j.getV();
    const double normalRelVelMag = autopas::utils::ArrayMath::dot(normalUnit, relVel);

    // Compute Forces
    // Compute normal forces
    const double normalContactFMag = computeNormalContactFMag(overlap, normalRelVelMag);
    const double normalVdWFMag = computeNormalVdWFMag(overlap, dist, sigma, epsilon6, cutoff);
    const double normalFMag = normalContactFMag + normalVdWFMag;
    const std::array<double, 3> normalF = autopas::utils::ArrayMath::mulScalar(normalUnit, normalFMag);

    // Compute tangential force
    const std::array<double, 3> tanF = computeTangentialForce(overlap, i, j, radiusIReduced, radiusJReduced, normalUnit,
                                                              normalRelVelMag, normalContactFMag);

    // Compute total force
    const std::array<double, 3> totalF = {0., 0., 0.}; // Todo: change

    // Apply forces
    i.addF(totalF);
    if (newton3) {
      j.subF(totalF);
    }

    // Compute Torques
    // Compute frictional torque
    const std::array<double, 3> frictionQI = computeFrictionTorqueI(overlap, radiusIReduced, normalUnit, tanF);
    const std::array<double, 3> rollingQI = computeRollingTorqueI(overlap, radiusIReduced, radiusJReduced, i, j, normalUnit, normalContactFMag);

    // Apply torques
    i.addTorque(frictionQI + rollingQI);
    if (newton3) {
      j.addTorque((frictionQI * (radiusJReduced / radiusIReduced)) - rollingQI);
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

    const auto *const __restrict angularVelXptr = soa.template begin<Particle::AttributeNames::angularVelX>();
    const auto *const __restrict angularVelYptr = soa.template begin<Particle::AttributeNames::angularVelY>();
    const auto *const __restrict angularVelZptr = soa.template begin<Particle::AttributeNames::angularVelZ>();

    SoAFloatPrecision *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    SoAFloatPrecision *const __restrict qXptr = soa.template begin<Particle::AttributeNames::torqueX>();
    SoAFloatPrecision *const __restrict qYptr = soa.template begin<Particle::AttributeNames::torqueY>();
    SoAFloatPrecision *const __restrict qZptr = soa.template begin<Particle::AttributeNames::torqueZ>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    // Initialize local variables for constants
    const SoAFloatPrecision cutoff = _cutoff;
    const SoAFloatPrecision constRadius = _radius;
    const SoAFloatPrecision const_sigma = _sigma;
    const SoAFloatPrecision const_epsilon6 = _epsilon6;

    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmas;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> epsilon6s;
    if constexpr (useMixing) {
      // Preload all sigma and epsilons for next vectorized region.
      // Not preloading and directly using the values, will produce worse results.
      sigmas.resize(soa.size());
      epsilon6s.resize(soa.size());
    }

    // Loop over Particles
    for (unsigned int i = 0; i < soa.size(); i++) {
      const auto ownedStateI = ownedStatePtr[i];
      if (ownedStateI == autopas::OwnershipState::dummy) {
        continue;
      }

      SoAFloatPrecision fxacc = 0.;
      SoAFloatPrecision fyacc = 0.;
      SoAFloatPrecision fzacc = 0.;

      SoAFloatPrecision qXacc = 0.;
      SoAFloatPrecision qYacc = 0.;
      SoAFloatPrecision qZacc = 0.;

      double radiusI = constRadius;

      if constexpr (useMixing) {
        radiusI = _PPLibrary->getRadius(typeptr[i]);
        for (unsigned int j = 0; j < soa.size(); ++j) {
          auto mixingData = _PPLibrary->getMixingData(typeptr[i], typeptr[j]);
          sigmas[j] = mixingData.sigma;
          epsilon6s[j] = mixingData.epsilon6;
        }
      }

      // TODO: parallelize with omp
      for (unsigned int j = i + 1; j < soa.size(); ++j) {
        // Compute necessary values for computations of forces
        const auto ownedStateJ = ownedStatePtr[j];

        const SoAFloatPrecision drx = xptr[i] - xptr[j];
        const SoAFloatPrecision dry = yptr[i] - yptr[j];
        const SoAFloatPrecision drz = zptr[i] - zptr[j];

        SoAFloatPrecision sigma = const_sigma;
        SoAFloatPrecision epsilon6 = const_epsilon6;
        double radiusJ = constRadius;

        if (useMixing) {
          sigma = sigmas[j];
          epsilon6 = epsilon6s[j];
          radiusJ = _PPLibrary->getRadius(typeptr[j]);
        }

        const SoAFloatPrecision distSquared = drx * drx + dry * dry + drz * drz;
        const SoAFloatPrecision dist = std::sqrt(distSquared);
        SoAFloatPrecision overlap = 2. * (radiusI + radiusJ) - dist;

        const SoAFloatPrecision invDist = 1. / dist;
        const SoAFloatPrecision invDistSquared = 1. / distSquared;
        const SoAFloatPrecision normalUnitX = drx * invDist;
        const SoAFloatPrecision normalUnitY = dry * invDist;
        const SoAFloatPrecision normalUnitZ = drz * invDist;

        const SoAFloatPrecision radiusIReduced = radiusI - overlap / 2.;
        const SoAFloatPrecision radiusJReduced = radiusJ - overlap / 2.;

        const SoAFloatPrecision relVelAngularIX =
            radiusIReduced * (normalUnitY * angularVelZptr[i] - normalUnitZ * angularVelYptr[i]);
        const SoAFloatPrecision relVelAngularIY =
            radiusIReduced * (normalUnitZ * angularVelXptr[i] - normalUnitX * angularVelZptr[i]);
        const SoAFloatPrecision relVelAngularIZ =
            radiusIReduced * (normalUnitX * angularVelYptr[i] - normalUnitY * angularVelXptr[i]);
        const SoAFloatPrecision relVelAngularJX =
            radiusJReduced * (normalUnitY * angularVelZptr[j] - normalUnitZ * angularVelYptr[j]);
        const SoAFloatPrecision relVelAngularJY =
            radiusJReduced * (normalUnitZ * angularVelXptr[j] - normalUnitX * angularVelZptr[j]);
        const SoAFloatPrecision relVelAngularJZ =
            radiusJReduced * (normalUnitX * angularVelYptr[j] - normalUnitY * angularVelXptr[j]);

        const SoAFloatPrecision relVelX = vxptr[i] - vxptr[j] + relVelAngularIX + relVelAngularJX;
        const SoAFloatPrecision relVelY = vyptr[i] - vyptr[j] + relVelAngularIY + relVelAngularJY;
        const SoAFloatPrecision relVelZ = vzptr[i] - vzptr[j] + relVelAngularIZ + relVelAngularJZ;

        const SoAFloatPrecision relVelDotNormalUnit =
            relVelX * normalUnitX + relVelY * normalUnitY + relVelZ * normalUnitZ;

        const SoAFloatPrecision tanRelVelX = relVelX - relVelDotNormalUnit * normalUnitX;
        const SoAFloatPrecision tanRelVelY = relVelY - relVelDotNormalUnit * normalUnitY;
        const SoAFloatPrecision tanRelVelZ = relVelZ - relVelDotNormalUnit * normalUnitZ;

        const SoAFloatPrecision invSigma = 1. / sigma;
        const SoAFloatPrecision lj2 = sigma * sigma * invDistSquared;
        const SoAFloatPrecision lj7 = lj2 * lj2 * lj2 * sigma * invDist;
        const SoAFloatPrecision ljCutoff2 = sigma * sigma / (cutoff * cutoff);
        const SoAFloatPrecision ljCutoff7 = ljCutoff2 * ljCutoff2 * ljCutoff2 * sigma / cutoff;

        // Mask away if distance is too large or overlap is non-positive or any particle is a dummy
        const bool cutOffMask = dist <= _cutoff and ownedStateJ != autopas::OwnershipState::dummy;
        const bool overlapIsPositive = overlap > 0;

        // Compute normal force as sum of contact force and long-range force (VdW).
        const SoAFloatPrecision normalContactFMag =
            _elasticStiffness * overlap - _normalViscosity * relVelDotNormalUnit;
        const SoAFloatPrecision normalVdWFMag = -1. * (not overlapIsPositive) * epsilon6 * invSigma * (lj7 - ljCutoff7);
        const SoAFloatPrecision normalFMag = overlapIsPositive * normalContactFMag + normalVdWFMag;

        const SoAFloatPrecision normalFX = normalFMag * normalUnitX;
        const SoAFloatPrecision normalFY = normalFMag * normalUnitY;
        const SoAFloatPrecision normalFZ = normalFMag * normalUnitZ;

        // Compute tangential force
        SoAFloatPrecision tanFX = -_frictionViscosity * tanRelVelX;
        SoAFloatPrecision tanFY = -_frictionViscosity * tanRelVelY;
        SoAFloatPrecision tanFZ = -_frictionViscosity * tanRelVelZ;

        const SoAFloatPrecision tanFMag = std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ);
        const SoAFloatPrecision coulombLimit =
            _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);

        if (tanFMag > coulombLimit) {
          const SoAFloatPrecision scale =
              _dynamicFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap) / tanFMag;
          tanFX *= scale;
          tanFY *= scale;
          tanFZ *= scale;
        }

        // Compute total force
        const SoAFloatPrecision totalFX = cutOffMask * (normalFX + overlapIsPositive * tanFX);
        const SoAFloatPrecision totalFY = cutOffMask * (normalFY + overlapIsPositive * tanFY);
        const SoAFloatPrecision totalFZ = cutOffMask * (normalFZ + overlapIsPositive * tanFZ);

        fxacc += totalFX;
        fyacc += totalFY;
        fzacc += totalFZ;

        // Apply total force
        if (newton3) {
          // only if we use newton 3 here, we want to
          fxptr[j] -= totalFX;
          fyptr[j] -= totalFY;
          fzptr[j] -= totalFZ;
        }

        // Compute torques
        // Compute frictional torque
        const SoAFloatPrecision frictionQIX = -radiusIReduced * (normalUnitY * tanFZ - normalUnitZ * tanFY);
        const SoAFloatPrecision frictionQIY = -radiusIReduced * (normalUnitZ * tanFX - normalUnitX * tanFZ);
        const SoAFloatPrecision frictionQIZ = -radiusIReduced * (normalUnitX * tanFY - normalUnitY * tanFX);

        qXacc += frictionQIX;
        qYacc += frictionQIY;
        qZacc += frictionQIZ;

        // Apply torques
        if (newton3) {
          qXptr[j] += (radiusJReduced / radiusIReduced) * frictionQIX;
          qYptr[j] += (radiusJReduced / radiusIReduced) * frictionQIY;
          qZptr[j] += (radiusJReduced / radiusIReduced) * frictionQIZ;
        }
      }  // end of j loop

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;

      qXptr[i] += qXacc;
      qYptr[i] += qYacc;
      qZptr[i] += qZacc;
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
  void setParticleProperties(SoAFloatPrecision epsilon6, SoAFloatPrecision sigma, SoAFloatPrecision radius) {
    _epsilon6 = epsilon6;
    _sigma = sigma;
    _radius = radius;
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 12>{
        Particle::AttributeNames::id,      Particle::AttributeNames::posX,    Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,    Particle::AttributeNames::forceX,  Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ,  Particle::AttributeNames::torqueX, Particle::AttributeNames::torqueY,
        Particle::AttributeNames::torqueZ, Particle::AttributeNames::typeId,  Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 12>{
        Particle::AttributeNames::id,      Particle::AttributeNames::posX,    Particle::AttributeNames::posY,
        Particle::AttributeNames::posZ,    Particle::AttributeNames::forceX,  Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ,  Particle::AttributeNames::torqueX, Particle::AttributeNames::torqueY,
        Particle::AttributeNames::torqueZ, Particle::AttributeNames::typeId,  Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 6>{
        Particle::AttributeNames::forceX,  Particle::AttributeNames::forceY,  Particle::AttributeNames::forceZ,
        Particle::AttributeNames::torqueX, Particle::AttributeNames::torqueY, Particle::AttributeNames::torqueZ};
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

    const auto *const __restrict angularVelXptr1 = soa1.template begin<Particle::AttributeNames::angularVelX>();
    const auto *const __restrict angularVelYptr1 = soa1.template begin<Particle::AttributeNames::angularVelY>();
    const auto *const __restrict angularVelZptr1 = soa1.template begin<Particle::AttributeNames::angularVelZ>();

    SoAFloatPrecision *const __restrict fxptr1 = soa1.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr1 = soa1.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr1 = soa1.template begin<Particle::AttributeNames::forceZ>();

    SoAFloatPrecision *const __restrict qXptr1 = soa1.template begin<Particle::AttributeNames::torqueX>();
    SoAFloatPrecision *const __restrict qYptr1 = soa1.template begin<Particle::AttributeNames::torqueY>();
    SoAFloatPrecision *const __restrict qZptr1 = soa1.template begin<Particle::AttributeNames::torqueZ>();

    // Initialize pointers to SoA Data for soa2
    const auto *const __restrict xptr2 = soa2.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr2 = soa2.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr2 = soa2.template begin<Particle::AttributeNames::posZ>();
    const auto *const __restrict ownedStatePtr2 = soa2.template begin<Particle::AttributeNames::ownershipState>();

    const auto *const __restrict vxptr2 = soa2.template begin<Particle::AttributeNames::velocityX>();
    const auto *const __restrict vyptr2 = soa2.template begin<Particle::AttributeNames::velocityY>();
    const auto *const __restrict vzptr2 = soa2.template begin<Particle::AttributeNames::velocityZ>();

    const auto *const __restrict angularVelXptr2 = soa2.template begin<Particle::AttributeNames::angularVelX>();
    const auto *const __restrict angularVelYptr2 = soa2.template begin<Particle::AttributeNames::angularVelY>();
    const auto *const __restrict angularVelZptr2 = soa2.template begin<Particle::AttributeNames::angularVelZ>();

    SoAFloatPrecision *const __restrict fxptr2 = soa2.template begin<Particle::AttributeNames::forceX>();
    SoAFloatPrecision *const __restrict fyptr2 = soa2.template begin<Particle::AttributeNames::forceY>();
    SoAFloatPrecision *const __restrict fzptr2 = soa2.template begin<Particle::AttributeNames::forceZ>();

    SoAFloatPrecision *const __restrict qXptr2 = soa2.template begin<Particle::AttributeNames::torqueX>();
    SoAFloatPrecision *const __restrict qYptr2 = soa2.template begin<Particle::AttributeNames::torqueY>();
    SoAFloatPrecision *const __restrict qZptr2 = soa2.template begin<Particle::AttributeNames::torqueZ>();

    [[maybe_unused]] auto *const __restrict typeptr1 = soa1.template begin<Particle::AttributeNames::typeId>();
    [[maybe_unused]] auto *const __restrict typeptr2 = soa2.template begin<Particle::AttributeNames::typeId>();

    // Initialize local variables for constants
    const SoAFloatPrecision cutoff = _cutoff;
    const SoAFloatPrecision radius = _radius;
    SoAFloatPrecision sigma = _sigma;
    SoAFloatPrecision epsilon6 = _epsilon6;

    // preload all sigma and epsilons for next vectorized region
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> sigmas;
    std::vector<SoAFloatPrecision, autopas::AlignedAllocator<SoAFloatPrecision>> epsilon6s;
    if constexpr (useMixing) {
      sigmas.resize(soa2.size());
      epsilon6s.resize(soa2.size());
    }

    // Loop over Particles in soa1
    for (unsigned int i = 0; i < soa1.size(); ++i) {
      const auto ownedStateI = ownedStatePtr1[i];
      if (ownedStateI == autopas::OwnershipState::dummy) continue;

      SoAFloatPrecision fxacc = 0.;
      SoAFloatPrecision fyacc = 0.;
      SoAFloatPrecision fzacc = 0.;

      SoAFloatPrecision qXacc = 0.;
      SoAFloatPrecision qYacc = 0.;
      SoAFloatPrecision qZacc = 0.;

      SoAFloatPrecision radiusI = radius;
      // preload all sigma and epsilons for next vectorized region
      if constexpr (useMixing) {
        radiusI = _PPLibrary->getRadius(typeptr1[i]);
        for (unsigned int j = 0; j < soa2.size(); ++j) {
          sigmas[j] = _PPLibrary->getMixingSigma(typeptr1[i], typeptr2[j]);
          epsilon6s[j] = _PPLibrary->getMixing6Epsilon(typeptr1[i], typeptr2[j]);
        }
      }

      // Loop over Particles in soa2
      // TODO: parallelize with omp
      for (unsigned int j = 0; j < soa2.size(); ++j) {
        const auto ownedStateJ = ownedStatePtr2[j];

        const SoAFloatPrecision drx = xptr1[i] - xptr2[j];
        const SoAFloatPrecision dry = yptr1[i] - yptr2[j];
        const SoAFloatPrecision drz = zptr1[i] - zptr2[j];

        const SoAFloatPrecision distSquared = drx * drx + dry * dry + drz * drz;
        const SoAFloatPrecision dist = std::sqrt(distSquared);
        SoAFloatPrecision radiusJ = radius;
        if (useMixing) {
          radiusJ = _PPLibrary->getRadius(typeptr2[j]);
          sigma = sigmas[j];
          epsilon6 = epsilon6s[j];
        }
        const SoAFloatPrecision overlap = 2. * (radiusI + radiusJ) - dist;

        const SoAFloatPrecision invDist = 1.0 / dist;
        const SoAFloatPrecision invDistSquared = 1.0 / distSquared;
        const SoAFloatPrecision normalUnitX = drx * invDist;
        const SoAFloatPrecision normalUnitY = dry * invDist;
        const SoAFloatPrecision normalUnitZ = drz * invDist;

        const SoAFloatPrecision radiusIReduced = radiusI - overlap / 2.;
        const SoAFloatPrecision radiusJReduced = radiusJ - overlap / 2.;

        const SoAFloatPrecision relVelAngularIX =
            radiusIReduced * (normalUnitY * angularVelZptr1[i] - normalUnitZ * angularVelYptr1[i]);
        const SoAFloatPrecision relVelAngularIY =
            radiusIReduced * (normalUnitZ * angularVelXptr1[i] - normalUnitX * angularVelZptr1[i]);
        const SoAFloatPrecision relVelAngularIZ =
            radiusIReduced * (normalUnitX * angularVelYptr1[i] - normalUnitY * angularVelXptr1[i]);
        const SoAFloatPrecision relVelAngularJX =
            radiusJReduced * (normalUnitY * angularVelZptr2[j] - normalUnitZ * angularVelYptr2[j]);
        const SoAFloatPrecision relVelAngularJY =
            radiusJReduced * (normalUnitZ * angularVelXptr2[j] - normalUnitX * angularVelZptr2[j]);
        const SoAFloatPrecision relVelAngularJZ =
            radiusJReduced * (normalUnitX * angularVelYptr2[j] - normalUnitY * angularVelXptr2[j]);

        const SoAFloatPrecision relVelX = vxptr1[i] - vxptr2[j] + relVelAngularIX + relVelAngularJX;
        const SoAFloatPrecision relVelY = vyptr1[i] - vyptr2[j] + relVelAngularIY + relVelAngularJY;
        const SoAFloatPrecision relVelZ = vzptr1[i] - vzptr2[j] + relVelAngularIZ + relVelAngularJZ;

        const SoAFloatPrecision relVelDotNormalUnit =
            relVelX * normalUnitX + relVelY * normalUnitY + relVelZ * normalUnitZ;

        const SoAFloatPrecision tanRelVelX = relVelX - relVelDotNormalUnit * normalUnitX;
        const SoAFloatPrecision tanRelVelY = relVelY - relVelDotNormalUnit * normalUnitY;
        const SoAFloatPrecision tanRelVelZ = relVelZ - relVelDotNormalUnit * normalUnitZ;

        const SoAFloatPrecision invSigma = 1. / sigma;
        const SoAFloatPrecision lj2 = sigma * sigma * invDistSquared;
        const SoAFloatPrecision lj7 = lj2 * lj2 * lj2 * sigma * invDist;
        const SoAFloatPrecision ljCutoff2 = sigma * sigma / (cutoff * cutoff);
        const SoAFloatPrecision ljCutoff7 = ljCutoff2 * ljCutoff2 * ljCutoff2 * sigma / cutoff;

        // Mask away if distance is too large or overlap is non-positive or any particle is a dummy
        const bool cutOffMask = dist <= _cutoff and ownedStateJ != autopas::OwnershipState::dummy;
        const bool overlapIsPositive = overlap > 0;

        // Compute normal force
        const SoAFloatPrecision normalContactFMag =
            _elasticStiffness * overlap - _normalViscosity * relVelDotNormalUnit;
        const SoAFloatPrecision normalVdWFMag = -1. * (not overlapIsPositive) * epsilon6 * invSigma * (lj7 - ljCutoff7);
        const SoAFloatPrecision normalFMag = overlapIsPositive * normalContactFMag + normalVdWFMag;

        const SoAFloatPrecision normalFX = normalFMag * normalUnitX;
        const SoAFloatPrecision normalFY = normalFMag * normalUnitY;
        const SoAFloatPrecision normalFZ = normalFMag * normalUnitZ;

        // Compute tangential force
        SoAFloatPrecision tanFX = -_frictionViscosity * tanRelVelX;  // TODO: add tangential spring
        SoAFloatPrecision tanFY = -_frictionViscosity * tanRelVelY;
        SoAFloatPrecision tanFZ = -_frictionViscosity * tanRelVelZ;

        const SoAFloatPrecision tanFMag = std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ);
        const SoAFloatPrecision coulombLimit =
            _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);

        if (tanFMag > coulombLimit) {
          const SoAFloatPrecision scale =
              _dynamicFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap) / tanFMag;
          tanFX *= scale;
          tanFY *= scale;
          tanFZ *= scale;
        }

        // Compute total force
        const SoAFloatPrecision totalFX = cutOffMask * (normalFX + overlapIsPositive * tanFX);
        const SoAFloatPrecision totalFY = cutOffMask * (normalFY + overlapIsPositive * tanFY);
        const SoAFloatPrecision totalFZ = cutOffMask * (normalFZ + overlapIsPositive * tanFZ);

        fxacc += totalFX;
        fyacc += totalFY;
        fzacc += totalFZ;

        // Apply total force
        if (newton3) {
          fxptr2[j] -= totalFX;
          fyptr2[j] -= totalFY;
          fzptr2[j] -= totalFZ;
        }

        // Compute torques
        // Compute frictional torque
        const SoAFloatPrecision frictionQIX = -radiusI * (normalUnitY * tanFZ - normalUnitZ * tanFY);
        const SoAFloatPrecision frictionQIY = -radiusIReduced * (normalUnitZ * tanFX - normalUnitX * tanFZ);
        const SoAFloatPrecision frictionQIZ = -radiusIReduced * (normalUnitX * tanFY - normalUnitY * tanFX);

        qXacc += frictionQIX;
        qYacc += frictionQIY;
        qZacc += frictionQIZ;

        if (newton3) {
          qXptr2[j] += (radiusJReduced / radiusIReduced) * frictionQIX;
          qYptr2[j] += (radiusJReduced / radiusIReduced) * frictionQIY;
          qZptr2[j] += (radiusJReduced / radiusIReduced) * frictionQIZ;
        }

      }  // end of j loop

      fxptr1[i] += fxacc;
      fyptr1[i] += fyacc;
      fzptr1[i] += fzacc;

      qXptr1[i] += qXacc;
      qYptr1[i] += qYacc;
      qZptr1[i] += qZacc;
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

  const double _cutoff;
  const double _elasticStiffness;
  const double _adhesiveStiffness;
  const double _frictionStiffness;
  const double _normalViscosity;
  const double _frictionViscosity;
  const double _rollingViscosity;
  const double _staticFrictionCoeff;
  const double _dynamicFrictionCoeff;
  const double _rollingFrictionCoeff;
  // not const because they might be reset through PPL
  double _epsilon6, _sigma, _radius = 0;

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

  // ----------------------------------------Helper Methods----------------------------------------

  std::tuple<double, double, double, double> computeMaterialProperties(const Particle &i, const Particle &j) const {
    if (!useMixing) {
      return {_sigma, _epsilon6, _radius, _radius};
    }
    double sigma = _PPLibrary->getMixingSigma(i.getTypeId(), j.getTypeId());
    double epsilon6 = _PPLibrary->getMixing6Epsilon(i.getTypeId(), j.getTypeId());
    double radiusI = _PPLibrary->getRadius(i.getTypeId());
    double radiusJ = _PPLibrary->getRadius(j.getTypeId());
    return {sigma, epsilon6, radiusI, radiusJ};
  }

  double computeNormalContactFMag(const double overlap, const double normalRelVelMag) {
    return overlap > 0 ? _elasticStiffness * overlap - _normalViscosity * normalRelVelMag : 0;
  }

  double computeNormalVdWFMag(const double overlap, const double dist, const double sigma, const double epsilon6,
                              const double cutoff) {
    if (overlap > 0) {
      return 0;
    }
    const double invSigma = 1. / sigma;
    const double lj2 = (sigma * sigma) / (dist * dist);
    const double lj7 = lj2 * lj2 * lj2 * (sigma / dist);
    const double ljCutoff2 = (sigma * sigma) / (cutoff * cutoff);
    const double ljCutoff7 = ljCutoff2 * ljCutoff2 * ljCutoff2 * (sigma / cutoff);
    return -epsilon6 * invSigma * (lj7 - ljCutoff7);
  }

  std::array<double, 3> computeTangentialForce(const double overlap, const Particle &i, const Particle &j,
                                               const double radiusIReduced, const double radiusJReduced,
                                               const std::array<double, 3> &normalUnit, const double normalRelVelMag,
                                               const double normalContactFMag) {
    using namespace autopas::utils::ArrayMath::literals;

    if (overlap <= 0) {
      return {0, 0, 0};
    }
    const std::array<double, 3> tanRelVel =
        i.getV() - j.getV() + autopas::utils::ArrayMath::cross(normalUnit * radiusIReduced, i.getAngularVel()) +
        autopas::utils::ArrayMath::cross(normalUnit * radiusJReduced, j.getAngularVel());
    const std::array<double, 3> normalRelVel = normalUnit * normalRelVelMag;
    const std::array<double, 3> tanVel = tanRelVel - normalRelVel;
    const double coulombLimit = _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);

    const std::array<double, 3> tanF = tanVel * (-_frictionViscosity);
    const double tanFMag = autopas::utils::ArrayMath::L2Norm(tanF);
    if (tanFMag > coulombLimit) {
      const std::array<double, 3> tanFUnit = tanF / tanFMag;
      const double scale = _dynamicFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);
      return tanFUnit * scale;
    } else {
      return tanF;
    }
  }

  std::array<double, 3> computeFrictionTorqueI(const double overlap, const double radiusIReduced,
                                               const std::array<double, 3> &normalUnit,
                                               const std::array<double, 3> &tanF) {
    using namespace autopas::utils::ArrayMath::literals;
    if (overlap <= 0) {
      return {0, 0, 0};
    }
    return autopas::utils::ArrayMath::cross(normalUnit * (-radiusIReduced), tanF);
  }

  std::array<double, 3> computeRollingTorqueI(const double overlap, const double radiusIReduced, const double radiusJReduced, const Particle &i, const Particle &j, cont std::array<double, 3> &normalUnit, const double normalContactFMag) {
    using namespace autopas::utils::ArrayMath::literals;
    if (overlap <= 0) {
      return {0, 0, 0};
    }
    const double radiusReduced = radiusIReduced * radiusJReduced / (radiusIReduced + radiusJReduced);
    const std::array<double, 3> rollingRelVel = -radiusReduced * (autopas::utils::ArrayMath::cross(normalUnit, i.getAngularVel()) - autopas::utils::ArrayMath::cross(normalUnit, j.getAngularVel()));
    const std::array<double, 3> rollingF = rollingRelVel * (-_rollingViscosity);
    const double rollingFMag = autopas::utils::ArrayMath::L2Norm(rollingF);

    const double coulombLimit = _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);
    if (rollingFMag > coulombLimit) {
      const std::array<double, 3> rollingFUnit = rollingF / rollingFMag;
      const double scale = _rollingFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);
      return rollingFUnit * scale;
    } else {
      return rollingF;
    }
  }

};
}  // namespace demLib