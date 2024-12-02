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
   * @param torsionViscosity
   * @param staticFrictionCoeff
   * @param dynamicFrictionCoeff
   * @param rollingFrictionCoeff
   * @param torsionFrictionCo
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit DEMFunctor(double cutoff, double elasticStiffness, double adhesiveStiffness, double frictionStiffness,
                      double normalViscosity, double frictionViscosity, double rollingViscosity,
                      double torsionViscosity, double staticFrictionCoeff, double dynamicFrictionCoeff,
                      double rollingFrictionCoeff, double torsionFrictionCoeff, void * /*dummy*/)
      : autopas::Functor<Particle,
                         DEMFunctor<Particle, useMixing, useNewton3, calculateGloabls, countFLOPs, relevantForTuning>>(
            cutoff),
        _radius{cutoff / 5.},  // default initialization, can be changed by ParticlePropertiesLibrary
        _cutoff{cutoff},
        _elasticStiffness{elasticStiffness},
        _adhesiveStiffness{adhesiveStiffness},
        _frictionStiffness{frictionStiffness},
        _normalViscosity{normalViscosity},
        _frictionViscosity{frictionViscosity},
        _rollingViscosity{rollingViscosity},
        _torsionViscosity{torsionViscosity},
        _staticFrictionCoeff{staticFrictionCoeff},
        _dynamicFrictionCoeff{dynamicFrictionCoeff},
        _rollingFrictionCoeff{rollingFrictionCoeff},
        _torsionFrictionCoeff{torsionFrictionCoeff},
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
  explicit DEMFunctor(double cutoff)
      : DEMFunctor(cutoff, 50., 5, 5, 5e-5, 1e-3, 2.5e-6, 2.5e-6, 25, 100, 25, 25, nullptr) {
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
      : DEMFunctor(cutoff, 50., 5, 5, 5e-5, 1e-3, 2.5e-6, 2.5e-6, 25, 100, 25, 25, nullptr) {
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
   * @param torsionViscosity
   * @param staticFrictionCoeff
   * @param dynamicFrictionCoeff
   * @param rollingFrictionCoeff
   * @param torsionFrictionCoeff
   */
  explicit DEMFunctor(double cutoff, double elasticStiffness, double adhesiveStiffness, double frictionStiffness,
                      double normalViscosity, double frictionViscosity, double rollingViscosity,
                      double torsionViscosity, double staticFrictionCoeff, double dynamicFrictionCoeff,
                      double rollingFrictionCoeff, double torsionFrictionCoeff)
      : DEMFunctor(cutoff, elasticStiffness, adhesiveStiffness, frictionStiffness, normalViscosity, frictionViscosity,
                   rollingViscosity, torsionViscosity, staticFrictionCoeff, dynamicFrictionCoeff, rollingFrictionCoeff,
                   torsionFrictionCoeff, nullptr) {
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
   * @param torsionViscosity
   * @param staticFrictionCoeff
   * @param dynamicFrictionCoeff
   * @param rollingFrictionCoeff
   * @param torsionFrictionCoeff
   * @param particlePropertiesLibrary
   */
  explicit DEMFunctor(double cutoff, double elasticStiffness, double adhesiveStiffness, double frictionStiffness,
                      double normalViscosity, double frictionViscosity, double rollingViscosity,
                      double torsionViscosity, double staticFrictionCoeff, double dynamicFrictionCoeff,
                      double rollingFrictionCoeff, double torsionFrictionCoeff,
                      ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : DEMFunctor(cutoff, elasticStiffness, adhesiveStiffness, frictionStiffness, normalViscosity, frictionViscosity,
                   rollingViscosity, torsionViscosity, staticFrictionCoeff, dynamicFrictionCoeff, rollingFrictionCoeff,
                   torsionFrictionCoeff, nullptr) {
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
    using namespace autopas::utils::ArrayMath;

    if (i.isDummy() or j.isDummy()) {
      return;
    }

    const auto threadnum = autopas::autopas_get_thread_num();

    if constexpr (countFLOPs) {
      ++_aosThreadDataFLOPs[threadnum].numDistCalls;
    }

    // Compute necessary values and check for distance and contact (3 + 6 = 9 FLOPS)
    const double cutoff = _cutoff;
    const std::array<double, 3> displacement = i.getR() - j.getR();
    const double dist = autopas::utils::ArrayMath::L2Norm(displacement);

    // Retrieve particle-specific parameters, calculate generally necessary values (2 FlOPS)
    const auto [radiusI, radiusJ] = retrieveRadii(i, j);
    const double overlap = radiusI + radiusJ - dist;
    const bool overlapIsPositive = overlap > 0;

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numContacts += (overlapIsPositive ? 1 : 0);
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += (overlapIsPositive and newton3 ? 1 : 0);
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += (overlapIsPositive and not newton3 ? 1 : 0);
    }

    if (dist > cutoff || overlap <= 0) return;  // VdW deactivated

    // (3 + 2 + 2 + 3 + 3 + 5 = 18 FLOPS)
    assert(dist != 0. && "Distance is zero, division by zero detected!");
    const std::array<double, 3> normalUnit = displacement / dist;
    const double radiusIReduced = radiusI - overlap / 2.;
    const double radiusJReduced = radiusJ - overlap / 2.;
    assert((radiusIReduced + radiusJReduced) != 0. && "Distance is zero, division by zero detected!");
    const double radiusReduced = radiusIReduced * radiusJReduced / (radiusIReduced + radiusJReduced);
    const std::array<double, 3> relVel = i.getV() - j.getV();
    const double relVelDotNormalUnit = dot(normalUnit, relVel);

    // Compute Forces
    // Compute normal forces
    const double normalContactFMag = _elasticStiffness * overlap - _normalViscosity * relVelDotNormalUnit;  // 3 FLOPS
    const double normalFMag = normalContactFMag;
    const std::array<double, 3> normalF = mulScalar(normalUnit, normalFMag);  // 3 FLOPS

    // Compute tangential force
    const std::array<double, 3> tanRelVel = relVel + cross(normalUnit * radiusIReduced, i.getAngularVel()) +
                                            cross(normalUnit * radiusJReduced, j.getAngularVel());          // 30 FLOPS
    const std::array<double, 3> normalRelVel = normalUnit * relVelDotNormalUnit;                            // 3 FLOPS
    const std::array<double, 3> tanVel = tanRelVel - normalRelVel;                                          // 3 FLOPS
    const double coulombLimit = _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);  // 3 FLOPS
    std::array<double, 3> tanF = tanVel * (-_frictionViscosity);                                            // 3 FLOPS
    const double tanFMag = L2Norm(tanF);                                                                    // 6 FLOPS

    if (tanFMag > coulombLimit) {
      // 3 + 3 + 3 = 9 FLOPS
      ++_aosThreadDataFLOPs[threadnum].numInnerIfTanFCalls;
      const std::array<double, 3> tanFUnit = tanF / tanFMag;
      const double scale = _dynamicFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);
      tanF = tanFUnit * scale;
    }

    // Compute total force
    const std::array<double, 3> totalF = normalF + tanF;  // 3 FLOPS

    // Apply forces
    i.addF(totalF);  // 3 FLOPS
    if (newton3) {
      j.subF(totalF);  // 3 FLOPS
    }

    // Compute Torques
    // Compute frictional torque
    const std::array<double, 3> frictionQI = cross(normalUnit * (-radiusIReduced), tanF);  // 3 + 9 = 12 FLOPS

    // Compute rolling torque
    const std::array<double, 3> rollingRelVel =
        (cross(normalUnit, i.getAngularVel()) - cross(normalUnit, j.getAngularVel())) *
        (-radiusReduced);                                                   // 9 + 9 + 3 + 3 = 24 FLOPS
    std::array<double, 3> rollingF = rollingRelVel * (-_rollingViscosity);  // 3 FLOPS
    const double rollingFMag = L2Norm(rollingF);                            // 6 FLOPS

    if (rollingFMag > coulombLimit) {  // 3 + 3 + 3 = 9 FLOPS
      ++_aosThreadDataFLOPs[threadnum].numInnerIfRollingQCalls;
      const std::array<double, 3> rollingFUnit = rollingF / rollingFMag;
      const double scale = _rollingFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);
      rollingF = rollingFUnit * scale;
    }
    const std::array<double, 3> rollingQI = cross(normalUnit * radiusReduced, rollingF);  // 3 + 9 = 12 FLOPS

    // Compute torsional torque
    const std::array<double, 3> torsionRelVel =
        normalUnit * (dot(normalUnit, i.getAngularVel()) - dot(normalUnit, j.getAngularVel())) *
        radiusReduced;                                                      // 5 + 5 + 1 + 1 + 3 = 15 FLOPS
    std::array<double, 3> torsionF = torsionRelVel * (-_torsionViscosity);  // 3 FLOPS
    const double torsionFMag = L2Norm(torsionF);                            // 6 FLOPS

    if (torsionFMag > coulombLimit) {  // 3 + 3 + 3 = 9 FLOPS
      ++_aosThreadDataFLOPs[threadnum].numInnerIfTorsionQCalls;
      const std::array<double, 3> torsionFUnit = torsionF / torsionFMag;
      const double scale = _torsionFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);
      torsionF = torsionFUnit * scale;
    }
    const std::array<double, 3> torsionQI = torsionF * radiusReduced;  // 3 = 3 FLOPS

    // Apply torques
    i.addTorque(frictionQI + rollingQI + torsionQI);  // 9 FLOPS
    if (newton3) {
      j.addTorque((frictionQI * (radiusJReduced / radiusIReduced)) - rollingQI - torsionQI);  // 10 FLOPS
    }
  }

  /**
   * @copydoc autopas::Functor::SoAFunctorSingle()
   */
  void SoAFunctorSingle(autopas::SoAView<SoAArraysType> soa, bool newton3) final {
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

    // Count FLOPS variables
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numContactsSum = 0;
    size_t numInnerIfTanFCallsSum = 0;
    size_t numInnerIfRollingQCallsSum = 0;
    size_t numInnerIfTorsionQCallsSum = 0;

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
      }

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, qXacc, qYacc, qZacc, numDistanceCalculationSum,                \
                               numKernelCallsN3Sum, numKernelCallsNoN3Sum, numContactsSum, numInnerIfTanFCallsSum, \
                               numInnerIfRollingQCallsSum, numInnerIfTorsionQCallsSum)
      for (unsigned int j = i + 1; j < soa.size(); ++j) {
        // Compute necessary values for computations of forces
        const auto ownedStateJ = ownedStatePtr[j];

        const SoAFloatPrecision drx = xptr[i] - xptr[j];
        const SoAFloatPrecision dry = yptr[i] - yptr[j];
        const SoAFloatPrecision drz = zptr[i] - zptr[j];

        double radiusJ = constRadius;

        if (useMixing) {
          radiusJ = _PPLibrary->getRadius(typeptr[j]);
        }

        const SoAFloatPrecision distSquared = drx * drx + dry * dry + drz * drz;
        const SoAFloatPrecision dist = std::sqrt(distSquared);
        SoAFloatPrecision overlap = (radiusI + radiusJ) - dist;

        // Mask away if distance is too large or overlap is non-positive or any particle is a dummy
        const bool cutOffMask = dist <= _cutoff and ownedStateJ != autopas::OwnershipState::dummy;
        const bool overlapIsPositive = overlap > 0;

        if constexpr (countFLOPs) {
          numContactsSum += overlapIsPositive ? 1 : 0;
          numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
          numKernelCallsN3Sum += (cutOffMask and overlapIsPositive ? 1 : 0);
        }

        if (not(cutOffMask and overlapIsPositive)) {
          continue;  // VdW deactivated
        }

        const SoAFloatPrecision invDist = 1. / dist;
        const SoAFloatPrecision normalUnitX = drx * invDist;
        const SoAFloatPrecision normalUnitY = dry * invDist;
        const SoAFloatPrecision normalUnitZ = drz * invDist;

        const SoAFloatPrecision radiusIReduced = radiusI - overlap / 2.;
        const SoAFloatPrecision radiusJReduced = radiusJ - overlap / 2.;

        const SoAFloatPrecision relVelX = vxptr[i] - vxptr[j];
        const SoAFloatPrecision relVelY = vyptr[i] - vyptr[j];
        const SoAFloatPrecision relVelZ = vzptr[i] - vzptr[j];

        const SoAFloatPrecision relVelDotNormalUnit =
            relVelX * normalUnitX + relVelY * normalUnitY + relVelZ * normalUnitZ;

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

        const SoAFloatPrecision tanRelVelX = relVelX + relVelAngularIX + relVelAngularJX;
        const SoAFloatPrecision tanRelVelY = relVelY + relVelAngularIY + relVelAngularJY;
        const SoAFloatPrecision tanRelVelZ = relVelY + relVelAngularIZ + relVelAngularJZ;

        const SoAFloatPrecision tanRelVelTotalX = tanRelVelX - relVelDotNormalUnit * normalUnitX;
        const SoAFloatPrecision tanRelVelTotalY = tanRelVelY - relVelDotNormalUnit * normalUnitY;
        const SoAFloatPrecision tanRelVelTotalZ = tanRelVelZ - relVelDotNormalUnit * normalUnitZ;

        // Compute normal force as sum of contact force and long-range force (VdW).
        const SoAFloatPrecision normalContactFMag =
            _elasticStiffness * overlap - _normalViscosity * relVelDotNormalUnit;
        const SoAFloatPrecision normalFMag = normalContactFMag;  // VdW deactivated

        const SoAFloatPrecision normalFX = normalFMag * normalUnitX;
        const SoAFloatPrecision normalFY = normalFMag * normalUnitY;
        const SoAFloatPrecision normalFZ = normalFMag * normalUnitZ;

        // Compute tangential force
        SoAFloatPrecision tanFX = -_frictionViscosity * tanRelVelTotalX;
        SoAFloatPrecision tanFY = -_frictionViscosity * tanRelVelTotalY;
        SoAFloatPrecision tanFZ = -_frictionViscosity * tanRelVelTotalZ;

        const SoAFloatPrecision tanFMag = std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ);
        const SoAFloatPrecision coulombLimit =
            _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);

        if (tanFMag > coulombLimit) {
          if constexpr (countFLOPs) {
            ++numInnerIfTanFCallsSum;
          }
          const SoAFloatPrecision scale =
              _dynamicFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap) / tanFMag;
          tanFX *= scale;
          tanFY *= scale;
          tanFZ *= scale;
        }

        // Compute total force
        const SoAFloatPrecision totalFX = (normalFX + tanFX);
        const SoAFloatPrecision totalFY = (normalFY + tanFY);
        const SoAFloatPrecision totalFZ = (normalFZ + tanFZ);

        // Apply total force
        fxacc += totalFX;
        fyacc += totalFY;
        fzacc += totalFZ;

        // newton3
        fxptr[j] -= totalFX;
        fyptr[j] -= totalFY;
        fzptr[j] -= totalFZ;

        // Compute torques
        // Compute frictional torque
        const SoAFloatPrecision frictionQIX = -radiusIReduced * (normalUnitY * tanFZ - normalUnitZ * tanFY);
        const SoAFloatPrecision frictionQIY = -radiusIReduced * (normalUnitZ * tanFX - normalUnitX * tanFZ);
        const SoAFloatPrecision frictionQIZ = -radiusIReduced * (normalUnitX * tanFY - normalUnitY * tanFX);

        // Compute rolling torque
        const SoAFloatPrecision radiusReduced = radiusIReduced * radiusJReduced / (radiusIReduced + radiusJReduced);
        const SoAFloatPrecision rollingRelVelX =
            -radiusReduced * (normalUnitY * (angularVelZptr[i] - angularVelZptr[j]) -
                              normalUnitZ * (angularVelYptr[i] - angularVelYptr[j]));
        const SoAFloatPrecision rollingRelVelY =
            -radiusReduced * (normalUnitZ * (angularVelXptr[i] - angularVelXptr[j]) -
                              normalUnitX * (angularVelZptr[i] - angularVelZptr[j]));
        const SoAFloatPrecision rollingRelVelZ =
            -radiusReduced * (normalUnitX * (angularVelYptr[i] - angularVelYptr[j]) -
                              normalUnitY * (angularVelXptr[i] - angularVelXptr[j]));

        SoAFloatPrecision rollingFX = rollingRelVelX * (-_rollingViscosity);
        SoAFloatPrecision rollingFY = rollingRelVelY * (-_rollingViscosity);
        SoAFloatPrecision rollingFZ = rollingRelVelZ * (-_rollingViscosity);

        const SoAFloatPrecision rollingFMag =
            std::sqrt(rollingFX * rollingFX + rollingFY * rollingFY + rollingFZ * rollingFZ);

        if (rollingFMag > coulombLimit) {
          if constexpr (countFLOPs) {
            ++numInnerIfRollingQCallsSum;
          }
          const SoAFloatPrecision scale =
              _rollingFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap) / rollingFMag;
          rollingFX *= scale;
          rollingFY *= scale;
          rollingFZ *= scale;
        }

        const SoAFloatPrecision rollingQIX = radiusReduced * (normalUnitY * rollingFZ - normalUnitZ * rollingFY);
        const SoAFloatPrecision rollingQIY = radiusReduced * (normalUnitZ * rollingFX - normalUnitX * rollingFZ);
        const SoAFloatPrecision rollingQIZ = radiusReduced * (normalUnitX * rollingFY - normalUnitY * rollingFX);

        // Compute torsion torque
        const SoAFloatPrecision torsionRelVelScalar =
            radiusReduced * ((angularVelXptr[i] - angularVelXptr[j]) * normalUnitX +
                             (angularVelYptr[i] - angularVelYptr[j]) * normalUnitY +
                             (angularVelZptr[i] - angularVelZptr[j]) * normalUnitZ);
        const SoAFloatPrecision torsionRelVelX = torsionRelVelScalar * normalUnitX;
        const SoAFloatPrecision torsionRelVelY = torsionRelVelScalar * normalUnitY;
        const SoAFloatPrecision torsionRelVelZ = torsionRelVelScalar * normalUnitZ;

        SoAFloatPrecision torsionFX = torsionRelVelX * (-_torsionViscosity);
        SoAFloatPrecision torsionFY = torsionRelVelY * (-_torsionViscosity);
        SoAFloatPrecision torsionFZ = torsionRelVelZ * (-_torsionViscosity);

        const SoAFloatPrecision torsionFMag =
            std::sqrt(torsionFX * torsionFX + torsionFY * torsionFY + torsionFZ * torsionFZ);

        if (torsionFMag > coulombLimit) {
          if constexpr (countFLOPs) {
            ++numInnerIfTorsionQCallsSum;
          }
          const SoAFloatPrecision scale =
              _torsionFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap) / torsionFMag;
          torsionFX *= scale;
          torsionFY *= scale;
          torsionFZ *= scale;
        }

        const SoAFloatPrecision torsionQIX = radiusReduced * torsionFX;
        const SoAFloatPrecision torsionQIY = radiusReduced * torsionFY;
        const SoAFloatPrecision torsionQIZ = radiusReduced * torsionFZ;

        // Apply torques
        qXacc += (frictionQIX + rollingQIX + torsionQIX);
        qYacc += (frictionQIY + rollingQIY + torsionQIY);
        qZacc += (frictionQIZ + rollingQIZ + torsionQIZ);

        // newton3
        qXptr[j] += ((radiusJReduced / radiusIReduced) * frictionQIX - rollingQIX - torsionQIX);
        qYptr[j] += ((radiusJReduced / radiusIReduced) * frictionQIY - rollingQIY - torsionQIY);
        qZptr[j] += ((radiusJReduced / radiusIReduced) * frictionQIZ - rollingQIZ - torsionQIZ);

      }  // end of j loop

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;

      qXptr[i] += qXacc;
      qYptr[i] += qYacc;
      qZptr[i] += qZacc;
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numContacts += numContactsSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfTanFCalls += numInnerIfTanFCallsSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfRollingQCalls += numInnerIfRollingQCallsSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfTorsionQCalls += numInnerIfTorsionQCallsSum;
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
    return std::array<typename Particle::AttributeNames, 18>{
        Particle::AttributeNames::id,          Particle::AttributeNames::posX,
        Particle::AttributeNames::posY,        Particle::AttributeNames::posZ,
        Particle::AttributeNames::velocityX,   Particle::AttributeNames::velocityY,
        Particle::AttributeNames::velocityZ,   Particle::AttributeNames::forceX,
        Particle::AttributeNames::forceY,      Particle::AttributeNames::forceZ,
        Particle::AttributeNames::angularVelX, Particle::AttributeNames::angularVelY,
        Particle::AttributeNames::angularVelZ, Particle::AttributeNames::torqueX,
        Particle::AttributeNames::torqueY,     Particle::AttributeNames::torqueZ,
        Particle::AttributeNames::typeId,      Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 18>{
        Particle::AttributeNames::id,          Particle::AttributeNames::posX,
        Particle::AttributeNames::posY,        Particle::AttributeNames::posZ,
        Particle::AttributeNames::velocityX,   Particle::AttributeNames::velocityY,
        Particle::AttributeNames::velocityZ,   Particle::AttributeNames::forceX,
        Particle::AttributeNames::forceY,      Particle::AttributeNames::forceZ,
        Particle::AttributeNames::angularVelX, Particle::AttributeNames::angularVelY,
        Particle::AttributeNames::angularVelZ, Particle::AttributeNames::torqueX,
        Particle::AttributeNames::torqueY,     Particle::AttributeNames::torqueZ,
        Particle::AttributeNames::typeId,      Particle::AttributeNames::ownershipState};
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
    /**
     * FLOP count:
     * DistCall: 11
     * KernelNoN3: 24 (NormalF) + 48 (TanF) + 6 (ApplyF) + 12 (FrictionQ) + 45 (RollingQ) + 27 (TorsionQ) + 9 (ApplyQ) = 171
     * KernelN3: KernelNoN3 + 3 (ApplyF) + 10 (ApplyQ) = 184
     * InnerIfTanFCall: 9
     * InnerIfRollingQCall: 9
     * InnerIfTorsionQCall: 9
     */
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

      // Inner-If calls
      const size_t numInnerIfTanFCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numInnerIfTanFCalls; });
      const size_t numInnerIfRollingQCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numInnerIfRollingQCalls; });
      const size_t numInnerIfTorsionQCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numInnerIfTorsionQCalls; });

      constexpr size_t numFLOPsPerDistanceCall = 11;
      constexpr size_t numFLOPsPerNoN3KernelCall = 171;
      constexpr size_t numFLOPsPerN3KernelCall = 184;
      constexpr size_t numFLOPsPerInnerIfTanFCall = 9;
      constexpr size_t numFLOPsPerInnerIfRollingQCall = 9;
      constexpr size_t numFLOPsPerInnerIfTorsionQCall = 9;

      return numDistCallsAcc * numFLOPsPerDistanceCall + numKernelCallsN3Acc * numFLOPsPerN3KernelCall +
             numKernelCallsNoN3Acc * numFLOPsPerNoN3KernelCall + numInnerIfTanFCallsAcc * numFLOPsPerInnerIfTanFCall +
             numInnerIfRollingQCallsAcc * numFLOPsPerInnerIfRollingQCall + numInnerIfTorsionQCallsAcc * numFLOPsPerInnerIfTorsionQCall;
    } else {
      // This is needed because this function still gets called with FLOP logging disabled, just nothing is done with it
      return std::numeric_limits<size_t>::max();
    }
  }

  [[nodiscard]] double getHitRate() const override {
    if constexpr (countFLOPs) {  // TODO: use "HitRate" only as variable to count number of contacts for now
      const size_t numContacts = std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                                                 [](size_t sum, const auto &data) { return sum + data.numContacts; });

      return static_cast<double>(numContacts);
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

    // Variables for counting FLOPS
    size_t numDistanceCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numContactsSum = 0;
    size_t numInnerIfTanFCallsSum = 0;
    size_t numInnerIfRollingQCallsSum = 0;
    size_t numInnerIfTorsionQCallsSum = 0;

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
      if constexpr (useMixing) {
        radiusI = _PPLibrary->getRadius(typeptr1[i]);
      }

      // Loop over Particles in soa2
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, qXacc, qYacc, qZacc, numDistanceCalculationSum,                \
                               numKernelCallsN3Sum, numKernelCallsNoN3Sum, numContactsSum, numInnerIfTanFCallsSum, \
                               numInnerIfRollingQCallsSum, numInnerIfTorsionQCallsSum)
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
        }
        const SoAFloatPrecision overlap = (radiusI + radiusJ) - dist;

        // Mask away if distance is too large or overlap is non-positive or any particle is a dummy
        const bool cutOffMask = dist <= _cutoff and ownedStateJ != autopas::OwnershipState::dummy;
        const bool overlapIsPositive = overlap > 0;

        if constexpr (countFLOPs) {
          numContactsSum += overlapIsPositive ? 1 : 0;
          numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
          if constexpr (newton3) {
            numKernelCallsN3Sum += (cutOffMask and overlapIsPositive ? 1 : 0);
          } else {
            numKernelCallsNoN3Sum += (cutOffMask and overlapIsPositive ? 1 : 0);
          }
        }

        if (not(cutOffMask and overlapIsPositive)) {
          continue;  // VdW deactivated
        }

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

        const SoAFloatPrecision relVelX = vxptr1[i] - vxptr2[j];
        const SoAFloatPrecision relVelY = vyptr1[i] - vyptr2[j];
        const SoAFloatPrecision relVelZ = vzptr1[i] - vzptr2[j];

        const SoAFloatPrecision relVelDotNormalUnit =
            relVelX * normalUnitX + relVelY * normalUnitY + relVelZ * normalUnitZ;

        const SoAFloatPrecision tanRelVelX = relVelX + relVelAngularIX + relVelAngularJX;
        const SoAFloatPrecision tanRelVelY = relVelY + relVelAngularIY + relVelAngularJY;
        const SoAFloatPrecision tanRelVelZ = relVelZ + relVelAngularIZ + relVelAngularJZ;

        const SoAFloatPrecision tanRelVelTotalX = tanRelVelX - relVelDotNormalUnit * normalUnitX;
        const SoAFloatPrecision tanRelVelTotalY = tanRelVelY - relVelDotNormalUnit * normalUnitY;
        const SoAFloatPrecision tanRelVelTotalZ = tanRelVelZ - relVelDotNormalUnit * normalUnitZ;

        // Compute normal force
        const SoAFloatPrecision normalContactFMag =
            _elasticStiffness * overlap - _normalViscosity * relVelDotNormalUnit;
        const SoAFloatPrecision normalFMag = normalContactFMag;  // VdW deactivated

        const SoAFloatPrecision normalFX = normalFMag * normalUnitX;
        const SoAFloatPrecision normalFY = normalFMag * normalUnitY;
        const SoAFloatPrecision normalFZ = normalFMag * normalUnitZ;

        // Compute tangential force
        SoAFloatPrecision tanFX = -_frictionViscosity * tanRelVelTotalX;
        SoAFloatPrecision tanFY = -_frictionViscosity * tanRelVelTotalY;
        SoAFloatPrecision tanFZ = -_frictionViscosity * tanRelVelTotalZ;

        const SoAFloatPrecision tanFMag = std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ);
        const SoAFloatPrecision coulombLimit =
            _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);

        if (tanFMag > coulombLimit) {
          if constexpr (countFLOPs) {
            ++numInnerIfTanFCallsSum;
          }
          const SoAFloatPrecision scale =
              _dynamicFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap) / tanFMag;
          tanFX *= scale;
          tanFY *= scale;
          tanFZ *= scale;
        }

        // Compute total force
        const SoAFloatPrecision totalFX = (normalFX + tanFX);
        const SoAFloatPrecision totalFY = (normalFY + tanFY);
        const SoAFloatPrecision totalFZ = (normalFZ + tanFZ);

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
        const SoAFloatPrecision frictionQIX = -radiusIReduced * (normalUnitY * tanFZ - normalUnitZ * tanFY);
        const SoAFloatPrecision frictionQIY = -radiusIReduced * (normalUnitZ * tanFX - normalUnitX * tanFZ);
        const SoAFloatPrecision frictionQIZ = -radiusIReduced * (normalUnitX * tanFY - normalUnitY * tanFX);

        // Compute rolling torque
        const SoAFloatPrecision radiusReduced = radiusIReduced * radiusJReduced / (radiusIReduced + radiusJReduced);
        const SoAFloatPrecision rollingRelVelX =
            -radiusReduced * (normalUnitY * (angularVelZptr1[i] - angularVelZptr2[j]) -
                              normalUnitZ * (angularVelYptr1[i] - angularVelYptr2[j]));
        const SoAFloatPrecision rollingRelVelY =
            -radiusReduced * (normalUnitZ * (angularVelXptr1[i] - angularVelXptr2[j]) -
                              normalUnitX * (angularVelZptr1[i] - angularVelZptr2[j]));
        const SoAFloatPrecision rollingRelVelZ =
            -radiusReduced * (normalUnitX * (angularVelYptr1[i] - angularVelYptr2[j]) -
                              normalUnitY * (angularVelXptr1[i] - angularVelXptr2[j]));

        SoAFloatPrecision rollingFX = rollingRelVelX * (-_rollingViscosity);
        SoAFloatPrecision rollingFY = rollingRelVelY * (-_rollingViscosity);
        SoAFloatPrecision rollingFZ = rollingRelVelZ * (-_rollingViscosity);

        const SoAFloatPrecision rollingFMag =
            std::sqrt(rollingFX * rollingFX + rollingFY * rollingFY + rollingFZ * rollingFZ);

        if (rollingFMag > coulombLimit) {
          if constexpr (countFLOPs) {
            ++numInnerIfRollingQCallsSum;
          }
          const SoAFloatPrecision scale =
              _rollingFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap) / rollingFMag;
          rollingFX *= scale;
          rollingFY *= scale;
          rollingFZ *= scale;
        }

        const SoAFloatPrecision rollingQIX = radiusReduced * (normalUnitY * rollingFZ - normalUnitZ * rollingFY);
        const SoAFloatPrecision rollingQIY = radiusReduced * (normalUnitZ * rollingFX - normalUnitX * rollingFZ);
        const SoAFloatPrecision rollingQIZ = radiusReduced * (normalUnitX * rollingFY - normalUnitY * rollingFX);

        // Compute torsion torque
        const SoAFloatPrecision torsionRelVelScalar =
            radiusReduced * ((angularVelXptr1[i] - angularVelXptr2[j]) * normalUnitX +
                             (angularVelYptr1[i] - angularVelYptr2[j]) * normalUnitY +
                             (angularVelZptr1[i] - angularVelZptr2[j]) * normalUnitZ);
        const SoAFloatPrecision torsionRelVelX = torsionRelVelScalar * normalUnitX;
        const SoAFloatPrecision torsionRelVelY = torsionRelVelScalar * normalUnitY;
        const SoAFloatPrecision torsionRelVelZ = torsionRelVelScalar * normalUnitZ;

        SoAFloatPrecision torsionFX = torsionRelVelX * (-_torsionViscosity);
        SoAFloatPrecision torsionFY = torsionRelVelY * (-_torsionViscosity);
        SoAFloatPrecision torsionFZ = torsionRelVelZ * (-_torsionViscosity);

        const SoAFloatPrecision torsionFMag =
            std::sqrt(torsionFX * torsionFX + torsionFY * torsionFY + torsionFZ * torsionFZ);

        if (torsionFMag > coulombLimit) {
          if constexpr (countFLOPs) {
            ++numInnerIfTorsionQCallsSum;
          }
          const SoAFloatPrecision scale =
              _torsionFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap) / torsionFMag;
          torsionFX *= scale;
          torsionFY *= scale;
          torsionFZ *= scale;
        }

        const SoAFloatPrecision torsionQIX = radiusReduced * torsionFX;
        const SoAFloatPrecision torsionQIY = radiusReduced * torsionFY;
        const SoAFloatPrecision torsionQIZ = radiusReduced * torsionFZ;

        // Apply torques
        qXacc += (frictionQIX + rollingQIX + torsionQIX);
        qYacc += (frictionQIY + rollingQIY + torsionQIY);
        qZacc += (frictionQIZ + rollingQIZ + torsionQIZ);

        if (newton3) {
          qXptr2[j] += ((radiusJReduced / radiusIReduced) * frictionQIX - rollingQIX - torsionQIX);
          qYptr2[j] += ((radiusJReduced / radiusIReduced) * frictionQIY - rollingQIY - torsionQIY);
          qZptr2[j] += ((radiusJReduced / radiusIReduced) * frictionQIZ - rollingQIZ - torsionQIZ);
        }
      }  // end of j loop

      fxptr1[i] += fxacc;
      fyptr1[i] += fyacc;
      fzptr1[i] += fzacc;

      qXptr1[i] += qXacc;
      qYptr1[i] += qYacc;
      qZptr1[i] += qZacc;
    }

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numContacts += numContactsSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfTanFCalls += numInnerIfTanFCallsSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfRollingQCalls += numInnerIfRollingQCallsSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfTorsionQCalls += numInnerIfTorsionQCallsSum;
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
      numContacts = 0;
      numKernelCallsNoN3 = 0;
      numKernelCallsN3 = 0;
      numDistCalls = 0;
      numInnerIfTanFCalls = 0;
      numInnerIfRollingQCalls = 0;
      numInnerIfTorsionQCalls = 0;
      // numGlobalCalcsNoN3 = 0;
      // numGlobalCalcsN3 = 0;
    }

    /**
     * Number of contacts.
     */
    size_t numContacts = 0;

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
     * Number of inner if tanF calls.
     */
    size_t numInnerIfTanFCalls = 0;

    /**
     * Number of inner if rollingQ calls.
     */
    size_t numInnerIfRollingQCalls = 0;

    /**
     * Number of inner if torsionQ calls.
     */
    size_t numInnerIfTorsionQCalls = 0;

    /**
     * Counter for the number of times the globals have been calculated. Excludes the special case that N3 is enabled
     * and we calculate globals for an owned-halo particle pair.
     */
    // size_t numGlobalCalcsN3 = 0;

    /**
     * Counter for the number of times the globals have been calculated. Excludes the special case that N3 is enabled
     * and we calculate globals for an owned-halo particle pair.
     */
    // size_t numGlobalCalcsNoN3 = 0;

   private:
    /**
     * dummy parameter to get the right size (64 bytes)
     */
    double __remainingTo64[(64 - 5 * sizeof(size_t)) / sizeof(size_t)];
  };

  // make sure of the size of AoSThreadDataGlobals and AoSThreadDataFLOPs
  static_assert(sizeof(AoSThreadDataGlobals) % 64 == 0, "AoSThreadDataGlobals has wrong size");
  // static_assert(sizeof(AoSThreadDataFLOPs) % 64 == 0, "AoSThreadDataFLOPs has wrong size");

  const double _cutoff;
  const double _elasticStiffness;
  const double _adhesiveStiffness;
  const double _frictionStiffness;
  const double _normalViscosity;
  const double _frictionViscosity;
  const double _rollingViscosity;
  const double _torsionViscosity;
  const double _staticFrictionCoeff;
  const double _dynamicFrictionCoeff;
  const double _rollingFrictionCoeff;
  const double _torsionFrictionCoeff;
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

  // ----------------------------------------Helper Methods------------------------------------------------------------

  std::tuple<double, double> retrieveRadii(const Particle &i, const Particle &j) const {
    if (!useMixing) {
      return {_radius, _radius};
    }
    double radiusI = _PPLibrary->getRadius(i.getTypeId());
    double radiusJ = _PPLibrary->getRadius(j.getTypeId());
    return {radiusI, radiusJ};
  }

  double computeNormalContactFMag(const double overlap, const double relVelNormalDotNormalUnit) {
    /**
     * If overlap <= 0: 0 FLOP
     * Else: 3 FLOPS
     */
    return overlap > 0 ? _elasticStiffness * overlap - _normalViscosity * relVelNormalDotNormalUnit : 0;
  }

  double computeNormalVdWFMag(const double overlap, const double dist, const double sigma, const double epsilon6,
                              const double cutoff) {
    /**
     * If overlap <= 0: 0 FLOP
     * Else: 19 FLOPS
     */
    if (overlap > 0) {
      return 0;
    }
    const double invSigma = 1. / sigma;                                             // 1 FLOP
    const double lj2 = (sigma * sigma) / (dist * dist);                             // 3 FLOPS
    const double lj7 = lj2 * lj2 * lj2 * (sigma / dist);                            // 4 FLOPS
    const double ljCutoff2 = (sigma * sigma) / (cutoff * cutoff);                   // 3 FLOPS
    const double ljCutoff7 = ljCutoff2 * ljCutoff2 * ljCutoff2 * (sigma / cutoff);  // 4 FLOPS
    return -epsilon6 * invSigma * (lj7 - ljCutoff7);                                // 4 FLOPS
  }

  std::array<double, 3> computeTangentialForce(const double overlap, const std::array<double, 3> relVel,
                                               const Particle &i, const Particle &j, const double radiusIReduced,
                                               const double radiusJReduced, const std::array<double, 3> &normalUnit,
                                               const double relVelNormalDotNormalUnit, const double normalContactFMag,
                                               const int threadnum) {
    /**
     * If overlap <= 0: 0 FLOP
     * Else:
     *  If tanFMag > coulombLimit: 60 FLOPS
     *  Else: 51 FLOPS
     */
    using namespace autopas::utils::ArrayMath::literals;
    if (overlap <= 0) {
      return {0, 0, 0};
    }
    const std::array<double, 3> tanRelVel =
        relVel + autopas::utils::ArrayMath::cross(normalUnit * radiusIReduced, i.getAngularVel()) +
        autopas::utils::ArrayMath::cross(normalUnit * radiusJReduced, j.getAngularVel());                   // 30 FLOPS
    const std::array<double, 3> normalRelVel = normalUnit * relVelNormalDotNormalUnit;                      // 3 FLOPS
    const std::array<double, 3> tanVel = tanRelVel - normalRelVel;                                          // 3 FLOPS
    const double coulombLimit = _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);  // 3 FLOPS

    const std::array<double, 3> tanF = tanVel * (-_frictionViscosity);  // 3 FLOPS
    const double tanFMag = autopas::utils::ArrayMath::L2Norm(tanF);     // 6 FLOPS
    if (tanFMag > coulombLimit) {
      // 3 + 3 + 3 = 9 FLOPS
      ++_aosThreadDataFLOPs[threadnum].numInnerIfTanFCalls;
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
    /**
     * If overlap <= 0: 0 FLOP
     * Else: 12 FLOPS
     */
    using namespace autopas::utils::ArrayMath::literals;
    if (overlap <= 0) {
      return {0, 0, 0};
    }
    return autopas::utils::ArrayMath::cross(normalUnit * (-radiusIReduced), tanF);  // 3 + 9 = 12 FLOPS
  }

  std::array<double, 3> computeRollingTorqueI(const double overlap, const double radiusReduced, const Particle &i,
                                              const Particle &j, const std::array<double, 3> &normalUnit,
                                              const double normalContactFMag, const int threadnum) {
    /**
     * If overlap <= 0: 0 FLOP
     * Else:
     *  If: rollingFMag > coulombLimit: 57 FLOPS
     *  Else: 48 FLOPS
     */
    using namespace autopas::utils::ArrayMath::literals;
    if (overlap <= 0) {
      return {0, 0, 0};
    }
    const std::array<double, 3> rollingRelVel = (autopas::utils::ArrayMath::cross(normalUnit, i.getAngularVel()) -
                                                 autopas::utils::ArrayMath::cross(normalUnit, j.getAngularVel())) *
                                                (-radiusReduced);            // 9 + 9 + 3 + 3 = 24 FLOPS
    std::array<double, 3> rollingF = rollingRelVel * (-_rollingViscosity);   // 3 FLOPS
    const double rollingFMag = autopas::utils::ArrayMath::L2Norm(rollingF);  // 6 FLOPS

    const double coulombLimit = _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);  // 3 FLOPS
    if (rollingFMag > coulombLimit) {  // 3 + 3 + 3 = 9 FLOPS
      ++_aosThreadDataFLOPs[threadnum].numInnerIfRollingQCalls;
      const std::array<double, 3> rollingFUnit = rollingF / rollingFMag;
      const double scale = _rollingFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);
      rollingF = rollingFUnit * scale;
    }
    return autopas::utils::ArrayMath::cross(normalUnit * radiusReduced, rollingF);  // 3 + 9 = 12 FLOPS
  }

  std::array<double, 3> computeTorsionTorqueI(const double overlap, const double radiusReduced, const Particle &i,
                                              const Particle &j, const std::array<double, 3> &normalUnit,
                                              const double normalContactFMag, const int threadnum) {
    /**
     * If overlap <= 0: 0 FLOP
     * Else:
     *  If: torsionFMag > coulombLimit: 34 FLOPS
     *  Else: 25 FLOPS
     */
    using namespace autopas::utils::ArrayMath::literals;
    using namespace autopas::utils::ArrayMath;
    if (overlap <= 0) {
      return {0, 0, 0};
    }
    const std::array<double, 3> torsionRelVel =
        normalUnit * (dot(normalUnit, i.getAngularVel()) - dot(normalUnit, j.getAngularVel())) *
        radiusReduced;                                                       // 3 + 3 + 1 + 3 = 10 FLOPS
    std::array<double, 3> torsionF = torsionRelVel * (-_torsionViscosity);   // 3 FLOPS
    const double torsionFMag = autopas::utils::ArrayMath::L2Norm(torsionF);  // 6 FLOPS

    const double coulombLimit = _staticFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);  // 3 FLOPS
    if (torsionFMag > coulombLimit) {  // 3 + 3 + 3 = 9 FLOPS
      ++_aosThreadDataFLOPs[threadnum].numInnerIfTorsionQCalls;
      const std::array<double, 3> torsionFUnit = torsionF / torsionFMag;
      const double scale = _torsionFrictionCoeff * (normalContactFMag + _adhesiveStiffness * overlap);
      torsionF = torsionFUnit * scale;
    }
    return torsionF * radiusReduced;  // 3 = 3 FLOPS
  }
};
}  // namespace demLib