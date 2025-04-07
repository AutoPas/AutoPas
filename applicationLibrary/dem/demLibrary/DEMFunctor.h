/**
 * @file DEMFunctor.h
 * @author Joon Kim
 * @date 27/03/2025
 */

#pragma once

#include "autopas/baseFunctors/PairwiseFunctor.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/StaticBoolSelector.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "DEMParameters.h"

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
    : public autopas::PairwiseFunctor<
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
   * @param demParameters Includes [elasticStiffness, normalViscosity, frictionViscosity, rollingViscosity, torsionViscosity, staticFrictionCoeff, dynamicFrictionCoeff, rollingFrictionCoeff, torsionFrictionCoeff, conductivity, heatGenerationFactor]
   * @note param dummy is unused, only there to make the signature different from the public constructor.
   */
  explicit DEMFunctor(double cutoff, DEMParameters &demParameters, void * /*dummy*/)
      : autopas::PairwiseFunctor<
            Particle, DEMFunctor<Particle, useMixing, useNewton3, calculateGloabls, countFLOPs, relevantForTuning>>(
            cutoff),
        _radius{cutoff / 5.},  // default initialization, can be changed by ParticlePropertiesLibrary
        _cutoff{cutoff},
        _elasticStiffness{demParameters.getElasticStiffness()},
        _normalViscosity{demParameters.getNormalViscosity()},
        _frictionViscosity{demParameters.getFrictionViscosity()},
        _rollingViscosity{demParameters.getRollingViscosity()},
        _torsionViscosity{demParameters.getTorsionViscosity()},
        _staticFrictionCoeff{demParameters.getStaticFrictionCoeff()},
        _dynamicFrictionCoeff{demParameters.getDynamicFrictionCoeff()},
        _rollingFrictionCoeff{demParameters.getRollingFrictionCoeff()},
        _torsionFrictionCoeff{demParameters.getTorsionFrictionCoeff()},
        _heatConductivity{demParameters.getHeatConductivity()},
        _heatGenerationFactor{demParameters.getHeatGenerationFactor()},
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
   * Minimal Constructor for Functor with mixing disabled.
   *
   * @note Only to be used with mixing == false;
   *
   * @param cutoff
   * @param demParameters
   */
  explicit DEMFunctor(double cutoff, DEMParameters &demParameters) : DEMFunctor(cutoff, demParameters, nullptr) {
    static_assert(not useMixing,
                  "Mixing without a ParticlePropertiesLibrary is not possible! Use a different constructor or set "
                  "mixing to false.");
  }

  /**
   * Constructor for Functor with mixing active. This functor takes a ParticlePropertiesLibrary to look up (mixed)
   * properties like radius.
   *
   * @param cutoff
   * @param demParameters
   * @param particlePropertiesLibrary
   */
  explicit DEMFunctor(double cutoff, DEMParameters &demParameters,
                      ParticlePropertiesLibrary<double, size_t> &particlePropertiesLibrary)
      : DEMFunctor(cutoff, demParameters, nullptr) {
    static_assert(useMixing,
                  "Not using Mixing but using a ParticlePropertiesLibrary is not allowed! Use a different constructor "
                  "or set mixing to true.");
    _PPLibrary = &particlePropertiesLibrary;
  }

  std::string getName() final { return "DEMFunctor"; }

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

    if (dist > cutoff) return;

    if constexpr (countFLOPs) {
      ++_aosThreadDataFLOPs[threadnum].numOverlapCalls;
    }

    // Retrieve particle-specific parameters, calculate generally necessary values (2 FlOPS)
    const auto [radiusI, radiusJ] = retrieveRadii(i, j);
    const double overlap = radiusI + radiusJ - dist;
    const bool overlapIsPositive = overlap > 0;

    if (not overlapIsPositive) return;  // VdW deactivated

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += (newton3 ? 1 : 0);
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += (not newton3 ? 1 : 0);
    }

    // (3 + 1 + 1 + 1 + 3 + 3 + 5 = 17 FLOPS)
    assert(dist != 0. && "Distance is zero, division by zero detected!");
    const std::array<double, 3> normalUnit = displacement / dist;
    const double overlapHalf = overlap / 2.;
    const double radiusIReduced = radiusI - overlapHalf;
    const double radiusJReduced = radiusJ - overlapHalf;
    assert((radiusIReduced + radiusJReduced) != 0. && "Distance is zero, division by zero detected!");
    const double radiusReduced = radiusIReduced * radiusJReduced / (radiusIReduced + radiusJReduced);
    const std::array<double, 3> relVel = i.getV() - j.getV();
    const double relVelDotNormalUnit = dot(normalUnit, relVel);

    // Compute Forces
    // Compute normal forces (3 + 3 = 6 FLOPS)
    const double normalContactFMag = _elasticStiffness * overlap - _normalViscosity * relVelDotNormalUnit;  // 3 FLOPS
    const double normalFMag = normalContactFMag;
    const std::array<double, 3> normalF = mulScalar(normalUnit, normalFMag);  // 3 FLOPS

    // Compute tangential force (30 + 3 + 3 + 1 + 3 + 6 = 46 FLOPS)
    const std::array<double, 3> tanRelVel = relVel + cross(normalUnit * radiusIReduced, i.getAngularVel()) +
                                            cross(normalUnit * radiusJReduced, j.getAngularVel());  // 30 FLOPS
    const std::array<double, 3> normalRelVel = normalUnit * relVelDotNormalUnit;                    // 3 FLOPS
    const std::array<double, 3> tanVel = tanRelVel - normalRelVel;                                  // 3 FLOPS
    const double coulombLimit = _staticFrictionCoeff * (normalContactFMag);                          // 1 FLOPS
    std::array<double, 3> tanF = tanVel * (-_frictionViscosity);                                     // 3 FLOPS
    const double tanFMag = L2Norm(tanF);                                                            // 6 FLOPS
    if (tanFMag > coulombLimit) {
      ++_aosThreadDataFLOPs[threadnum].numInnerIfTanFCalls;  // -> additional 5 FLOPS
      const double scale = _dynamicFrictionCoeff * (normalContactFMag) / tanFMag;  // 2 FLOPS
      tanF = tanF * scale;                                                                                  // 3 FLOPS
    }

    // Compute total force (3 FLOPS)
    const std::array<double, 3> totalF = normalF + tanF;  // 3 FLOPS

    // Apply forces (if newton3: 6 FLOPS, if not: 3 FLOPS)
    i.addF(totalF);  // 3 FLOPS
    if (newton3) {
      j.subF(totalF);  // 3 FLOPS
    }

    // Compute Torques
    // Compute frictional torque (12 FLOPS)
    const std::array<double, 3> frictionQI = cross(normalUnit * (-radiusIReduced), tanF);  // 3 + 9 = 12 FLOPS

    // Compute rolling torque (15 + 3 + 6 + 12 = 36 FLOPS)
    const std::array<double, 3> rollingRelVel =
        cross(normalUnit, (i.getAngularVel() - j.getAngularVel())) * (-radiusReduced);  // 15 FLOPS
    std::array<double, 3> rollingF = rollingRelVel * (-_rollingViscosity);               // 3 FLOPS
    const double rollingFMag = L2Norm(rollingF);                                        // 6 FLOPS
    if (rollingFMag > coulombLimit) {
      ++_aosThreadDataFLOPs[threadnum].numInnerIfRollingQCalls;  // -> additional 5 FLOPS
      const double scale =
          _rollingFrictionCoeff * (normalContactFMag) / rollingFMag;  // 2 FLOPS
      rollingF = rollingF * scale;
    }
    const std::array<double, 3> rollingQI = cross(normalUnit * radiusReduced, rollingF);  // 3 + 9 = 12 FLOPS

    // Compute torsional torque (11 + 3 + 6 + 3 = 23 FLOPS)
    const std::array<double, 3> torsionRelVel = normalUnit * dot(normalUnit, (i.getAngularVel() - j.getAngularVel())) *
                                                radiusReduced;             // 1 + 3 + 5 + 2 = 11 FLOPS
    std::array<double, 3> torsionF = torsionRelVel * (-_torsionViscosity);  // 3 FLOPS
    const double torsionFMag = L2Norm(torsionF);                           // 6 FLOPS
    if (torsionFMag > coulombLimit) {
      ++_aosThreadDataFLOPs[threadnum].numInnerIfTorsionQCalls;  // -> additional 5 FLOPS
      const double scale =
          _torsionFrictionCoeff * (normalContactFMag) / torsionFMag;  // 2 FLOPS
      torsionF = torsionF * scale;
    }
    const std::array<double, 3> torsionQI = torsionF * radiusReduced;  // 3 FLOPS

    // Apply torques (if newton3: 19 FLOPS, if not: 9 FLOPS)
    i.addTorque(frictionQI + rollingQI + torsionQI);  // 9 FLOPS
    if (newton3) {
      j.addTorque((frictionQI * (radiusJReduced / radiusIReduced)) - rollingQI - torsionQI);  // 10 FLOPS
    }

    // Heat Transfer (2 + 4 + 2 + 1 = 9 FLOPS, for newton3 additional 1 FLOP)
    const double geometricMeanRadius = std::sqrt(radiusI * radiusJ);
    const double conductance =
        2. * _heatConductivity * std::pow((3. * _elasticStiffness * overlap * geometricMeanRadius) / 4., 1. / 3.);
    const double heatFluxI = conductance * (j.getTemperature() - i.getTemperature());

    // Heat Generation (12 + 1 = 13 FLOPS, for newton3 additional 1 FLOP)
    const double heatFluxGenerated = _heatGenerationFactor * L2Norm(tanF) * L2Norm(tanVel);

    i.addHeatFlux(heatFluxI + heatFluxGenerated);
    if (newton3) {
      j.addHeatFlux(-heatFluxI + heatFluxGenerated);
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

    const auto *const __restrict temperaturePtr = soa.template begin<Particle::AttributeNames::temperature>();
    SoAFloatPrecision *const __restrict heatFluxPtr = soa.template begin<Particle::AttributeNames::heatFlux>();

    [[maybe_unused]] auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    // Initialize local variables for constants
    const SoAFloatPrecision cutoff = _cutoff;
    const SoAFloatPrecision constRadius = _radius;

    // Count FLOPS variables
    size_t numDistanceCalculationSum = 0;
    size_t numOverlapCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
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

      SoAFloatPrecision heatFluxAcc = 0.;

      const double radiusI = _PPLibrary->getRadius(typeptr[i]);

#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, qXacc, qYacc, qZacc, heatFluxAcc, numDistanceCalculationSum,     \
                               numOverlapCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, \
                               numInnerIfTanFCallsSum, numInnerIfRollingQCallsSum, numInnerIfTorsionQCallsSum)
      for (unsigned int j = i + 1; j < soa.size(); ++j) {
        // Compute necessary values for computations of forces
        const auto ownedStateJ = ownedStatePtr[j];

        const SoAFloatPrecision drx = xptr[i] - xptr[j];
        const SoAFloatPrecision dry = yptr[i] - yptr[j];
        const SoAFloatPrecision drz = zptr[i] - zptr[j];

        const double radiusJ = _PPLibrary->getRadius(typeptr[j]);

        const SoAFloatPrecision distSquared = drx * drx + dry * dry + drz * drz;
        const SoAFloatPrecision dist = std::sqrt(distSquared);
        SoAFloatPrecision overlap = (radiusI + radiusJ) - dist;

        // Mask away if distance is too large or overlap is non-positive or any particle is a dummy
        const bool cutOffMask = dist <= _cutoff and ownedStateJ != autopas::OwnershipState::dummy;
        const bool overlapIsPositive = overlap > 0;

        if constexpr (countFLOPs) {
          numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
          numOverlapCalculationSum += (cutOffMask ? 1 : 0);
          numKernelCallsN3Sum += (cutOffMask and overlapIsPositive ? 1 : 0);
        }

        if (not(cutOffMask and overlapIsPositive)) {
          continue;
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
        const SoAFloatPrecision tanRelVelZ = relVelZ + relVelAngularIZ + relVelAngularJZ;

        const SoAFloatPrecision tanRelVelTotalX = tanRelVelX - relVelDotNormalUnit * normalUnitX;
        const SoAFloatPrecision tanRelVelTotalY = tanRelVelY - relVelDotNormalUnit * normalUnitY;
        const SoAFloatPrecision tanRelVelTotalZ = tanRelVelZ - relVelDotNormalUnit * normalUnitZ;

        // Compute normal force as sum of contact force.
        const SoAFloatPrecision normalContactFMag = _elasticStiffness * overlap - _normalViscosity * relVelDotNormalUnit;

        const SoAFloatPrecision normalFX = normalContactFMag * normalUnitX;
        const SoAFloatPrecision normalFY = normalContactFMag * normalUnitY;
        const SoAFloatPrecision normalFZ = normalContactFMag * normalUnitZ;

        // Compute tangential force
        SoAFloatPrecision tanFX = -_frictionViscosity * tanRelVelTotalX;
        SoAFloatPrecision tanFY = -_frictionViscosity * tanRelVelTotalY;
        SoAFloatPrecision tanFZ = -_frictionViscosity * tanRelVelTotalZ;

        const SoAFloatPrecision tanFMag = std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ);
        const SoAFloatPrecision coulombLimit = _staticFrictionCoeff * normalContactFMag;
        if (tanFMag > coulombLimit) {
          if constexpr (countFLOPs) {
            ++numInnerIfTanFCallsSum;
          }
          const SoAFloatPrecision scale =
              _dynamicFrictionCoeff * normalContactFMag / tanFMag;
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
              _rollingFrictionCoeff * normalContactFMag / rollingFMag;
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
              _torsionFrictionCoeff * normalContactFMag / torsionFMag;
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

        // Compute Heat Transfer (2 + 3 + 2 + 2 = 9 FLOPS)
        const SoAFloatPrecision geometricMeanRadius = std::sqrt(radiusI * radiusJ);
        const SoAFloatPrecision conductance =
            2. * _heatConductivity * std::pow((3. * _elasticStiffness * overlap * geometricMeanRadius) / 4., 1. / 3.);
        const SoAFloatPrecision heatFluxI = conductance * (temperaturePtr[j] - temperaturePtr[i]);

        // Compute Heat Generation (12 + 1 = 13 FLOPS)
        const SoAFloatPrecision heatFluxGenerated =
            _heatGenerationFactor * std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ) *
            std::sqrt(tanRelVelTotalX * tanRelVelTotalX + tanRelVelTotalY * tanRelVelTotalY +
                      tanRelVelTotalZ * tanRelVelTotalZ);

        // Apply Heat Flux
        heatFluxAcc += (heatFluxI + heatFluxGenerated);
        heatFluxPtr[j] += (-heatFluxI + heatFluxGenerated);  // "newton 3"
      }  // end of j loop

      fxptr[i] += fxacc;
      fyptr[i] += fyacc;
      fzptr[i] += fzacc;

      qXptr[i] += qXacc;
      qYptr[i] += qYacc;
      qZptr[i] += qZacc;

      heatFluxPtr[i] += heatFluxAcc;
    }
    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numOverlapCalls += numOverlapCalculationSum;
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
   * @param radius
   */
  void setParticleProperties(SoAFloatPrecision radius) {
    _radius = radius;
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr()
   */
  constexpr static auto getNeededAttr() {
    return std::array<typename Particle::AttributeNames, 20>{
        Particle::AttributeNames::id,          Particle::AttributeNames::posX,
        Particle::AttributeNames::posY,        Particle::AttributeNames::posZ,
        Particle::AttributeNames::velocityX,   Particle::AttributeNames::velocityY,
        Particle::AttributeNames::velocityZ,   Particle::AttributeNames::forceX,
        Particle::AttributeNames::forceY,      Particle::AttributeNames::forceZ,
        Particle::AttributeNames::angularVelX, Particle::AttributeNames::angularVelY,
        Particle::AttributeNames::angularVelZ, Particle::AttributeNames::torqueX,
        Particle::AttributeNames::torqueY,     Particle::AttributeNames::torqueZ,
        Particle::AttributeNames::temperature, Particle::AttributeNames::heatFlux,
        Particle::AttributeNames::typeId,      Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getNeededAttr(std::false_type)
   */
  constexpr static auto getNeededAttr(std::false_type) {
    return std::array<typename Particle::AttributeNames, 20>{
        Particle::AttributeNames::id,          Particle::AttributeNames::posX,
        Particle::AttributeNames::posY,        Particle::AttributeNames::posZ,
        Particle::AttributeNames::velocityX,   Particle::AttributeNames::velocityY,
        Particle::AttributeNames::velocityZ,   Particle::AttributeNames::forceX,
        Particle::AttributeNames::forceY,      Particle::AttributeNames::forceZ,
        Particle::AttributeNames::angularVelX, Particle::AttributeNames::angularVelY,
        Particle::AttributeNames::angularVelZ, Particle::AttributeNames::torqueX,
        Particle::AttributeNames::torqueY,     Particle::AttributeNames::torqueZ,
        Particle::AttributeNames::temperature, Particle::AttributeNames::heatFlux,
        Particle::AttributeNames::typeId,      Particle::AttributeNames::ownershipState};
  }

  /**
   * @copydoc autopas::Functor::getComputedAttr()
   */
  constexpr static auto getComputedAttr() {
    return std::array<typename Particle::AttributeNames, 7>{
        Particle::AttributeNames::forceX,  Particle::AttributeNames::forceY,  Particle::AttributeNames::forceZ,
        Particle::AttributeNames::torqueX, Particle::AttributeNames::torqueY, Particle::AttributeNames::torqueZ,
        Particle::AttributeNames::heatFlux};
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
    if (calculateGloabls) {}
  }

  [[nodiscard]] size_t getNumFLOPs() const override {
    /**
     * FLOP count:
     * KernelNoN3: 155 + 9 + 13 = 177;
     * KernelN3: 155 + 13 + 10 + 14 = 192;
     */
    if constexpr (countFLOPs) {
      const size_t numDistCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numDistCalls; });
      const size_t numOverlapCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numOverlapCalls; });
      const size_t numKernelCallsN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsN3; });
      const size_t numKernelCallsNoN3Acc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numKernelCallsNoN3; });
      const size_t numInnerIfTanFCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numInnerIfTanFCalls; });
      const size_t numInnerIfRollingQCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numInnerIfRollingQCalls; });
      const size_t numInnerIfTorsionQCallsAcc =
          std::accumulate(_aosThreadDataFLOPs.begin(), _aosThreadDataFLOPs.end(), 0ul,
                          [](size_t sum, const auto &data) { return sum + data.numInnerIfTorsionQCalls; });

      constexpr size_t numFLOPsPerDistanceCall = 9;
      constexpr size_t numFLOPsPerOverlapCall = 2;
      constexpr size_t numFLOPsPerNoN3KernelCall = 177;
      constexpr size_t numFLOPsPerN3KernelCall = 192;
      constexpr size_t numFLOPsPerInnerIfTanF = 5;
      constexpr size_t numFLOPsPerInnerIfRollingQ = 5;
      constexpr size_t numFLOPsPerInnerIfTorsionQ = 5;

      return numDistCallsAcc * numFLOPsPerDistanceCall + numOverlapCallsAcc * numFLOPsPerOverlapCall +
             numKernelCallsN3Acc * numFLOPsPerN3KernelCall + numKernelCallsNoN3Acc * numFLOPsPerNoN3KernelCall +
             numInnerIfTanFCallsAcc * numFLOPsPerInnerIfTanF + numInnerIfRollingQCallsAcc * numFLOPsPerInnerIfRollingQ +
             numInnerIfTorsionQCallsAcc * numFLOPsPerInnerIfTorsionQ;
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

    const auto *const __restrict temperaturePtr1 = soa1.template begin<Particle::AttributeNames::temperature>();
    SoAFloatPrecision *const __restrict heatFluxPtr1 = soa1.template begin<Particle::AttributeNames::heatFlux>();

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

    const auto *const __restrict temperaturePtr2 = soa2.template begin<Particle::AttributeNames::temperature>();
    SoAFloatPrecision *const __restrict heatFluxPtr2 = soa2.template begin<Particle::AttributeNames::heatFlux>();

    auto *const __restrict typeptr1 = soa1.template begin<Particle::AttributeNames::typeId>();
    auto *const __restrict typeptr2 = soa2.template begin<Particle::AttributeNames::typeId>();

    // Initialize local variables for constants
    const SoAFloatPrecision cutoff = _cutoff;
    const SoAFloatPrecision radius = _radius;

    // Variables for counting FLOPS
    size_t numDistanceCalculationSum = 0;
    size_t numOverlapCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
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

      SoAFloatPrecision heatFluxAcc = 0.;

      SoAFloatPrecision radiusI = _PPLibrary->getRadius(typeptr1[i]);

      // Loop over Particles in soa2
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, qXacc, qYacc, qZacc, heatFluxAcc, numDistanceCalculationSum,     \
                               numOverlapCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, \
                               numInnerIfTanFCallsSum, numInnerIfRollingQCallsSum, numInnerIfTorsionQCallsSum)
      for (unsigned int j = 0; j < soa2.size(); ++j) {
        const auto ownedStateJ = ownedStatePtr2[j];

        const SoAFloatPrecision drx = xptr1[i] - xptr2[j];
        const SoAFloatPrecision dry = yptr1[i] - yptr2[j];
        const SoAFloatPrecision drz = zptr1[i] - zptr2[j];

        const SoAFloatPrecision distSquared = drx * drx + dry * dry + drz * drz;
        const SoAFloatPrecision dist = std::sqrt(distSquared);
        SoAFloatPrecision radiusJ = _PPLibrary->getRadius(typeptr2[j]);
        const SoAFloatPrecision overlap = (radiusI + radiusJ) - dist;

        // Mask away if distance is too large or overlap is non-positive or any particle is a dummy
        const bool cutOffMask = dist <= _cutoff and ownedStateJ != autopas::OwnershipState::dummy;
        const bool overlapIsPositive = overlap > 0;

        if constexpr (countFLOPs) {
          numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
          numOverlapCalculationSum += (cutOffMask ? 1 : 0);
          if constexpr (newton3) {
            numKernelCallsN3Sum += (cutOffMask and overlapIsPositive ? 1 : 0);
          } else {
            numKernelCallsNoN3Sum += (cutOffMask and overlapIsPositive ? 1 : 0);
          }
        }

        if (not(cutOffMask and overlapIsPositive)) {
          continue;
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
        const SoAFloatPrecision normalContactFMag = _elasticStiffness * overlap - _normalViscosity * relVelDotNormalUnit;

        const SoAFloatPrecision normalFX = normalContactFMag * normalUnitX;
        const SoAFloatPrecision normalFY = normalContactFMag * normalUnitY;
        const SoAFloatPrecision normalFZ = normalContactFMag * normalUnitZ;

        // Compute tangential force
        const SoAFloatPrecision tanFUnitMag = std::sqrt(
            tanRelVelTotalX * tanRelVelTotalX + tanRelVelTotalY * tanRelVelTotalY + tanRelVelTotalZ * tanRelVelTotalZ);
        const SoAFloatPrecision tanFScale =
            _dynamicFrictionCoeff * normalContactFMag / tanFUnitMag;

        SoAFloatPrecision tanFX = -_frictionViscosity * tanRelVelTotalX;
        SoAFloatPrecision tanFY = -_frictionViscosity * tanRelVelTotalY;
        SoAFloatPrecision tanFZ = -_frictionViscosity * tanRelVelTotalZ;

        const SoAFloatPrecision tanFMag = std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ);
        const SoAFloatPrecision coulombLimit = _staticFrictionCoeff * normalContactFMag;

        if (tanFMag > coulombLimit) {
          if constexpr (countFLOPs) {
            ++numInnerIfTanFCallsSum;
          }
          const SoAFloatPrecision scale =
              _dynamicFrictionCoeff * normalContactFMag / tanFMag;
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
              _rollingFrictionCoeff * normalContactFMag / rollingFMag;
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
              _torsionFrictionCoeff * normalContactFMag / torsionFMag;
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

        // Compute heat flux
        const SoAFloatPrecision geomMeanRadius = std::sqrt(radiusIReduced * radiusJReduced);
        const SoAFloatPrecision conductance =
            2. * _heatConductivity * std::pow((3. * _elasticStiffness * overlap * geomMeanRadius) / 4., 1. / 3.);
        const SoAFloatPrecision heatFluxI = conductance * (temperaturePtr2[j] - temperaturePtr1[i]);

        const SoAFloatPrecision heatFluxGenerated =
            _heatGenerationFactor * std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ) *
            std::sqrt(tanRelVelTotalX * tanRelVelTotalX + tanRelVelTotalY * tanRelVelTotalY +
                      tanRelVelTotalZ * tanRelVelTotalZ);

        // Apply heat flux
        heatFluxAcc += (heatFluxI);
        if (newton3) {
          heatFluxPtr2[j] += ((-heatFluxI));
        }
      }  // end of j loop

      fxptr1[i] += fxacc;
      fyptr1[i] += fyacc;
      fzptr1[i] += fzacc;

      qXptr1[i] += qXacc;
      qYptr1[i] += qYacc;
      qZptr1[i] += qZacc;

      heatFluxPtr1[i] += heatFluxAcc;
    }

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numOverlapCalls += numOverlapCalculationSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfTanFCalls += numInnerIfTanFCallsSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfRollingQCalls += numInnerIfRollingQCallsSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfTorsionQCalls += numInnerIfTorsionQCallsSum;
    }
  }

  template <bool newton3>
  void SoAFunctorVerletImpl(autopas::SoAView<SoAArraysType> soa, const size_t indexFirst,
                            const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList) {
    const auto *const __restrict xptr = soa.template begin<Particle::AttributeNames::posX>();
    const auto *const __restrict yptr = soa.template begin<Particle::AttributeNames::posY>();
    const auto *const __restrict zptr = soa.template begin<Particle::AttributeNames::posZ>();

    const auto *const __restrict vxptr = soa.template begin<Particle::AttributeNames::velocityX>();
    const auto *const __restrict vyptr = soa.template begin<Particle::AttributeNames::velocityY>();
    const auto *const __restrict vzptr = soa.template begin<Particle::AttributeNames::velocityZ>();

    const auto *const __restrict wxptr = soa.template begin<Particle::AttributeNames::angularVelX>();
    const auto *const __restrict wyptr = soa.template begin<Particle::AttributeNames::angularVelY>();
    const auto *const __restrict wzptr = soa.template begin<Particle::AttributeNames::angularVelZ>();

    auto *const __restrict fxptr = soa.template begin<Particle::AttributeNames::forceX>();
    auto *const __restrict fyptr = soa.template begin<Particle::AttributeNames::forceY>();
    auto *const __restrict fzptr = soa.template begin<Particle::AttributeNames::forceZ>();

    auto *const __restrict qxptr = soa.template begin<Particle::AttributeNames::torqueX>();
    auto *const __restrict qyptr = soa.template begin<Particle::AttributeNames::torqueY>();
    auto *const __restrict qzptr = soa.template begin<Particle::AttributeNames::torqueZ>();

    const auto *const __restrict temperaturePtr = soa.template begin<Particle::AttributeNames::temperature>();
    auto *const __restrict heatFluxPtr = soa.template begin<Particle::AttributeNames::heatFlux>();

    auto *const __restrict typeptr = soa.template begin<Particle::AttributeNames::typeId>();

    const auto *const __restrict ownedStatePtr = soa.template begin<Particle::AttributeNames::ownershipState>();

    // countFLOPs counters
    size_t numDistanceCalculationSum = 0;
    size_t numOverlapCalculationSum = 0;
    size_t numKernelCallsN3Sum = 0;
    size_t numKernelCallsNoN3Sum = 0;
    size_t numInnerIfTanFCallsSum = 0;
    size_t numInnerIfRollingQCallsSum = 0;
    size_t numInnerIfTorsionQCallsSum = 0;

    SoAFloatPrecision fxacc = 0;
    SoAFloatPrecision fyacc = 0;
    SoAFloatPrecision fzacc = 0;

    SoAFloatPrecision qxacc = 0;
    SoAFloatPrecision qyacc = 0;
    SoAFloatPrecision qzacc = 0;

    SoAFloatPrecision heatFluxAcc = 0;

    const size_t neighborListSize = neighborList.size();
    const size_t *const __restrict neighborListPtr = neighborList.data();

    // checks whether particle i is owned.
    const auto ownedStateI = ownedStatePtr[indexFirst];
    if (ownedStateI == autopas::OwnershipState::dummy) {
      return;
    }

    const auto threadnum = autopas::autopas_get_thread_num();

    constexpr size_t vecsize = 12;  // hyper-parameter

    size_t joff = 0;

    // retrieve the radius of particle i
    SoAFloatPrecision radiusI = _PPLibrary->getRadius(typeptr[indexFirst]);

    // if the size of the verlet list is larger than the given size vecsize,
    // we will use a vectorized version.
    if (neighborListSize >= vecsize) {
      alignas(64) std::array<SoAFloatPrecision, vecsize> xtmp, ytmp, ztmp, vxtmp, vytmp, vztmp, wxtmp, wytmp, wztmp,
          rtmp, temperaturetmp, xArr, yArr, zArr, vxArr, vyArr, vzArr, wxArr, wyArr, wzArr, fxArr{0}, fyArr{0},
          fzArr{0}, qxArr{0}, qyArr{0}, qzArr{0}, rArr, temperatureArr, heatFluxArr{0};
      alignas(64) std::array<autopas::OwnershipState, vecsize> ownedStateArr{};

      // broadcast of the values of particle i
      for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
        xtmp[tmpj] = xptr[indexFirst];
        ytmp[tmpj] = yptr[indexFirst];
        ztmp[tmpj] = zptr[indexFirst];

        vxtmp[tmpj] = vxptr[indexFirst];
        vytmp[tmpj] = vyptr[indexFirst];
        vztmp[tmpj] = vzptr[indexFirst];

        wxtmp[tmpj] = wxptr[indexFirst];
        wytmp[tmpj] = wyptr[indexFirst];
        wztmp[tmpj] = wzptr[indexFirst];

        rtmp[tmpj] = radiusI;

        temperaturetmp[tmpj] = temperaturePtr[indexFirst];
      }

      // loop over the verlet list from 0 to x*vecsize
      for (; joff < neighborListSize - vecsize + 1; joff += vecsize) {
        // in each iteration we calculate the interactions of particle i with
        // vecsize particles in the neighborlist of particle i starting at
        // particle joff

        // gather values of particle j
#pragma omp simd safelen(vecsize)
        for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
          xArr[tmpj] = xptr[neighborListPtr[joff + tmpj]];
          yArr[tmpj] = yptr[neighborListPtr[joff + tmpj]];
          zArr[tmpj] = zptr[neighborListPtr[joff + tmpj]];

          vxArr[tmpj] = vxptr[neighborListPtr[joff + tmpj]];
          vyArr[tmpj] = vyptr[neighborListPtr[joff + tmpj]];
          vzArr[tmpj] = vzptr[neighborListPtr[joff + tmpj]];

          wxArr[tmpj] = wxptr[neighborListPtr[joff + tmpj]];
          wyArr[tmpj] = wyptr[neighborListPtr[joff + tmpj]];
          wzArr[tmpj] = wzptr[neighborListPtr[joff + tmpj]];

          rArr[tmpj] = _PPLibrary->getRadius(typeptr[neighborListPtr[joff + tmpj]]);

          temperatureArr[tmpj] = temperaturePtr[neighborListPtr[joff + tmpj]];

          ownedStateArr[tmpj] = ownedStatePtr[neighborListPtr[joff + tmpj]];
        }

        // do omp simd with reduction of the interaction
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc, qxacc, qyacc, qzacc, heatFluxAcc, numDistanceCalculationSum,     \
                               numOverlapCalculationSum, numKernelCallsN3Sum, numKernelCallsNoN3Sum, \
                               numInnerIfTanFCallsSum, numInnerIfRollingQCallsSum, numInnerIfTorsionQCallsSum)       \
    safelen(vecsize)
        for (size_t j = 0; j < vecsize; j++) {
          // --------------------------------- main - calculation ------------------------------
          const auto ownedStateJ = ownedStateArr[j];

          const SoAFloatPrecision drx = xtmp[j] - xArr[j];
          const SoAFloatPrecision dry = ytmp[j] - yArr[j];
          const SoAFloatPrecision drz = ztmp[j] - zArr[j];

          const SoAFloatPrecision drx2 = drx * drx;
          const SoAFloatPrecision dry2 = dry * dry;
          const SoAFloatPrecision drz2 = drz * drz;

          const SoAFloatPrecision dr2 = drx2 + dry2 + drz2;
          const SoAFloatPrecision dist = std::sqrt(dr2);
          const SoAFloatPrecision overlap = (rtmp[j] + rArr[j]) - dist;

          // Mask away if distance is too large or overlap is non-positive or any particle is a dummy
          const bool cutOffMask = dist <= _cutoff and ownedStateJ != autopas::OwnershipState::dummy;
          const bool overlapIsPositive = overlap > 0;

          if constexpr (countFLOPs) {
            numDistanceCalculationSum += ownedStateJ != autopas::OwnershipState::dummy ? 1 : 0;
            numOverlapCalculationSum += (cutOffMask ? 1 : 0);
            if constexpr (newton3) {
              numKernelCallsN3Sum += (cutOffMask and overlapIsPositive ? 1 : 0);
            } else {
              numKernelCallsNoN3Sum += (cutOffMask and overlapIsPositive ? 1 : 0);
            }
          }

          if (not(cutOffMask and overlapIsPositive)) {
            continue;
          }

          const SoAFloatPrecision invDist = 1.0 / dist;
          const SoAFloatPrecision normalUnitX = drx * invDist;
          const SoAFloatPrecision normalUnitY = dry * invDist;
          const SoAFloatPrecision normalUnitZ = drz * invDist;

          const SoAFloatPrecision radiusIReduced = rtmp[j] - overlap / 2.;
          const SoAFloatPrecision radiusJReduced = rArr[j] - overlap / 2.;

          const SoAFloatPrecision relVelX = vxtmp[j] - vxArr[j];
          const SoAFloatPrecision relVelY = vytmp[j] - vyArr[j];
          const SoAFloatPrecision relVelZ = vztmp[j] - vzArr[j];

          const SoAFloatPrecision relVelDotNormalUnit =
              relVelX * normalUnitX + relVelY * normalUnitY + relVelZ * normalUnitZ;

          const SoAFloatPrecision relVelAngularIX = radiusIReduced * (normalUnitY * wztmp[j] - normalUnitZ * wytmp[j]);
          const SoAFloatPrecision relVelAngularIY = radiusIReduced * (normalUnitZ * wxtmp[j] - normalUnitX * wztmp[j]);
          const SoAFloatPrecision relVelAngularIZ = radiusIReduced * (normalUnitX * wytmp[j] - normalUnitY * wxtmp[j]);

          const SoAFloatPrecision relVelAngularJX = radiusJReduced * (normalUnitY * wzArr[j] - normalUnitZ * wyArr[j]);
          const SoAFloatPrecision relVelAngularJY = radiusJReduced * (normalUnitZ * wxArr[j] - normalUnitX * wzArr[j]);
          const SoAFloatPrecision relVelAngularJZ = radiusJReduced * (normalUnitX * wyArr[j] - normalUnitY * wxArr[j]);

          const SoAFloatPrecision tanRelVelX = relVelX + relVelAngularIX + relVelAngularJX;
          const SoAFloatPrecision tanRelVelY = relVelY + relVelAngularIY + relVelAngularJY;
          const SoAFloatPrecision tanRelVelZ = relVelZ + relVelAngularIZ + relVelAngularJZ;

          const SoAFloatPrecision tanRelVelTotalX = tanRelVelX - relVelDotNormalUnit * normalUnitX;
          const SoAFloatPrecision tanRelVelTotalY = tanRelVelY - relVelDotNormalUnit * normalUnitY;
          const SoAFloatPrecision tanRelVelTotalZ = tanRelVelZ - relVelDotNormalUnit * normalUnitZ;

          // Compute normal force as sum of contact force.
          const SoAFloatPrecision normalContactFMag =
              _elasticStiffness * overlap - _normalViscosity * relVelDotNormalUnit;

          const SoAFloatPrecision normalFX = normalContactFMag * normalUnitX;
          const SoAFloatPrecision normalFY = normalContactFMag * normalUnitY;
          const SoAFloatPrecision normalFZ = normalContactFMag * normalUnitZ;

          // Compute tangential force
          SoAFloatPrecision tanFX = -_frictionViscosity * tanRelVelTotalX;
          SoAFloatPrecision tanFY = -_frictionViscosity * tanRelVelTotalY;
          SoAFloatPrecision tanFZ = -_frictionViscosity * tanRelVelTotalZ;

          const SoAFloatPrecision tanFMag = std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ);
          const SoAFloatPrecision coulombLimit = _staticFrictionCoeff * normalContactFMag;

          if (tanFMag > coulombLimit) {
            if constexpr (countFLOPs) {
              ++numInnerIfTanFCallsSum;
            }
            const SoAFloatPrecision scale =
                _dynamicFrictionCoeff * normalContactFMag / tanFMag;
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
          if (newton3) {  // Only assignments here, as they will be subtracted later
            fxArr[j] = totalFX;
            fyArr[j] = totalFY;
            fzArr[j] = totalFZ;
          }

          // Compute torques
          // Compute frictional torque
          const SoAFloatPrecision frictionQIX = -radiusIReduced * (normalUnitY * tanFZ - normalUnitZ * tanFY);
          const SoAFloatPrecision frictionQIY = -radiusIReduced * (normalUnitZ * tanFX - normalUnitX * tanFZ);
          const SoAFloatPrecision frictionQIZ = -radiusIReduced * (normalUnitX * tanFY - normalUnitY * tanFX);

          // Compute rolling torque
          const SoAFloatPrecision radiusReduced = radiusIReduced * radiusJReduced / (radiusIReduced + radiusJReduced);
          const SoAFloatPrecision rollingRelVelX =
              -radiusReduced * (normalUnitY * (wztmp[j] - wzArr[j]) - normalUnitZ * (wytmp[j] - wyArr[j]));
          const SoAFloatPrecision rollingRelVelY =
              -radiusReduced * (normalUnitZ * (wxtmp[j] - wxArr[j]) - normalUnitX * (wztmp[j] - wzArr[j]));
          const SoAFloatPrecision rollingRelVelZ =
              -radiusReduced * (normalUnitX * (wytmp[j] - wyArr[j]) - normalUnitY * (wxtmp[j] - wxArr[j]));

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
                _rollingFrictionCoeff * (normalContactFMag) / rollingFMag;
            rollingFX *= scale;
            rollingFY *= scale;
            rollingFZ *= scale;
          }

          const SoAFloatPrecision rollingQIX = radiusReduced * (normalUnitY * rollingFZ - normalUnitZ * rollingFY);
          const SoAFloatPrecision rollingQIY = radiusReduced * (normalUnitZ * rollingFX - normalUnitX * rollingFZ);
          const SoAFloatPrecision rollingQIZ = radiusReduced * (normalUnitX * rollingFY - normalUnitY * rollingFX);

          // Compute torsion torque
          const SoAFloatPrecision torsionRelVelScalar =
              radiusReduced * ((wxtmp[j] - wxArr[j]) * normalUnitX + (wytmp[j] - wyArr[j]) * normalUnitY +
                               (wztmp[j] - wzArr[j]) * normalUnitZ);
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
                _torsionFrictionCoeff * (normalContactFMag) / torsionFMag;
            torsionFX *= scale;
            torsionFY *= scale;
            torsionFZ *= scale;
          }

          const SoAFloatPrecision torsionQIX = radiusReduced * torsionFX;
          const SoAFloatPrecision torsionQIY = radiusReduced * torsionFY;
          const SoAFloatPrecision torsionQIZ = radiusReduced * torsionFZ;

          // Apply torques
          qxacc += (frictionQIX + rollingQIX + torsionQIX);
          qyacc += (frictionQIY + rollingQIY + torsionQIY);
          qzacc += (frictionQIZ + rollingQIZ + torsionQIZ);
          if (newton3) {
            qxArr[j] = ((radiusJReduced / radiusIReduced) * frictionQIX - rollingQIX - torsionQIX);
            qyArr[j] = ((radiusJReduced / radiusIReduced) * frictionQIY - rollingQIY - torsionQIY);
            qzArr[j] = ((radiusJReduced / radiusIReduced) * frictionQIZ - rollingQIZ - torsionQIZ);
          }

          // Compute heat flux
          const SoAFloatPrecision geomMeanRadius = std::sqrt(radiusIReduced * radiusJReduced);
          const SoAFloatPrecision conductance =
              2. * _heatConductivity * std::pow((3. * _elasticStiffness * overlap * geomMeanRadius) / 4., 1. / 3.);
          const SoAFloatPrecision heatFluxI = conductance * (temperatureArr[j] - temperaturetmp[j]);

          const SoAFloatPrecision heatFluxGenerated =
              _heatGenerationFactor * std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ) *
              std::sqrt(tanRelVelTotalX * tanRelVelTotalX + tanRelVelTotalY * tanRelVelTotalY +
                        tanRelVelTotalZ * tanRelVelTotalZ);

          // Apply heat flux
          heatFluxAcc += (heatFluxI + heatFluxGenerated);
          if (newton3) {
            heatFluxArr[j] = (-heatFluxI + heatFluxGenerated);
          }
        }  // end of j loop

        // scatter the forces to where they belong, this is only needed for newton3
        if (newton3) {
#pragma omp simd safelen(vecsize)
          for (size_t tmpj = 0; tmpj < vecsize; tmpj++) {
            const size_t j = neighborListPtr[joff + tmpj];
            fxptr[j] -= fxArr[tmpj];
            fyptr[j] -= fyArr[tmpj];
            fzptr[j] -= fzArr[tmpj];

            // Caution: plus!
            qxptr[j] += qxArr[tmpj];
            qyptr[j] += qyArr[tmpj];
            qzptr[j] += qzArr[tmpj];

            heatFluxPtr[j] += heatFluxArr[tmpj];
          }
        }  // end of newton3 force / torque application
      }  // end of joff loop
    }  // end of vectorized part

    // loop over the remaining particles in the verlet list (without optimization)
    for (size_t jNeighIndex = joff; jNeighIndex < neighborListSize; ++jNeighIndex) {
      size_t j = neighborList[jNeighIndex];
      if (indexFirst == j) continue;

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
      const SoAFloatPrecision dist = std::sqrt(dr2);

      if constexpr (countFLOPs) {
        ++numDistanceCalculationSum;
      }

      if (dist > _cutoff) {
        continue;
      }

      const SoAFloatPrecision radiusJ = _PPLibrary->getRadius(typeptr[j]);
      SoAFloatPrecision overlap = (radiusI + radiusJ) - dist;
      const bool overlapIsPositive = overlap > 0;

      if constexpr (countFLOPs) {
        ++numOverlapCalculationSum;
      }

      if (not overlapIsPositive) {
        continue;
      }

      if constexpr (countFLOPs) {
        if constexpr (newton3) {
          ++numKernelCallsN3Sum;
        } else {
          ++numKernelCallsNoN3Sum;
        }
      }

      const SoAFloatPrecision invDist = 1.0 / dist;
      const SoAFloatPrecision normalUnitX = drx * invDist;
      const SoAFloatPrecision normalUnitY = dry * invDist;
      const SoAFloatPrecision normalUnitZ = drz * invDist;

      const SoAFloatPrecision radiusIReduced = radiusI - overlap / 2.;
      const SoAFloatPrecision radiusJReduced = radiusJ - overlap / 2.;

      const SoAFloatPrecision relVelX = vxptr[indexFirst] - vxptr[j];
      const SoAFloatPrecision relVelY = vyptr[indexFirst] - vyptr[j];
      const SoAFloatPrecision relVelZ = vzptr[indexFirst] - vzptr[j];

      const SoAFloatPrecision relVelDotNormalUnit =
          relVelX * normalUnitX + relVelY * normalUnitY + relVelZ * normalUnitZ;

      const SoAFloatPrecision relVelAngularIX =
          radiusIReduced * (normalUnitY * wzptr[indexFirst] - normalUnitZ * wyptr[indexFirst]);
      const SoAFloatPrecision relVelAngularIY =
          radiusIReduced * (normalUnitZ * wxptr[indexFirst] - normalUnitX * wzptr[indexFirst]);
      const SoAFloatPrecision relVelAngularIZ =
          radiusIReduced * (normalUnitX * wyptr[indexFirst] - normalUnitY * wxptr[indexFirst]);

      const SoAFloatPrecision relVelAngularJX = radiusJReduced * (normalUnitY * wzptr[j] - normalUnitZ * wyptr[j]);
      const SoAFloatPrecision relVelAngularJY = radiusJReduced * (normalUnitZ * wxptr[j] - normalUnitX * wzptr[j]);
      const SoAFloatPrecision relVelAngularJZ = radiusJReduced * (normalUnitX * wyptr[j] - normalUnitY * wxptr[j]);

      const SoAFloatPrecision tanRelVelX = relVelX + relVelAngularIX + relVelAngularJX;
      const SoAFloatPrecision tanRelVelY = relVelY + relVelAngularIY + relVelAngularJY;
      const SoAFloatPrecision tanRelVelZ = relVelZ + relVelAngularIZ + relVelAngularJZ;

      const SoAFloatPrecision tanRelVelTotalX = tanRelVelX - relVelDotNormalUnit * normalUnitX;
      const SoAFloatPrecision tanRelVelTotalY = tanRelVelY - relVelDotNormalUnit * normalUnitY;
      const SoAFloatPrecision tanRelVelTotalZ = tanRelVelZ - relVelDotNormalUnit * normalUnitZ;

      // Compute normal force as sum of contact force
      const SoAFloatPrecision normalContactFMag = _elasticStiffness * overlap - _normalViscosity * relVelDotNormalUnit;

      const SoAFloatPrecision normalFX = normalContactFMag * normalUnitX;
      const SoAFloatPrecision normalFY = normalContactFMag * normalUnitY;
      const SoAFloatPrecision normalFZ = normalContactFMag * normalUnitZ;

      // Compute tangential force
      SoAFloatPrecision tanFX = -_frictionViscosity * tanRelVelTotalX;
      SoAFloatPrecision tanFY = -_frictionViscosity * tanRelVelTotalY;
      SoAFloatPrecision tanFZ = -_frictionViscosity * tanRelVelTotalZ;

      const SoAFloatPrecision tanFMag = std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ);
      const SoAFloatPrecision coulombLimit = _staticFrictionCoeff * (normalContactFMag);

      if (tanFMag > coulombLimit) {
        if constexpr (countFLOPs) {
          ++numInnerIfTanFCallsSum;
        }
        const SoAFloatPrecision scale = _dynamicFrictionCoeff * normalContactFMag / tanFMag;
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
      if (newton3) {
        fxptr[j] -= totalFX;
        fyptr[j] -= totalFY;
        fzptr[j] -= totalFZ;
      }

      // Compute torques
      // Compute frictional torque
      const SoAFloatPrecision frictionQIX = -radiusIReduced * (normalUnitY * tanFZ - normalUnitZ * tanFY);
      const SoAFloatPrecision frictionQIY = -radiusIReduced * (normalUnitZ * tanFX - normalUnitX * tanFZ);
      const SoAFloatPrecision frictionQIZ = -radiusIReduced * (normalUnitX * tanFY - normalUnitY * tanFX);

      // Compute rolling torque
      const SoAFloatPrecision radiusReduced = radiusIReduced * radiusJReduced / (radiusIReduced + radiusJReduced);
      const SoAFloatPrecision rollingRelVelX = -radiusReduced * (normalUnitY * (wzptr[indexFirst] - wzptr[j]) -
                                                                 normalUnitZ * (wyptr[indexFirst] - wyptr[j]));
      const SoAFloatPrecision rollingRelVelY = -radiusReduced * (normalUnitZ * (wxptr[indexFirst] - wxptr[j]) -
                                                                 normalUnitX * (wzptr[indexFirst] - wzptr[j]));
      const SoAFloatPrecision rollingRelVelZ = -radiusReduced * (normalUnitX * (wyptr[indexFirst] - wyptr[j]) -
                                                                 normalUnitY * (wxptr[indexFirst] - wxptr[j]));

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
            _rollingFrictionCoeff * (normalContactFMag) / rollingFMag;
        rollingFX *= scale;
        rollingFY *= scale;
        rollingFZ *= scale;
      }

      const SoAFloatPrecision rollingQIX = radiusReduced * (normalUnitY * rollingFZ - normalUnitZ * rollingFY);
      const SoAFloatPrecision rollingQIY = radiusReduced * (normalUnitZ * rollingFX - normalUnitX * rollingFZ);
      const SoAFloatPrecision rollingQIZ = radiusReduced * (normalUnitX * rollingFY - normalUnitY * rollingFX);

      // Compute torsion torque
      const SoAFloatPrecision torsionRelVelScalar =
          radiusReduced * ((wxptr[indexFirst] - wxptr[j]) * normalUnitX + (wyptr[indexFirst] - wyptr[j]) * normalUnitY +
                           (wzptr[indexFirst] - wzptr[j]) * normalUnitZ);
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
            _torsionFrictionCoeff * normalContactFMag / torsionFMag;
        torsionFX *= scale;
        torsionFY *= scale;
        torsionFZ *= scale;
      }

      const SoAFloatPrecision torsionQIX = radiusReduced * torsionFX;
      const SoAFloatPrecision torsionQIY = radiusReduced * torsionFY;
      const SoAFloatPrecision torsionQIZ = radiusReduced * torsionFZ;

      // Apply torques
      qxacc += (frictionQIX + rollingQIX + torsionQIX);
      qyacc += (frictionQIY + rollingQIY + torsionQIY);
      qzacc += (frictionQIZ + rollingQIZ + torsionQIZ);
      if (newton3) {
        qxptr[j] += (radiusJReduced / radiusIReduced) * frictionQIX - rollingQIX - torsionQIX;
        qyptr[j] += (radiusJReduced / radiusIReduced) * frictionQIY - rollingQIY - torsionQIY;
        qzptr[j] += (radiusJReduced / radiusIReduced) * frictionQIZ - rollingQIZ - torsionQIZ;
      }

      // Compute heat flux
      const SoAFloatPrecision geomMeanRadius = std::sqrt(radiusIReduced * radiusJReduced);
      const SoAFloatPrecision conductance =
          2. * _heatConductivity * std::pow((3. * _elasticStiffness * overlap * geomMeanRadius) / 4., 1. / 3.);
      const SoAFloatPrecision heatFluxI = conductance * (temperaturePtr[j] - temperaturePtr[indexFirst]);

      const SoAFloatPrecision heatFluxGenerated =
          _heatGenerationFactor * std::sqrt(tanFX * tanFX + tanFY * tanFY + tanFZ * tanFZ) *
          std::sqrt(tanRelVelTotalX * tanRelVelTotalX + tanRelVelTotalY * tanRelVelTotalY +
                    tanRelVelTotalZ * tanRelVelTotalZ);

      // Apply heat flux
      heatFluxAcc += (heatFluxI + heatFluxGenerated);
      if (newton3) {
        heatFluxPtr[j] += (-heatFluxI + heatFluxGenerated);
      }

    }  // end of jNeighIndex loop

    if (fxacc != 0 or fyacc != 0 or fzacc != 0 or qxacc != 0 or qyacc != 0 or qzacc != 0 or heatFluxAcc != 0) {
      fxptr[indexFirst] += fxacc;
      fyptr[indexFirst] += fyacc;
      fzptr[indexFirst] += fzacc;

      qxptr[indexFirst] += qxacc;
      qyptr[indexFirst] += qyacc;
      qzptr[indexFirst] += qzacc;

      heatFluxPtr[indexFirst] += heatFluxAcc;
    }

    if constexpr (countFLOPs) {
      _aosThreadDataFLOPs[threadnum].numDistCalls += numDistanceCalculationSum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsNoN3 += numKernelCallsNoN3Sum;
      _aosThreadDataFLOPs[threadnum].numKernelCallsN3 += numKernelCallsN3Sum;
      _aosThreadDataFLOPs[threadnum].numInnerIfTanFCalls += numInnerIfTanFCallsSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfRollingQCalls += numInnerIfRollingQCallsSum;
      _aosThreadDataFLOPs[threadnum].numInnerIfTorsionQCalls += numInnerIfTorsionQCallsSum;
      _aosThreadDataFLOPs[threadnum].numOverlapCalls += numOverlapCalculationSum;
    }
  }

  /**
   * This class stores internal data for global calculations for each thread. Make sure that this data has proper size,
   * i.e. k*64 Bytes!
   *
   * No global calculation has been done in this functor, probably to be extended!
   */
  class AoSThreadDataGlobals {
   public:
    AoSThreadDataGlobals() : __remainingTo64{} {}

    void setZero() {}

   private:
    // dummy parameter to get the right size (64 bytes)
    double __remainingTo64[(64 - 4 * sizeof(double)) / sizeof(double)];
  };

  /**
   * This class stores internal data for FLOP counters for each thread. Make sure that this data has proper size, i.e.
   * k*64 Bytes!
   * The FLOP count and HitRate are not counted/calculated directly, but through helper counters (numKernelCallsNoN3,
   * numKernelCallsN3, numDistCalls, numGlobalCalcs, etc.) to reduce computational cost in the functors themselves and to
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
      numOverlapCalls = 0;
      numInnerIfRollingQCalls = 0;
      numInnerIfTorsionQCalls = 0;
      numInnerIfTanFCalls = 0;
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
     * Number of Overlap calls.
     */
    size_t numOverlapCalls = 0;

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
    size_t numInnerIfTanFCalls = 0;

   private:
    /**
     * dummy parameter to get the right size (64 bytes)
     */
    double __remainingTo64[(64 - 5 * sizeof(size_t)) / sizeof(size_t)];
  };

  const double _cutoff;
  const double _elasticStiffness;
  const double _normalViscosity;
  const double _frictionViscosity;
  const double _rollingViscosity;
  const double _torsionViscosity;
  const double _staticFrictionCoeff;
  const double _dynamicFrictionCoeff;
  const double _rollingFrictionCoeff;
  const double _torsionFrictionCoeff;
  const double _heatConductivity;
  const double _heatGenerationFactor;
  // not const because it might be reset through PPL
  double _radius = 0;

  ParticlePropertiesLibrary<SoAFloatPrecision, size_t> *_PPLibrary = nullptr;

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
};
}  // namespace demLib