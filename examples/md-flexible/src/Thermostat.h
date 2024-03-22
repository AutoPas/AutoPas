/**
 * @file Thermostat.h
 * @author N. Fottner
 * @date 27/8/19
 */

#pragma once
#include <cstdlib>

#include "TypeDefinitions.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapMPI.h"
#include "autopas/utils/WrapOpenMP.h"

/**
 * Thermostat to adjust the Temperature of the Simulation.
 */
namespace Thermostat {
/**
 * Calculates temperature of system.
 * Assuming dimension-less units and Boltzmann constant = 1.
 * @tparam AutoPasTemplate Type of AutoPas Object (no pointer)
 * @tparam ParticlePropertiesLibraryTemplate Type of ParticlePropertiesLibrary Object (no pointer)
 * @param autopas
 * @param particlePropertiesLibrary
 * @return Temperature of system.
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
double calcTemperature(const AutoPasTemplate &autopas, ParticlePropertiesLibraryTemplate &particlePropertiesLibrary) {
  // kinetic energy times 2
  double kineticEnergyMul2 = 0;
  AUTOPAS_OPENMP(parallel reduction(+ : kineticEnergyMul2) default(none) shared(autopas, particlePropertiesLibrary))
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    const auto vel = iter->getV();
#if MD_FLEXIBLE_MODE == MULTISITE
    const auto angVel = iter->getAngularVel();
#endif
    kineticEnergyMul2 +=
        particlePropertiesLibrary.getMolMass(iter->getTypeId()) * autopas::utils::ArrayMath::dot(vel, vel);
#if MD_FLEXIBLE_MODE == MULTISITE
    kineticEnergyMul2 += autopas::utils::ArrayMath::dot(particlePropertiesLibrary.getMomentOfInertia(iter->getTypeId()),
                                                        autopas::utils::ArrayMath::mul(angVel, angVel));
#endif
  }
  // md-flexible's molecules have 3 DoF for translational velocity and optionally 3 additional rotational DoF
  constexpr unsigned int degreesOfFreedom {
#if MD_FLEXIBLE_MODE == MULTISITE
    6
#else
    3
#endif
  };
  return kineticEnergyMul2 / (autopas.getNumberOfParticles() * degreesOfFreedom);
}

/**
 * Calculates temperature of system, for each component separately.
 *
 * Kinetic Energy for each molecule is
 *    1/2 * mass * dot(vel, vel) + 1/2 Sum_{0 <= i < 3} MoI_i * angVel_i^2
 * where MoI is the diagonal Moment of Inertia. This formula comes from Rapport, The Art of MD, equation (8.2.34).
 *
 * The second term is only applied for Multi-Site MD.
 *
 * Assuming dimension-less units and Boltzmann constant = 1.
 *
 * @warning This function is only usable in multi-site mode, where mixing parameters are stored as per site types
 *
 * @tparam AutoPasTemplate Type of AutoPas Object (no pointer)
 * @tparam ParticlePropertiesLibraryTemplate Type of ParticlePropertiesLibrary Object (no pointer)
 * @param autopas
 * @param particlePropertiesLibrary ignored if using single-site mode
 * @return map of: particle typeID -> temperature for this type
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
auto calcTemperatureComponent(const AutoPasTemplate &autopas,
                              ParticlePropertiesLibraryTemplate &particlePropertiesLibrary) {
#if MD_FLEXIBLE_MODE == SINGLESITE
  autopas::utils::ExceptionHandler::exception("calcTemperatureComponent is no longer supported for single-site molecules!");
#endif

  using autopas::utils::ArrayMath::dot;
  using namespace autopas::utils::ArrayMath::literals;

  // map of: particle typeID -> kinetic energy times 2 for this type
  std::map<size_t, double> kineticEnergyMul2Map;
  // map of: particle typeID -> number of particles of this type
  std::map<size_t, size_t> numParticleMap;

  for (int typeID = 0; typeID < particlePropertiesLibrary.getNumberRegisteredSiteTypes(); typeID++) {
    kineticEnergyMul2Map[typeID] = 0.;
    numParticleMap[typeID] = 0ul;
  }

  AUTOPAS_OPENMP(parallel) {
    // create aggregators for each thread
    std::map<size_t, double> kineticEnergyMul2MapThread;
    std::map<size_t, size_t> numParticleMapThread;
    for (int typeID = 0; typeID < particlePropertiesLibrary.getNumberRegisteredSiteTypes(); typeID++) {
      kineticEnergyMul2MapThread[typeID] = 0.;
      numParticleMapThread[typeID] = 0ul;
    }
    // parallel iterators
    for (auto iter = autopas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      const auto &vel = iter->getV();
      kineticEnergyMul2MapThread.at(iter->getTypeId()) +=
          particlePropertiesLibrary.getMolMass(iter->getTypeId()) * dot(vel, vel);
#if MD_FLEXIBLE_MODE == MULTISITE
      // add contribution from angular momentum
      const auto &angVel = iter->getAngularVel();
      kineticEnergyMul2MapThread.at(iter->getTypeId()) +=
          dot(particlePropertiesLibrary.getMomentOfInertia(iter->getTypeId()), angVel * angVel);
#endif
      numParticleMapThread.at(iter->getTypeId())++;
    }
    // manual reduction
    AUTOPAS_OPENMP(critical) {
      for (int typeID = 0; typeID < particlePropertiesLibrary.getNumberRegisteredSiteTypes(); typeID++) {
        kineticEnergyMul2Map[typeID] += kineticEnergyMul2MapThread[typeID];
        numParticleMap[typeID] += numParticleMapThread[typeID];
      }
    }
  }
  // md-flexible's molecules have 3 DoF for translational velocity and optionally 3 additional rotational DoF
#if MD_FLEXIBLE_MODE == MULTISITE
  constexpr unsigned int degreesOfFreedom{6};
#else
  constexpr unsigned int degreesOfFreedom{3};
#endif

  for (int typeID = 0; typeID < particlePropertiesLibrary.getNumberRegisteredSiteTypes(); typeID++) {
    // workaround for MPICH: send and receive buffer must not be the same.
    autopas::AutoPas_MPI_Allreduce(AUTOPAS_MPI_IN_PLACE, &kineticEnergyMul2Map[typeID], 1, AUTOPAS_MPI_DOUBLE,
                                   AUTOPAS_MPI_SUM, AUTOPAS_MPI_COMM_WORLD);

    autopas::AutoPas_MPI_Allreduce(AUTOPAS_MPI_IN_PLACE, &numParticleMap[typeID], 1, AUTOPAS_MPI_UNSIGNED_LONG,
                                   AUTOPAS_MPI_SUM, AUTOPAS_MPI_COMM_WORLD);
  }

  auto kineticEnergyAndParticleMaps = std::make_tuple(kineticEnergyMul2Map.begin(), numParticleMap.begin());

  for (auto [kineticEnergyMapIter, numParticleMapIter] = kineticEnergyAndParticleMaps;
       kineticEnergyMapIter != kineticEnergyMul2Map.end(); ++kineticEnergyMapIter, ++numParticleMapIter) {
    // The calculation below assumes that the Boltzmann constant is 1.
    kineticEnergyMapIter->second /= static_cast<double>(numParticleMapIter->second) * degreesOfFreedom;
  }
  return kineticEnergyMul2Map;
}

/**
 * Adds brownian motion to the given system.
 *
 * This is achieved by each degree-of-freedom for the Kinetic Energy being sampled via the normal distribution and
 * scaled appropriately, as determined by the equipartition theorem.
 *
 * For multi-site MD we assume that the kinetic energy of the system can be split equally into translational and
 * rotational kinetic energies.
 *
 * For translational velocity, each degree-of-freedom is sampled via the normal distribution and scaled equally,
 * dependant on the temperature and mass.
 *
 * For angular velocity, each degree-of-freedom is sampled via the normal distribution and scaled dependent on the
 * temperature and the corresponding component of the diagonalized Moment-Of-Inertia.
 *
 * In all cases, we assume a Boltzmann constant of 1.
 *
 * @tparam AutoPasTemplate Type of AutoPas Object (no pointer)
 * @tparam ParticlePropertiesLibraryTemplate Type of ParticlePropertiesLibrary Object (no pointer)
 * @param autopas
 * @param particlePropertiesLibrary ignored if using single-site mode
 * @param targetTemperature temperature of the system after applying the function on a system with temperature = 0.
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void addBrownianMotion(AutoPasTemplate &autopas, ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
                       const double targetTemperature) {
  using namespace autopas::utils::ArrayMath::literals;

#if MD_FLEXIBLE_MODE == SINGLESITE
  AUTOPAS_OPENMP(parallel default(none) shared(autopas, targetTemperature)) {
    // we use a constant seed for repeatability.
    // we need one random engine and distribution per thread
    std::default_random_engine randomEngine(42 + autopas::autopas_get_thread_num());
    std::normal_distribution<double> normalDistribution{0, 1};

    // This needlessly requires a sqrt per molecule -> If we move this to object creation we can do this once per object,
    for (auto mol = autopas.begin(autopas::IteratorBehavior::owned); mol.isValid(); ++mol) {
      const std::array<double, 3> normal3DVecTranslational = {
          normalDistribution(randomEngine), normalDistribution(randomEngine), normalDistribution(randomEngine)};
      const auto scale = std::sqrt(targetTemperature / mol->getMass());
      mol->addV(normal3DVecTranslational * scale);
    }
  }
#elif MD_FLEXIBLE_MODE == MULTISITE // If Multisite, use old handling with PPL for this
  // Generate map(s) of molecule type Id to scaling factors
  std::map<size_t, double> translationalVelocityScale;
  std::map<size_t, std::array<double, 3>> rotationalVelocityScale;

  for (int typeID = 0; typeID < particlePropertiesLibrary.getNumberRegisteredSiteTypes(); typeID++) {
    translationalVelocityScale.emplace(typeID,
                                       std::sqrt(targetTemperature / particlePropertiesLibrary.getMolMass(typeID)));

    const auto momentOfInertia = particlePropertiesLibrary.getMomentOfInertia(typeID);
    const std::array<double, 3> scale{std::sqrt(targetTemperature / momentOfInertia[0]),
                                      std::sqrt(targetTemperature / momentOfInertia[1]),
                                      std::sqrt(targetTemperature / momentOfInertia[2])};
    rotationalVelocityScale.emplace(typeID, scale);
  }

  AUTOPAS_OPENMP(parallel default(none) shared(autopas, translationalVelocityScale, rotationalVelocityScale)) {
    // we use a constant seed for repeatability.
    // we need one random engine and distribution per thread
    std::default_random_engine randomEngine(42 + autopas::autopas_get_thread_num());
    std::normal_distribution<double> normalDistribution{0, 1};
    for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
      const std::array<double, 3> normal3DVecTranslational = {
          normalDistribution(randomEngine), normalDistribution(randomEngine), normalDistribution(randomEngine)};
      iter->addV(normal3DVecTranslational * translationalVelocityScale[iter->getTypeId()]);
      const std::array<double, 3> normal3DVecRotational = {
          normalDistribution(randomEngine), normalDistribution(randomEngine), normalDistribution(randomEngine)};
      iter->addAngularVel(normal3DVecRotational * rotationalVelocityScale[iter->getTypeId()]);
    }
  }

#endif


}

/**
 * Scales velocity of particles towards a given temperature. For Multi-site simulations, angular velocity is also
 * scaled.
 * @tparam AutoPasTemplate Type of AutoPas Object (no pointer)
 * @tparam ParticlePropertiesLibraryTemplate Type of ParticlePropertiesLibrary Object (no pointer)
 * @param autopas
 * @param particlePropertiesLibrary ignored if using single-site mode
 * @param targetTemperature
 * @param deltaTemperature Maximum temperature change.
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void apply(AutoPasTemplate &autopas, ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
           const double targetTemperature, const double deltaTemperature) {
  using namespace autopas::utils::ArrayMath::literals;
  using autopas::utils::ArrayMath::dot;

#if MD_FLEXIBLE_MODE == SINGLESITE
  // get total number of particles
  double numParticles = autopas.getNumberOfParticles(autopas::IteratorBehavior::owned);
  autopas::AutoPas_MPI_Allreduce(AUTOPAS_MPI_IN_PLACE, &numParticles, 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_SUM, AUTOPAS_MPI_COMM_WORLD);

  // calculate current temperature:
  double kineticEnergy = 0.;
  AUTOPAS_OPENMP(parallel default(none) reduction( + : kineticEnergy) shared(autopas)) {
    for (auto mol = autopas.begin(autopas::IteratorBehavior::owned); mol.isValid(); ++mol) {
      const auto &vel = mol->getV();
      const auto &mass = mol->getMass();
      kineticEnergy += mass * dot(vel, vel);
    }
  }
  autopas::AutoPas_MPI_Allreduce(AUTOPAS_MPI_IN_PLACE, &kineticEnergy, 1, AUTOPAS_MPI_DOUBLE, AUTOPAS_MPI_SUM, AUTOPAS_MPI_COMM_WORLD);

  // Assume Boltzmann constant is 1. Assume 3 DoF (needs to be modified to 6 for multi-site)
  const auto currentTemperature = kineticEnergy / (numParticles * 3);

  // make sure we work with a positive delta
  const double absoluteDeltaTemperature = std::abs(deltaTemperature);

  // If the current temperature is within absoluteDeltaTemperature of target temperature,
  //      set immediate target temperature to target temperature.
  //
  // Else, set immediate target temperature to absoluteDeltaTemperature towards the target temperature from the
  // current temperature.
  const auto immediateTargetTemperature =
      currentTemperature < targetTemperature
          ? std::min(currentTemperature + absoluteDeltaTemperature, targetTemperature)
          : std::max(currentTemperature - absoluteDeltaTemperature, targetTemperature);

  const auto scalingFactor = std::sqrt(immediateTargetTemperature / currentTemperature);

  // Apply scaling
  AUTOPAS_OPENMP(parallel default(none) shared(autopas, scalingFactor))
  for (auto mol = autopas.begin(autopas::IteratorBehavior::owned); mol.isValid(); ++mol) {
    mol->setV(mol->getV() * scalingFactor);
  }

#elif MD_FLEXIBLE_MODE == MULTISITE
  const auto currentTemperatureMap = calcTemperatureComponent(autopas, particlePropertiesLibrary);

  // make sure we work with a positive delta
  const double absoluteDeltaTemperature = std::abs(deltaTemperature);

  std::remove_const_t<decltype(currentTemperatureMap)> scalingMap;

  for (const auto &[particleTypeID, currentTemperature] : currentTemperatureMap) {
    // If the current temperature is within absoluteDeltaTemperature of target temperature,
    //      set immediate target temperature to target temperature.
    //
    // Else, set immediate target temperature to absoluteDeltaTemperature towards the target temperature from the
    // current temperature.

    const auto immediateTargetTemperature =
        currentTemperature < targetTemperature
            ? std::min(currentTemperature + absoluteDeltaTemperature, targetTemperature)
            : std::max(currentTemperature - absoluteDeltaTemperature, targetTemperature);
    // Determine a scaling factor for each particle type.
    scalingMap[particleTypeID] = std::sqrt(immediateTargetTemperature / currentTemperature);
  }

  // Scale velocities (and angular velocities) with the scaling map
  AUTOPAS_OPENMP(parallel default(none) shared(autopas, scalingMap))
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    iter->setV(iter->getV() * scalingMap[iter->getTypeId()]);
#if MD_FLEXIBLE_MODE == MULTISITE
    iter->setAngularVel(iter->getAngularVel() * scalingMap[iter->getTypeId()]);
#endif
  }
#endif
}
}  // namespace Thermostat
