/**
 * @file Thermostat.h
 * @author N. Fottner
 * @date 27/8/19
 */

#pragma once
#include <cstdlib>
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"

/**
 * Thermostat to adjust the Temperature of the Simulation
 */
namespace Thermostat {

namespace {
/**
 * Add a random velocity according to the Maxwell-Boltzmann distribution to the particle.
 *
 * @param p The particle to initialize.
 * @param averageVelocity Average velocity per dimension to be added.
 */
void maxwellBoltzmannDistribution(autopas::Particle &p, const double averageVelocity) {
  std::default_random_engine randomEngine(42);  // constant seed for repeatability
  std::normal_distribution<double> normalDistribution{0, 1};

  p.setV(autopas::utils::ArrayMath::addScalar(p.getV(), averageVelocity * normalDistribution(randomEngine)));
}
}


/**
 * Calculates temperature of system.
 * @tparam AutoPasTemplate Type of AutoPas Object (no pointer)
 * @tparam ParticlePropertiesLibraryTemplate Type of ParticlePropertiesLibrary Object (no pointer)
 * @param autopas
 * @param particlePropertiesLibrary
 * @return Temperature of system.
 */
template<class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
double calcTemperature(AutoPasTemplate &autopas, ParticlePropertiesLibraryTemplate &particlePropertiesLibrary) {
  // kinetic energy times 2
  double kineticEnergyMul2 = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : kineticEnergyMul2)
#endif
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    auto vel = iter->getV();
    kineticEnergyMul2 +=
        particlePropertiesLibrary.getMass(iter->getTypeId()) * autopas::utils::ArrayMath::dot(vel, vel);
  }
  // AutoPas works always on 3 dimensions
  constexpr unsigned int dimensions{3};
  return kineticEnergyMul2 / (autopas.getNumberOfParticles() * dimensions);
}

/**
 * Adds brownian motion to the given system.
 *
 * If useCurrentTemp is set to true the factor of the brownian motion is calculated per particle based on its mass and
 * the system's temperature. Otherwise a constant factor of 0.1 is used.
 * Set this to false if the system is initialized without velocities.
 *
 * @tparam AutoPasTemplate Type of AutoPas Object (no pointer)
 * @tparam ParticlePropertiesLibraryTemplate Type of ParticlePropertiesLibrary Object (no pointer)
 * @param autopas
 * @param particlePropertiesLibrary
 * @param useCurrentTemp
 */
template<class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void addBrownianMotion(AutoPasTemplate &autopas,
                       ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
                       bool useCurrentTemp) {
  // factors for the brownian motion per particle type.
  std::map<size_t, double> factors;
  if (useCurrentTemp) {
    // brownian motion with disturbance depending on current temperature and mass
    constexpr double boltzmannConst = 1.38064852e-23;  // Boltzmann constant in J/K
    double currentTempMulKB = calcTemperature(autopas, particlePropertiesLibrary) * boltzmannConst;
    for (auto typeID : particlePropertiesLibrary.getTypes()) {
      factors.emplace(typeID, std::sqrt(currentTempMulKB / particlePropertiesLibrary.getMass(typeID)));
    }
  } else {
    // simple version of brownian motion with constant disturbance for all particles
    for (auto typeID : particlePropertiesLibrary.getTypes()) {
      factors.emplace(typeID, 0.1);
    }
  }
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    maxwellBoltzmannDistribution(*iter, factors[iter->getTypeId()]);
  }
}

/**
 * Scales velocity of particles towards a gived temperature.
 * @tparam AutoPasTemplate Type of AutoPas Object (no pointer)
 * @tparam ParticlePropertiesLibraryTemplate Type of ParticlePropertiesLibrary Object (no pointer)
 * @param autopas
 * @param particlePropertiesLibrary
 * @param targetTemperature
 * @param deltaTemperature Maximum temperature change.
 */
template<class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void apply(AutoPasTemplate &autopas,
           ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
           double targetTemperature,
           double deltaTemperature) {
  double currentTemperature = calcTemperature(autopas, particlePropertiesLibrary);
  double nextTargetTemperature = currentTemperature + deltaTemperature;
  // check if we are already in the vicinity of our target or if we still need full steps
  if (std::abs(nextTargetTemperature) > std::abs(targetTemperature)) {
    nextTargetTemperature = targetTemperature;
  }
  double scaling = std::sqrt(nextTargetTemperature / currentTemperature);
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    iter->setV(autopas::utils::ArrayMath::mulScalar(iter->getV(), scaling));
  }
}
};
