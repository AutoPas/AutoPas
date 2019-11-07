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
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
class Thermostat {
 public:
  using ThermostatFloatType = typename ParticlePropertiesLibraryTemplate::ParticlePropertiesLibraryFloatType;
  using ThermostatIntType = typename ParticlePropertiesLibraryTemplate::ParticlePropertiesLibraryIntType;

  /**
   * Constructor if target Temperature is specified
   * @param tInit initialTemperature of System
   * @param tTarget target Temperature of System
   * @param deltaTemp
   * @param particlePropertiesLibrary
   */
  Thermostat(ThermostatFloatType tInit, ThermostatFloatType tTarget, ThermostatFloatType deltaTemp,
             const ParticlePropertiesLibraryTemplate &particlePropertiesLibrary);

  /**
   * Default Destructor
   */
  ~Thermostat() = default;

  /**
   * Copy Constructor
   * @param ThermostatCopy
   */
  Thermostat(const Thermostat &ThermostatCopy) = default;

  /**
   * Copy Assignment Constructor
   * @param thermo
   * @return
   */
  Thermostat &operator=(const Thermostat &thermo) = default;

  /**
   * Adds brownian motion.
   * @param ps Particle system to be initialized.
   */

  /**
   * Adds brownian motion to the given system.
   *
   * If useCurrentTemp is set to true the factor of the brownian motion is calculated per particle based on its mass and
   * the system's temperature. Otherwise a constant factor of 0.1 is used.
   * Set this to false if the system is initialized without velocities.
   *
   * @param autopas
   * @param useCurrentTemp
   */
  void addBrownianMotion(AutoPasTemplate &autopas, bool useCurrentTemp);

  /**
   * Scales velocity of particles to reach desired temperature.
   * @param ps Particle system
   */
  void apply(AutoPasTemplate &autopas);

  /**
   * Calculates temperature of system.
   * @param c Particle System.
   * @return Temperature of system.
   */
  ThermostatFloatType calcTemperature(AutoPasTemplate &autopas);

 private:
  /**
   * Add a random velocity according to the Maxwell-Boltzmann distribution to the
   * particles, with a given mean velocity.
   *
   * @param p The particle to initialize.
   * @param factor
   */
  void maxwellBoltzmannDistribution(autopas::Particle &p, const ThermostatFloatType factor) {
    std::default_random_engine randomEngine(42);  // constant seed for repeatability
    std::normal_distribution<double> normalDistribution{0, 1};

    p.setV(autopas::ArrayMath::addScalar(p.getV(), factor * normalDistribution(randomEngine)));
  }

  /**
   * Initial temperature + deltaTemp * simstep / n_thermostat.
   */
  ThermostatFloatType _tInit;

  /**
   * Target temperature.
   */
  const ThermostatFloatType _tTarget;

  /**
   * Temperature difference per thermostat application until t_target is reached.
   */
  const ThermostatFloatType _deltaTemp;

  /**
   * ParticlePropertiesLibrary to access Mass of Particles
   */
  ParticlePropertiesLibraryTemplate _particlePropertiesLibrary;
};
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void Thermostat<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::addBrownianMotion(AutoPasTemplate &autopas,
                                                                                       bool useCurrentTemp) {
  // factors for the brownian motion per particle type.
  std::map<size_t, ThermostatFloatType> factors;
  if (useCurrentTemp) {
    // brownian motion with disturbance depending on current temperature and mass
    constexpr ThermostatFloatType boltzmannConst = 1.38064852e-23;  // Boltzmann constant in J/K
    ThermostatFloatType currentTempMulKB = calcTemperature(autopas) * boltzmannConst;
    for (auto typeID : _particlePropertiesLibrary.getTypes()) {
      factors.emplace(typeID, std::sqrt(currentTempMulKB / _particlePropertiesLibrary.getMass(typeID)));
    }
  } else {
    // simple version of brownian motion with constant disturbance for all particles
    for (auto typeID : _particlePropertiesLibrary.getTypes()) {
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

template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void Thermostat<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::apply(AutoPasTemplate &autopas) {
  ThermostatFloatType temp = calcTemperature(autopas);
  ThermostatFloatType scaling;
  if (_tInit == _tTarget) {
    scaling = std::sqrt(_tTarget / temp);
  } else {
    if (std::abs(_tInit + _deltaTemp) > std::abs(_tTarget)) {
      _tInit = _tTarget;
    } else {
      _tInit = _tInit + _deltaTemp;
    }
    scaling = std::sqrt(_tInit / temp);
  }
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    iter->setV(autopas::ArrayMath::mulScalar(iter->getV(), scaling));
  }
}

template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
typename ParticlePropertiesLibraryTemplate::ParticlePropertiesLibraryFloatType
Thermostat<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::calcTemperature(AutoPasTemplate &autopas) {
  // kinetic energy times 2
  ThermostatFloatType kineticEnergyMul2 = 0;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : kineticEnergyMul2)
#endif
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    auto vel = iter->getV();
    kineticEnergyMul2 += _particlePropertiesLibrary.getMass(iter->getTypeId()) * autopas::ArrayMath::dot(vel, vel);
  }
  // AutoPas works always on 3 dimensions
  constexpr unsigned int dimensions{3};
  return kineticEnergyMul2 / (autopas.getNumberOfParticles() * dimensions);
}

template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
Thermostat<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::Thermostat(
    ThermostatFloatType tInit, ThermostatFloatType tTarget, ThermostatFloatType deltaTemp,
    const ParticlePropertiesLibraryTemplate &particlePropertiesLibrary)
    : _tInit(tInit), _tTarget(tTarget), _deltaTemp(deltaTemp), _particlePropertiesLibrary(particlePropertiesLibrary) {}
