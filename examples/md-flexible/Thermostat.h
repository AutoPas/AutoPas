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
 * WIP
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
class Thermostat {
 public:
  using ThermostatFloatType = typename ParticlePropertiesLibraryTemplate::ParticlePropertiesLibraryFloatType;
  using ThermostatIntType = typename ParticlePropertiesLibraryTemplate::ParticlePropertiesLibraryIntType;

  /**
   * Constructor if target Temperature is specified
   * @param tInit initialTemperature of System
   * @param initBm initialization with BM or other Formula
   * @param tTarget target Temperature of System
   * @param deltaTemp
   * @param particlePropertiesLibrary
   * */
  //@todo FRAGE : add Parsing option to switch beetween _initParticlesOnlyWithBrownianMotion on & off when initThermo is spezified as
  // true
  Thermostat(ThermostatFloatType tInit, bool initBm, ThermostatFloatType tTarget, ThermostatFloatType deltaTemp,
             const ParticlePropertiesLibraryTemplate &particlePropertiesLibrary);

  /**
   * Default Destructor
   * */
  ~Thermostat() = default;

  /**
   * Copy Constructor
   * @param ThermostatCopy
   * */
  Thermostat(const Thermostat &ThermostatCopy) = default;

  /**
   * Copy Assignment Constructor
   * @param thermo
   * @return
   * */
  Thermostat &operator=(const Thermostat &thermo) = default;

  /**
   * Initializes the Simulation according to brownian movement.
   * @param ps Particle system to be initialized.
   */
  void initialize(AutoPasTemplate &autopas);

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

  /**
   * Add a random velocity according to the Maxwell-Boltzmann distribution to the
   * particles, with a given mean velocity.
   * code taken from the MolSim Course
   *
   * @param p The particle to initialize.
   * @param factor
   */
  void maxwellBoltzmannDistribution(autopas::Particle &p, const ThermostatFloatType factor) {
    p.setV(autopas::ArrayMath::addScalar(p.getV(), factor * gaussDeviate()));
  }

  /**
   * helper function for MaxwellBoltzmannDistribution().
   * Generates a gauss deviate, i.e. values according to the normal distribution.
   * code taken from:
   * Griebel et. al.: Numerical Simulation in Molecular Dynamics, p. 427
   */
  static ThermostatFloatType gaussDeviate() {
    ThermostatFloatType a1, a2, s, r, b1;
    static bool iset = false;
    static ThermostatFloatType b2;

    if (!iset) {
      do {
        a1 = 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
        a2 = 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
        r = a1 * a1 + a2 * a2;
      } while (r >= 1.0);
      s = sqrt(-2.0 * log(r) / r);
      b1 = a1 * s;
      b2 = a2 * s;
      iset = true;
      return b1;
    } else {
      iset = false;
      return b2;
    }
  }

 private:
  /**
   * Initial temperature + deltaTemp * simstep / n_thermostat.
   */
  ThermostatFloatType _tInit;

  /**
   * Target temperature.
   */
  const ThermostatFloatType _tTarget;

  /**
   * Specifies, how the velocity values will be initialized: With Brownian Motion or according to the temperature
   */
  const bool _initParticlesOnlyWithBrownianMotion;

  /**
   * Temperature difference per thermostat application until t_target is reached.
   */
  const ThermostatFloatType _deltaTemp;

  /**
   * ParticlePropertiesLibrary to access Mass of Particles
   * */
  ParticlePropertiesLibraryTemplate _particlePropertiesLibrary;
};
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void Thermostat<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::initialize(AutoPasTemplate &autopas) {
  constexpr ThermostatFloatType boltzmannConst = 1.38064852 * std::pow(10, -23);  // Boltzmann constant in J/K
  if (_initParticlesOnlyWithBrownianMotion) {
// we consider: "However, the initialization with the Brownian Motion should be optional." means factor=0.1
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
    for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
      maxwellBoltzmannDistribution(*iter, 0.1);
    }
  } else {
    ThermostatFloatType currentTempMulKB = calcTemperature(autopas) * boltzmannConst;
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
    for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
      ThermostatFloatType factor =
          std::sqrt(currentTempMulKB / _particlePropertiesLibrary.getMass(iter->getTypeId()));
      maxwellBoltzmannDistribution(*iter, factor);
    }
  }
}

template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void Thermostat<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::apply(AutoPasTemplate &autopas) {
  ThermostatFloatType temp = calcTemperature(autopas);
  static ThermostatFloatType scaling;
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
  return kineticEnergyMul2 / (autopas.getNumberOfParticles() * 3);  // AutoPas works always on 3 dimensions
}

template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
Thermostat<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::Thermostat(
    ThermostatFloatType tInit, bool initBm, ThermostatFloatType tTarget, ThermostatFloatType deltaTemp,
    const ParticlePropertiesLibraryTemplate &particlePropertiesLibrary)
    : _tInit(tInit),
      _tTarget(tTarget),
      _initParticlesOnlyWithBrownianMotion(initBm),
      _deltaTemp(deltaTemp),
      _particlePropertiesLibrary(particlePropertiesLibrary) {}
