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
#include "autopas/utils/WrapOpenMP.h"

/**
 * Thermostat to adjust the Temperature of the Simulation.
 */
namespace Thermostat {

namespace {
/**
 * Add a random velocity according to the Maxwell-Boltzmann distribution to the particle.
 *
 * @param p The particle to initialize.
 * @param averageVelocity Average velocity per dimension to be added.
 * @param randomEngine Random engine used for the generation of the velocity.
 * @param normalDistribution Distribution used for constructing the maxwell boltzmann distribution.
 */
 KOKKOS_INLINE_FUNCTION
void addMaxwellBoltzmannDistributedVelocity(ParticleType &p, const double averageVelocity,
                                            std::default_random_engine &randomEngine,
                                            std::normal_distribution<double> &normalDistribution) {
  // when adding independent normally distributed values to all velocity components
  // the velocity change is maxwell boltzmann distributed
  std::array<double, 3> randomVelocity{};
  for (double &v : randomVelocity) {
    v = averageVelocity * normalDistribution(randomEngine);
  }
  p.setV(autopas::utils::ArrayMath::add(p.getV(), randomVelocity));
}
}  // namespace

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
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(+ : kineticEnergyMul2) default(none) shared(autopas, particlePropertiesLibrary)
#endif
  autopas.reduceParallel([&] (ParticleType &p, double reductionValue) {
    auto vel = p.getV();
    reductionValue +=
        particlePropertiesLibrary.getMass(p.getTypeId()) * autopas::utils::ArrayMath::dot(vel, vel);
  }, kineticEnergyMul2, autopas::IteratorBehavior::ownedOrHalo);

  // AutoPas works always on 3 dimensions
  constexpr unsigned int dimensions{3};
  return kineticEnergyMul2 / (autopas.getNumberOfParticles() * dimensions);
}

/**
 * Calculates temperature of system, for each component separately.
 * Assuming dimension-less units and Boltzmann constant = 1.
 * @tparam AutoPasTemplate Type of AutoPas Object (no pointer)
 * @tparam ParticlePropertiesLibraryTemplate Type of ParticlePropertiesLibrary Object (no pointer)
 * @param autopas
 * @param particlePropertiesLibrary
 * @return map of: particle typeID -> temperature for this type
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
auto calcTemperatureComponent(const AutoPasTemplate &autopas,
                              ParticlePropertiesLibraryTemplate &particlePropertiesLibrary) {
  // map of: particle typeID -> kinetic energy times 2 for this type
  std::map<size_t, double> kineticEnergyMul2Map;
  // map of: particle typeID -> number of particles of this type
  std::map<size_t, size_t> numParticleMap;

  //TODO lgaertner sucks ass... improve with allowing multiple reductions per loop?
  for (const auto &typeID : particlePropertiesLibrary.getTypes()) {
    double kemRedVal = 0.;
    size_t npmRedVal = 0ul;
    autopas.reduce([&] (ParticleType &p, double redVal) {
        if (p.getTypeId() == typeID) {
            auto vel = p.getV();
            redVal +=
                particlePropertiesLibrary.getMass(typeID) * autopas::utils::ArrayMath::dot(vel, vel);
          }
        }, kemRedVal, autopas::IteratorBehavior::ownedOrHalo);

    autopas.reduce([&] (ParticleType &p, size_t redVal) {
      if (p.getTypeId() == typeID) {
        redVal += 1;
      }
    }, npmRedVal, autopas::IteratorBehavior::ownedOrHalo);

    kineticEnergyMul2Map[typeID] = 0.;
    numParticleMap[typeID] = 0ul;
  }
  // AutoPas works always on 3 dimensions
  constexpr unsigned int dimensions{3};

  auto kineticEnergyAndParticleMaps = std::make_tuple(kineticEnergyMul2Map.begin(), numParticleMap.begin());

  for (auto [kinEIter, numParIter] = kineticEnergyAndParticleMaps; kinEIter != kineticEnergyMul2Map.end();
       ++kinEIter, ++numParIter) {
    kinEIter->second /= numParIter->second * dimensions;
  }
  return kineticEnergyMul2Map;
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
 * @param targetTemperature temperature of the system after applying the function on a system with temperature = 0.
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void addBrownianMotion(AutoPasTemplate &autopas, ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
                       const double targetTemperature) {
  // factors for the brownian motion per particle type.
  std::map<size_t, double> factors;
  // brownian motion with disturbance depending on current temperature and mass
  for (auto typeID : particlePropertiesLibrary.getTypes()) {
    factors.emplace(typeID, std::sqrt(targetTemperature / particlePropertiesLibrary.getMass(typeID)));
  }
//#ifdef AUTOPAS_OPENMP
//#pragma omp parallel default(none) shared(autopas, factors)
//#endif
  {
    // we use a constant seed for repeatability.
    // we need one random engine and distribution per thread
    std::default_random_engine randomEngine(42 + autopas::autopas_get_thread_num()); // TODO lgaertner ?? why what how?
    std::normal_distribution<double> normalDistribution{0, 1};
    autopas.forEachParallel([&] (ParticleType &p) {
//      TODO lgaertner: fix random.h access for gpu usage
      addMaxwellBoltzmannDistributedVelocity(p, factors[p.getTypeId()], randomEngine, normalDistribution);
    }, autopas::IteratorBehavior::ownedOrHalo);
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
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void apply(AutoPasTemplate &autopas, ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
           const double targetTemperature, const double deltaTemperature) {
  auto currentTemperatureMap = calcTemperatureComponent(autopas, particlePropertiesLibrary);
  // make sure we work with a positive delta
  const double deltaTemperaturePositive = std::abs(deltaTemperature);
  decltype(currentTemperatureMap) scalingMap;

  for (const auto &[particleTypeID, currentTemperature] : currentTemperatureMap) {
    double nextTargetTemperature;
    // check if we are already in the vicinity of our target or if we still need full steps
    if (std::abs(currentTemperature - targetTemperature) < std::abs(deltaTemperature)) {
      nextTargetTemperature = targetTemperature;
    } else {
      // make sure we scale in the right direction
      nextTargetTemperature = currentTemperature < targetTemperature ? currentTemperature + deltaTemperaturePositive
                                                                     : currentTemperature - deltaTemperaturePositive;
    }
    scalingMap[particleTypeID] = std::sqrt(nextTargetTemperature / currentTemperature);
  }
  autopas.forEachParallel([&] (ParticleType &p) {
    p.setV(autopas::utils::ArrayMath::mulScalar(p.getV(), scalingMap[p.getTypeId()]));
  });
}
}  // namespace Thermostat
