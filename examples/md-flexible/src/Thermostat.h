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

#ifdef AUTOPAS_ENABLE_KOKKOS
#include <Kokkos_DualView.hpp>
#include <Kokkos_Random.hpp>
#endif

/**
 * Thermostat to adjust the Temperature of the Simulation.
 */
namespace Thermostat {

#ifdef AUTOPAS_ENABLE_KOKKOS
// TODO: it might make sense to outsource this to a common location to avoid redefining this over and over again
#ifdef KOKKOS_ENABLE_CUDA
using DeviceSpace = Kokkos::CudaSpace;
constexpr bool ForEachHostFlag = false;
#else
using DeviceSpace = Kokkos::HostSpace;
constexpr bool ForEachHostFlag = true;
#endif

template <class ParticleType>
struct CalcTemperatureFunctor {
  KOKKOS_INLINE_FUNCTION
  void operator()(int i, const autopas::utilsKokkos::KokkosStorage<ParticleType>& storage, double& localKinetic) const {
    const auto velX = storage.template operator()<ParticleType::AttributeNames::velocityX, ForEachHostFlag>(i);
    const auto velY = storage.template operator()<ParticleType::AttributeNames::velocityY, ForEachHostFlag>(i);
    const auto velZ = storage.template operator()<ParticleType::AttributeNames::velocityZ, ForEachHostFlag>(i);

    const auto mass = storage.template operator()<ParticleType::AttributeNames::mass, ForEachHostFlag>(i);

    localKinetic += mass * (velX * velX + velY * velY + velZ * velZ);
  }
};

template <class ParticleType>
struct CalcTemperatureComponentFunctor {
  using KineticEnergyDualViewType = Kokkos::DualView<double*, DeviceSpace::device_type, Kokkos::MemoryTraits<Kokkos::Atomic>>;
  using NumParticleDualViewType = Kokkos::DualView<size_t*, DeviceSpace::device_type, Kokkos::MemoryTraits<Kokkos::Atomic>>;

  KineticEnergyDualViewType _kineticEnergyMul2Map;
  NumParticleDualViewType _numParticleMap;

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, const autopas::utilsKokkos::KokkosStorage<ParticleType>& storage) const {
    const auto velX = storage.template operator()<ParticleType::AttributeNames::velocityX, ForEachHostFlag>(i);
    const auto velY = storage.template operator()<ParticleType::AttributeNames::velocityY, ForEachHostFlag>(i);
    const auto velZ = storage.template operator()<ParticleType::AttributeNames::velocityZ, ForEachHostFlag>(i);

    const auto mass = storage.template operator()<ParticleType::AttributeNames::mass, ForEachHostFlag>(i);
    const auto typeId = storage.template operator()<ParticleType::AttributeNames::typeId, ForEachHostFlag>(i);

    _kineticEnergyMul2Map.view_device()(typeId) += mass * (velX * velX + velY * velY + velZ * velZ);
    _numParticleMap.view_device()(typeId) += 1;
  }
};

template <class ParticleType>
struct BrownianMotionFunctor {
  using FloatPrecision = ParticleType::ParticleSoAFloatPrecision;
  using RandomPool = Kokkos::Random_XorShift64_Pool<>;

  RandomPool _randomEngine;
  double _targetTemperature;

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, const autopas::utilsKokkos::KokkosStorage<ParticleType>& storage) const {
    auto generator = _randomEngine.get_state();

    const auto mass = storage.template operator()<ParticleType::AttributeNames::mass, ForEachHostFlag>(i);
    double velScale = Kokkos::sqrt(_targetTemperature / mass);

    const FloatPrecision velIncrementX = generator.normal(0, 1) * velScale;
    const FloatPrecision velIncrementY = generator.normal(0, 1) * velScale;
    const FloatPrecision velIncrementZ = generator.normal(0, 1) * velScale;

    storage.template operator()<ParticleType::AttributeNames::velocityX, ForEachHostFlag>(i) += velIncrementX;
    storage.template operator()<ParticleType::AttributeNames::velocityY, ForEachHostFlag>(i) += velIncrementY;
    storage.template operator()<ParticleType::AttributeNames::velocityZ, ForEachHostFlag>(i) += velIncrementZ;

    _randomEngine.free_state(generator);
  }
};

template <class ParticleType>
struct ApplyFunctor {
  Kokkos::DualView<double*, DeviceSpace::device_type, Kokkos::MemoryTraits<Kokkos::Atomic>> _scalingMap;

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, const autopas::utilsKokkos::KokkosStorage<ParticleType>& storage) const {
    const auto typeId = storage.template operator()<ParticleType::AttributeNames::typeId, ForEachHostFlag>(i);
    const auto scaling = _scalingMap.view_device()(typeId);

    storage.template operator()<ParticleType::AttributeNames::velocityX, ForEachHostFlag>(i) *= scaling;
    storage.template operator()<ParticleType::AttributeNames::velocityY, ForEachHostFlag>(i) *= scaling;
    storage.template operator()<ParticleType::AttributeNames::velocityZ, ForEachHostFlag>(i) *= scaling;
  }
};

#endif

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
  bool containerAllowsKokkos = autopas.containerAllowsKokkos();
  if (!containerAllowsKokkos) {
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
  } else {
#ifdef AUTOPAS_ENABLE_KOKKOS
    CalcTemperatureFunctor<ParticleType> functor {};
    autopas.template reduceKokkos<DeviceSpace::execution_space, double, Kokkos::Sum<double>>(functor, kineticEnergyMul2, autopas::IteratorBehavior::ownedOrHalo,  "mdFlexible::Thermostat::calcTemperature"); // TODO: check iterator behavior
#endif
    // TODO: throw exception
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
 * @tparam AutoPasTemplate Type of AutoPas Object (no pointer)
 * @tparam ParticlePropertiesLibraryTemplate Type of ParticlePropertiesLibrary Object (no pointer)
 * @param autopas
 * @param particlePropertiesLibrary
 * @return map of: particle typeID -> temperature for this type
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
auto calcTemperatureComponent(AutoPasTemplate &autopas,
                              ParticlePropertiesLibraryTemplate &particlePropertiesLibrary) {
  using autopas::utils::ArrayMath::dot;
  using namespace autopas::utils::ArrayMath::literals;

  const auto numberComponents =
#if MD_FLEXIBLE_MODE == SINGLESITE
      particlePropertiesLibrary.getNumberRegisteredSiteTypes();
#elif MD_FLEXIBLE_MODE == MULTISITE
      particlePropertiesLibrary.getNumberRegisteredMolTypes();
#endif

  std::map<size_t, double> kineticEnergyMul2Map;
  std::map<size_t, size_t> numParticleMap;

  for (int typeID = 0; typeID < numberComponents; typeID++) {
    kineticEnergyMul2Map.at(typeID) = 0;
    numParticleMap.at(typeID) = 0;
  }

  bool containerAllowsKokkos = autopas.containerAllowsKokkos();
  if (!containerAllowsKokkos) {
    AUTOPAS_OPENMP(parallel) {
      // create aggregators for each thread
      std::map<size_t, double> kineticEnergyMul2MapThread;
      std::map<size_t, size_t> numParticleMapThread;
      for (int typeID = 0; typeID < numberComponents; typeID++) {
        kineticEnergyMul2MapThread[typeID] = 0.;
        numParticleMapThread[typeID] = 0ul;
      }
      // parallel iterators
      for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
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
        for (int typeID = 0; typeID < numberComponents; typeID++) {
          kineticEnergyMul2Map.at(typeID) += kineticEnergyMul2MapThread[typeID];
          numParticleMap.at(typeID) += numParticleMapThread[typeID];
        }
      }
    }
  } else {
#ifdef AUTOPAS_ENABLE_KOKKOS
    // TODO: Kokkos version of kinetic energy maps
    //CalcTemperatureComponentFunctor<ParticleType> functor {kineticEnergyMul2Map, numParticleMap};
    //autopas.template forEachKokkos<DeviceSpace::execution_space>(functor, autopas::IteratorBehavior::owned, "mdFlexible::Thermostat::calcTemperatureComponent"); // TODO: check which iterator behavior to use
    //kineticEnergyMul2Map.modify<DeviceSpace::execution_space>();
    //numParticleMap.modify<DeviceSpace::execution_space>();
#endif
    // TODO: throw exception
  }
  // md-flexible's molecules have 3 DoF for translational velocity and optionally 3 additional rotational DoF
#if MD_FLEXIBLE_MODE == MULTISITE
  constexpr unsigned int degreesOfFreedom{6};
#else
  constexpr unsigned int degreesOfFreedom{3};
#endif

    for (int typeID = 0; typeID < numberComponents; typeID++) {
      // workaround for MPICH: send and receive buffer must not be the same.
      autopas::AutoPas_MPI_Allreduce(AUTOPAS_MPI_IN_PLACE, &kineticEnergyMul2Map.at(typeID), 1, AUTOPAS_MPI_DOUBLE,
                                     AUTOPAS_MPI_SUM, AUTOPAS_MPI_COMM_WORLD);

      autopas::AutoPas_MPI_Allreduce(AUTOPAS_MPI_IN_PLACE, &numParticleMap.at(typeID), 1, AUTOPAS_MPI_UNSIGNED_LONG,
                                     AUTOPAS_MPI_SUM, AUTOPAS_MPI_COMM_WORLD);
    }

    for (int typeID = 0; typeID < numberComponents; typeID++) {
      kineticEnergyMul2Map.at(typeID) /= static_cast<double>(numParticleMap.at(typeID)) * degreesOfFreedom;
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
 * @param particlePropertiesLibrary
 * @param targetTemperature temperature of the system after applying the function on a system with temperature = 0.
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void addBrownianMotion(AutoPasTemplate &autopas, ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
                       const double targetTemperature) {
  using namespace autopas::utils::ArrayMath::literals;
  // Generate map(s) of molecule type Id to scaling factors
  std::map<size_t, double> translationalVelocityScale;
  std::map<size_t, std::array<double, 3>> rotationalVelocityScale;

  const auto numberComponents =
#if MD_FLEXIBLE_MODE == SINGLESITE
      particlePropertiesLibrary.getNumberRegisteredSiteTypes();
#elif MD_FLEXIBLE_MODE == MULTISITE
      particlePropertiesLibrary.getNumberRegisteredMolTypes();
#endif

  for (int typeID = 0; typeID < numberComponents; typeID++) {
    translationalVelocityScale.emplace(typeID,
                                       std::sqrt(targetTemperature / particlePropertiesLibrary.getMolMass(typeID)));
#if MD_FLEXIBLE_MODE == MULTISITE
    const auto momentOfInertia = particlePropertiesLibrary.getMomentOfInertia(typeID);
    const std::array<double, 3> scale{std::sqrt(targetTemperature / momentOfInertia[0]),
                                      std::sqrt(targetTemperature / momentOfInertia[1]),
                                      std::sqrt(targetTemperature / momentOfInertia[2])};
    rotationalVelocityScale.emplace(typeID, scale);
#endif
  }
  bool containerAllowsKokkos = autopas.containerAllowsKokkos();
  if (!containerAllowsKokkos) {
    AUTOPAS_OPENMP(parallel default(none) shared(autopas, translationalVelocityScale, rotationalVelocityScale)) {
      // we use a constant seed for repeatability.
      // we need one random engine and distribution per thread
      std::default_random_engine randomEngine(42 + autopas::autopas_get_thread_num());
      std::normal_distribution<double> normalDistribution{0, 1};
      for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
        const std::array<double, 3> normal3DVecTranslational = {
          normalDistribution(randomEngine), normalDistribution(randomEngine), normalDistribution(randomEngine)};
        auto velIncrement = normal3DVecTranslational * translationalVelocityScale[iter->getTypeId()];
        iter->addV({
          static_cast<ParticleType::ParticleSoAFloatPrecision>(velIncrement.at(0)),
          static_cast<ParticleType::ParticleSoAFloatPrecision>(velIncrement.at(1)),
          static_cast<ParticleType::ParticleSoAFloatPrecision>(velIncrement.at(2))
        });
#if MD_FLEXIBLE_MODE == MULTISITE
        const std::array<double, 3> normal3DVecRotational = {
          normalDistribution(randomEngine), normalDistribution(randomEngine), normalDistribution(randomEngine)};
        iter->addAngularVel(normal3DVecRotational * rotationalVelocityScale[iter->getTypeId()]);
#endif
      }
    }
  } else {
#ifdef AUTOPAS_ENABLE_KOKKOS
    Kokkos::Random_XorShift64_Pool<> random_engine (42);
    BrownianMotionFunctor<ParticleType> functor(random_engine, targetTemperature);
    autopas.template forEachKokkos<DeviceSpace::execution_space>(functor, autopas::IteratorBehavior::ownedOrHalo, "mdFlexible::Thermostat::addBrownianMotion"); // TODO: check iterator behavior
#endif
    // TODO: throw exception
  }
}

/**
 * Scales velocity of particles towards a given temperature. For Multi-site simulations, angular velocity is also
 * scaled.
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
  using namespace autopas::utils::ArrayMath::literals;

  AutoPasLog(DEBUG, "Applying Thermostat");

  const auto currentTemperatureMap = calcTemperatureComponent(autopas, particlePropertiesLibrary);

  // make sure we work with a positive delta
  const double absoluteDeltaTemperature = std::abs(deltaTemperature);

  std::remove_const_t<decltype(currentTemperatureMap)> scalingMap;

  for (int typeId = 0; typeId < currentTemperatureMap.size(); typeId++) {
    // If the current temperature is within absoluteDeltaTemperature of target temperature,
    //      set immediate target temperature to target temperature.
    //
    // Else, set immediate target temperature to absoluteDeltaTemperature towards the target temperature from the
    // current temperature.

    double currentTemperature = currentTemperatureMap.at(typeId);

    const auto immediateTargetTemperature =
        currentTemperature < targetTemperature
            ? std::min(currentTemperature + absoluteDeltaTemperature, targetTemperature)
            : std::max(currentTemperature - absoluteDeltaTemperature, targetTemperature);
    // Determine a scaling factor for each particle type.
    scalingMap.at(typeId) = std::sqrt(immediateTargetTemperature / currentTemperature);

    AutoPasLog(DEBUG, "Current temperature of typeID {}: {}", typeId, currentTemperature);
    AutoPasLog(DEBUG, "Temperature of typeID {} after application of thermostat: {}", typeId,
               immediateTargetTemperature);
  }

  if (!autopas.containerAllowsKokkos()) {
    // Scale velocities (and angular velocities) with the scaling map
    AUTOPAS_OPENMP(parallel default(none) shared(autopas, scalingMap))
    for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
      std::array newVel = iter->getV() * scalingMap.at(iter->getTypeId());
      iter->setV(newVel);
#if MD_FLEXIBLE_MODE == MULTISITE
      iter->setAngularVel(iter->getAngularVel() * scalingMap[iter->getTypeId()]);
#endif
    }
  } else {
#ifdef AUTOPAS_ENABLE_KOKKOS
    // TODO: Kokkos version of scalingMap
    // scalingMap.template sync<DeviceSpace::execution_space>();
    // ApplyFunctor<ParticleType> functor {scalingMap};
    // autopas.template forEachKokkos<DeviceSpace::execution_space>(functor, autopas::IteratorBehavior::owned, "mdFlexible::Thermostat::apply"); // TODO: decide iterator behavior, figure out how to handle scalingMap
#endif
    // TODO: throw exception
  }
}
}  // namespace Thermostat
