/**
 * @file TimeDiscretization.h
 * @author N. Fottner
 * @date 13/05/19
 */
#pragma once
#include <chrono>
#include <fstream>
#include <iostream>
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/MemoryProfiler.h"

template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
class TimeDiscretization {
 public:
  using TimeDiscretizationFloatType = typename ParticlePropertiesLibraryTemplate::ParticlePropertiesLibraryFloatType;
  using TimeDiscretizationIntType = typename ParticlePropertiesLibraryTemplate::ParticlePropertiesLibraryIntType;

  /**
   * Constructor
   * */
  explicit TimeDiscretization(TimeDiscretizationFloatType particleDeltaT,
                              ParticlePropertiesLibraryTemplate &particlePropertiesLibrary)
      : particle_delta_t(particleDeltaT), _particlePropertiesLibrary(particlePropertiesLibrary){};

  /**
   * Default Destructor
   * */
  virtual ~TimeDiscretization() = default;

  /**
   * Calculate the new Position for every Praticle using the Iterator and the Störmer-Verlet Algorithm
   * @param autopas
   * @return time for the calculation in microseconds
   */
  long CalculateX(AutoPasTemplate &autopas);

  /**
   * Calculate the new Velocity for every Praticle using the Iterator and the Störmer-Verlet Algorithm
   * @param autopas
   * @return time for the calculation in microseconds
   */
  long CalculateV(AutoPasTemplate &autopas);
  /**
   * Getter for particleDeltaT
   * @return particle_delta_t
   */
  TimeDiscretizationFloatType getParticleDeltaT() const;

 private:
  /**
   * Duration of a timestep
   */
  TimeDiscretizationFloatType particle_delta_t;
  /**
   * ParticlePropertiesLibrary to access Mass of Particles
   */
  ParticlePropertiesLibraryTemplate _particlePropertiesLibrary;
};

template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
long TimeDiscretization<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::CalculateX(AutoPasTemplate &autopas) {
  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
  startCalc = std::chrono::high_resolution_clock::now();
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto v = iter->getV();
    auto m = _particlePropertiesLibrary.getMass(iter->getTypeId());
    auto f = iter->getF();
    iter->setOldF(f);
    iter->setF({0., 0., 0.});
    v = autopas::ArrayMath::mulScalar(v, particle_delta_t);
    f = autopas::ArrayMath::mulScalar(f, (particle_delta_t * particle_delta_t / (2 * m)));
    auto newR = autopas::ArrayMath::add(v, f);
    iter->addR(newR);
  }
  stopCalc = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
}

template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
long TimeDiscretization<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::CalculateV(AutoPasTemplate &autopas) {
  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
  startCalc = std::chrono::high_resolution_clock::now();
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(autopas::IteratorBehavior::ownedOnly); iter.isValid(); ++iter) {
    auto m = _particlePropertiesLibrary.getMass(iter->getTypeId());
    auto force = iter->getF();
    auto old_force = iter->getOldf();
    auto newV = autopas::ArrayMath::mulScalar((autopas::ArrayMath::add(force, old_force)), particle_delta_t / (2 * m));
    iter->addV(newV);
  }
  stopCalc = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
}

template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
typename TimeDiscretization<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::TimeDiscretizationFloatType
TimeDiscretization<AutoPasTemplate, ParticlePropertiesLibraryTemplate>::getParticleDeltaT() const {
  return particle_delta_t;
}