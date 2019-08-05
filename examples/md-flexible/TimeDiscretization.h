//
// Created by nicola on 13.05.19.
//

#pragma once
#include <chrono>
#include <fstream>
#include <iostream>
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/MemoryProfiler.h"

template <class AutoPasTemplate, typename floatType, typename intType>
class TimeDiscretization {
 public:
  explicit TimeDiscretization(floatType particleDeltaT,
                              ParticlePropertiesLibrary<floatType, intType> &particlePropertiesLibrary);

  virtual ~TimeDiscretization() = default;

  /**
   * Calculate the new Position for every Praticle using the Iterator and the Störmer-Verlet Algorithm
   * @param autopas
   * @return time for the calculation in microseconds
   */
  long VSCalculateX(AutoPasTemplate &autopas);

  /**
   * Calculate the new Velocity for every Praticle using the Iterator and the Störmer-Verlet Algorithm
   * @param autopas
   * @return time for the calculation in microseconds
   */
  long VSCalculateV(AutoPasTemplate &autopas);
  /**Getter for particleDeltaT
   * @return particle_delta_t
   * */
  floatType getParticleDeltaT() const;

 private:
  /**  Duration of a timestep
   * */
  floatType particle_delta_t;
  ParticlePropertiesLibrary<floatType, intType> _particlePropertiesLibrary;
};

template <class AutoPasTemplate, typename floatType, typename intType>
long TimeDiscretization<AutoPasTemplate, floatType, intType>::VSCalculateX(AutoPasTemplate &autopas) {
  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
  startCalc = std::chrono::high_resolution_clock::now();
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    auto v = iter->getV();
    auto m = _particlePropertiesLibrary.getMass(iter->getTypeId());
    auto f = iter->getF();
    iter->setOldf(f);
    v = autopas::ArrayMath::mulScalar(v, particle_delta_t);
    f = autopas::ArrayMath::mulScalar(f, (particle_delta_t * particle_delta_t / (2 * m)));
    auto newR = autopas::ArrayMath::add(v, f);
    iter->addR(newR);
  }
  stopCalc = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
}

template <class AutoPasTemplate, typename floatType, typename intType>
long TimeDiscretization<AutoPasTemplate, floatType, intType>::VSCalculateV(AutoPasTemplate &autopas) {
  std::chrono::high_resolution_clock::time_point startCalc, stopCalc;
  startCalc = std::chrono::high_resolution_clock::now();
#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(); iter.isValid(); ++iter) {
    auto m = _particlePropertiesLibrary.getMass(iter->getTypeId());
    auto force = iter->getF();
    auto old_force = iter->getOldf();
    auto newV = autopas::ArrayMath::mulScalar((autopas::ArrayMath::add(force, old_force)), particle_delta_t / (2 * m));
    iter->addV(newV);
  }
  stopCalc = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(stopCalc - startCalc).count();
}

template <class AutoPasTemplate, typename floatType, typename intType>
TimeDiscretization<AutoPasTemplate, floatType, intType>::TimeDiscretization(
    floatType particleDeltaT, ParticlePropertiesLibrary<floatType, intType> &particlePropertiesLibrary)
    : particle_delta_t(particleDeltaT), _particlePropertiesLibrary(particlePropertiesLibrary) {}

template <class AutoPasTemplate, typename floatType, typename intType>
floatType TimeDiscretization<AutoPasTemplate, floatType, intType>::getParticleDeltaT() const {
  return particle_delta_t;
}