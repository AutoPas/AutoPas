/**
 * @file TimeDiscretization.h
 * @author N. Fottner
 * @date 13/05/19
 */
#pragma once
#include "autopas/AutoPas.h"
#include "autopas/utils/ArrayMath.h"

/**
 * Functions for updating velocities and positions as simulation time progresses.
 */
namespace TimeDiscretization {

/**
 * Calculate and update the position for every particle using the Störmer-Verlet Algorithm.
 * @param autopas
 * @param particlePropertiesLibrary
 * @param deltaT time step width
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void calculatePositions(AutoPasTemplate &autopas, const ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
                        const double deltaT) {
  // helper declarations for operations with vector
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  std::cout << "Taddel 0" << std::endl;
  for (auto iter = autopas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
  	std::cout << "Krabs 0" << std::endl;
    auto v = iter->getV();
  	std::cout << "Krabs 1" << std::endl;
    auto m = particlePropertiesLibrary.getMass(iter->getTypeId());
  	std::cout << "Krabs 2" << std::endl;
    auto f = iter->getF();
  	std::cout << "Krabs 3" << std::endl;
    iter->setOldF(f);
  	std::cout << "Krabs 4" << std::endl;
    iter->setF({0., 0., 0.});
  	std::cout << "Krabs 5" << std::endl;
    v = mulScalar(v, deltaT);
  	std::cout << "Krabs 6" << std::endl;
    f = mulScalar(f, (deltaT * deltaT / (2 * m)));
  	std::cout << "Krabs 7" << std::endl;
    auto newR = add(v, f);
  	std::cout << "Krabs 8" << std::endl;
    iter->addR(newR);
  	std::cout << "Krabs 9" << std::endl;
  }
  std::cout << "Taddel 1" << std::endl;
}

/**
 * Calculate and update the velocity for every particle using the the Störmer-Verlet Algorithm.
 * @param autopas
 * @param particlePropertiesLibrary
 * @param deltaT time step width
 */
template <class AutoPasTemplate, class ParticlePropertiesLibraryTemplate>
void calculateVelocities(AutoPasTemplate &autopas, const ParticlePropertiesLibraryTemplate &particlePropertiesLibrary,
                         const double deltaT) {
  // helper declarations for operations with vector
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autopas.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    auto m = particlePropertiesLibrary.getMass(iter->getTypeId());
    auto force = iter->getF();
    auto oldForce = iter->getOldf();
    auto newV = mulScalar((add(force, oldForce)), deltaT / (2 * m));
    iter->addV(newV);
  }
}

}  // namespace TimeDiscretization
