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

  // TODO: SWiMM
  // FIXME: Calculate the new positions for the owned particles.
  //  Hint: You may use the Störmer-Verlet approach to update the positions (r) where:
  //        r_new = r + v * deltaT + deltaT^2 * f/(2*m)
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

  // TODO: SWiMM
  // FIXME: Calculate the new velocities for the owned particles.
  //  Hint: You may use the Störmer-Verlet approach to update the positions (r) where:
  //        v_new = v + deltaT * (f_old + f) / (2*m)
}

};  // namespace TimeDiscretization
