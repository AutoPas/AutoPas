/**
 * @file TimeDiscretization.h
 * @author N. Fottner
 * @date 13/05/19
 */
#pragma once

#include "TypeDefinitions.h"
#include "autopas/AutoPasDecl.h"
#include "src/configuration/MDFlexConfig.h"

/**
 * Functions for updating velocities and positions as simulation time progresses.
 */
namespace TimeDiscretization {
/**
 * Calculate and update the position for every particle using the Störmer-Verlet Algorithm.
 * @param autoPasContainer The container for which to update the positions.
 * @param particlePropertiesLibrary The particle properties library for the particles in the container.
 * @param deltaT The time step width.
 * @param globalForce Base force value to which every particle is reset.
 */
template <class ParticleClass>
void calculatePositions(autopas::AutoPas<ParticleClass> &autoPasContainer,
                        const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                        const std::array<double, 3> &globalForce);

/**
 * Calculate and update the quaternion for every particle
 * @tparam ParticleClass
 * @param autoPasContainer
 * @param particlePropertiesLibrary
 * @param deltaT
 * @param globalForce
 */
template <class ParticleClass>
void calculateQuaternions(autopas::AutoPas<ParticleClass> &autoPasContainer,
                          const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                          const std::array<double, 3> &globalForce);

/**
 * Calculate and update the velocity for every particle using the the Störmer-Verlet Algorithm.
 * @param autoPasContainer The container for which to update the velocities.
 * @param particlePropertiesLibrary The particle properties library for the particles in the container.
 * @param deltaT The time step width.
 */
template <class ParticleClass>
void calculateVelocities(autopas::AutoPas<ParticleClass> &autoPasContainer,
                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT);

}  // namespace TimeDiscretization
