/**
 * @file TimeDiscretization.h
 * @author J. Körner
 * @date 20/05/21
 */
#pragma once

#include "TypeDefinitions.h"
#include "autopas/AutoPas.h"
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
 */
void calculatePositions(autopas::AutoPas<ParticleType> &autoPasContainer,
                        const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT);

/**
 * Calculate and update the velocity for every particle using the the Störmer-Verlet Algorithm.
 * @param autoPasContainer The container for which to update the velocities.
 * @param particlePropertiesLibrary The particle properties library for the particles in the container.
 * @param deltaT The time step width.
 * @param useThermmostat Decides if a thermostat should be used.
 * @param targetTemperature The target temperature used for the thermostat. If useThermostat is false,
 * 	this value does not matter.
 */
void calculateVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                         const bool &useThermostat, const double &targetTemperature);

/**
 * Calculates the pairwise forces between particles in an autopas container.
 * @param autoPasContainer The container for which to update the pairwise forces.
 * @param particlePropertiesLibrary The particle properties library for the particles in the container.
 * @param deltaT The time step width.
 *	@param funtorOption The functor on which to base the force calculations.
 * @param wasTuninIteration Tells the user if the current iteration of force calculations was a tuning iteration.
 */
void calculatePairwiseForces(autopas::AutoPas<ParticleType> &autoPasContainer,
                             ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                             MDFlexConfig::FunctorOption functorOption, bool &wasTuningIteration);

/**
 * Adds global forces to the particles in the container.
 * @param autoPasContainer The container for which to update the particle forces.
 * @param globalForce The global force which will be applied to each particle in the container.
 */
void calculateGlobalForces(autopas::AutoPas<ParticleType> &autoPasContainer, std::array<double, 3> &globalForce);
}  // namespace TimeDiscretization
