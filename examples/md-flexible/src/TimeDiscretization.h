/**
 * @file TimeDiscretization.h
 * @author N. Fottner & S. Newcome
 * @date 13/05/19
 */
#pragma once

#include "TypeDefinitions.h"
#include "autopas/AutoPasDecl.h"
#include "src/configuration/MDFlexConfig.h"
#include "autopas/utils/Quaternion.h"

#include "TimeDiscretizationImpl.h"

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
 * Calculate and update the quaternion for every particle. Uses the rotational velocity-verlet algorithm as described by
 * Rozmanov, 2010, Robust rotational-velocity-Verlet integration methods (method A); with slight adaptations to account
 * for md-flexible primarily using (angular) velocities rather than (angular) momentums. Code lines are commented with
 * references to corresponding equations within the paper. Throws error unless ParticleClass is specialised to a
* rotational molecule, i.e. MulticenteredMoleculeLJ.
 * @tparam ParticleClass
 * @param autoPasContainer
 * @param particlePropertiesLibrary
 * @param deltaT
 * @param globalForce
 */
//template <class ParticleClass>
//void calculateQuaternions(autopas::AutoPas<ParticleClass> &autoPasContainer,
//                          const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
//                          const std::array<double, 3> &globalForce);

/**
 * Calculate and update the velocity for every particle using the the Störmer-Verlet Algorithm.
 * @param autoPasContainer The container for which to update the velocities.
 * @param particlePropertiesLibrary The particle properties library for the particles in the container.
 * @param deltaT The time step width.
 */
template <class ParticleClass>
void calculateVelocities(autopas::AutoPas<ParticleClass> &autoPasContainer,
                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT);

/**
 * Calculate and update the angular velocity for every particle. Throws error unless ParticleClass is specialised to a
 * rotational molecule, i.e. MulticenteredMoleculeLJ.
 * @tparam ParticleClass
 * @param autoPasContainer
 * @param particlePropertiesLibrary
 * @param deltaT
 * @param globalForce
 */
template <class ParticleClass>
void calculateAngularVelocities(autopas::AutoPas<ParticleClass> &autoPasContainer,
                                const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT);


}  // namespace TimeDiscretization
