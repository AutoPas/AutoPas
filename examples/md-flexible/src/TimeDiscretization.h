/**
 * @file TimeDiscretization.h
 * @author N. Fottner
 * @date 13/05/19
 */
#pragma once

#include "TypeDefinitions.h"
#include "autopas/AutoPasDecl.h"
#include "autopas/utils/Quaternion.h"
#include "src/configuration/MDFlexConfig.h"

/**
 * Functions for updating velocities and positions as simulation time progresses.
 */
namespace TimeDiscretization {
/**
 * Calculate and update the position for every particle using the Störmer-Verlet Algorithm.
 *
 * Specifically, the formula for this is
 *      x_{n+1} = x_n + delta_t * v_n + delta_t^2 / ( 2 * mass) * f_n
 *                      {   velTerm }   {        forceTerm          }
 *
 * In addition, pushes the force stored in the force vector to the old force vector and sets the force vector to the
 * global force in preparation for the calculate forces stage.
 *
 * When using the Bundling Molecule approach this function does NOT update the position of the sites based on the molecule movement.
 * Therefore this still has to be done lateron. The most reasonable place to do this is after the computation of the new rotation.
 *
 * @param autoPasContainer The container for which to update the positions.
 * @param particlePropertiesLibrary The particle properties library for the particles in the container.
 * @param deltaT The time step width.
 * @param globalForce Base force value to which every particle is reset.
 * @param fastParticlesThrow When true throws an exception if particles moved too far for verlet technique
 * (>skin/2/rebuildFrequency).
 */
#if not defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH) or MD_FLEXIBLE_MODE!=MULTISITE
void calculatePositionsAndResetForces(autopas::AutoPas<ParticleType> &autoPasContainer,
                                      const ParticlePropertiesLibraryType &particlePropertiesLibrary,
                                      const double &deltaT, const std::array<double, 3> &globalForce,
                                      bool fastParticlesThrow);
#else
void calculatePositionsAndResetForces(autopas::AutoPas<ParticleType> &autoPasContainer, MoleculeContainer& moleculeContainer,
                                      const ParticlePropertiesLibraryType &particlePropertiesLibrary,
                                      const double &deltaT, const std::array<double, 3> &globalForce,
                                      bool fastParticlesThrow);
#endif

/**
 * Calculate and update the quaternion for every particle. Uses the rotational velocity-verlet algorithm as described by
 * Rozmanov, 2010, Robust rotational-velocity-Verlet integration methods (https://doi.org/10.1103/PhysRevE.81.056706)
 * (method A); with slight adaptations to account for md-flexible primarily using (angular) velocities rather than
 * (angular) momentums. Code lines are commented with references to corresponding equations within the paper.
 *
 * In addition, resets the torques to that determined by the global force only.
 *
 * @note Throws error if md-flexible is compiled without multi-site support.
 *
 * @param autoPasContainer
 * @param particlePropertiesLibrary
 * @param deltaT
 * @param globalForce
 */
#if not defined MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH or MD_FLEXIBLE_MODE!=MULTISITE
void calculateQuaternionsAndResetTorques(autopas::AutoPas<ParticleType> &autoPasContainer,
                                         const ParticlePropertiesLibraryType &particlePropertiesLibrary,
                                         const double &deltaT, const std::array<double, 3> &globalForce);
#else
void calculateQuaternionsAndResetTorques(autopas::AutoPas<ParticleType> &autoPasContainer, MoleculeContainer& moleculeContainer,
                                         const ParticlePropertiesLibraryType &particlePropertiesLibrary,
                                         const double &deltaT, const std::array<double, 3> &globalForce);
#endif

/**
 * Calculate and update the velocity for every particle using the the Störmer-Verlet Algorithm.
 *
 * Specifically
 *      v_{n+1} = v_n + delta_t / (2 * mass) * (F_n + F_{n-1})
 *
 * @param autoPasContainer The container for which to update the velocities.
 * @param particlePropertiesLibrary The particle properties library for the particles in the container.
 * @param deltaT The time step width.
 */
#if not defined(MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH) or MD_FLEXIBLE_MODE!=MULTISITE
void calculateVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT);
#else
void calculateVelocities(autopas::AutoPas<ParticleType> &autoPasContainer, MoleculeContainer& moleculeContainer,
                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT);
#endif

/**
 * Calculate and update the angular velocity for every particle.
 *
 * @note Throws error if md-flexible is compiled without multi-site support.
 *
 * @param autoPasContainer
 * @param particlePropertiesLibrary
 * @param deltaT
 */
#if not defined MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH or MD_FLEXIBLE_MODE!=MULTISITE
void calculateAngularVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
                                const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT);
#else
void calculateAngularVelocities(autopas::AutoPas<ParticleType> &autoPasContainer, MoleculeContainer& moleculeContainer,
                                const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT);
#endif

}  // namespace TimeDiscretization
