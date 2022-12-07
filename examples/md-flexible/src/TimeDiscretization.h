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

#include "TimeDiscretization.cpp"

/**
 * Functions for updating velocities and positions as simulation time progresses.
 */
namespace TimeDiscretization {
/**
 * Calculate and update the position for every particle using the Störmer-Verlet Algorithm.
 * In addition, pushes the force stored in the force vector to the old force vector and sets the force vector to the
 * global force in preparation for the calculate forces stage.
 *
 * Specifically, the formula for this is
 *      x_{n+1} = x_n + delta_t * v_n + delta_t^2 / ( 2 * mass) * f_n
 *                      {   velTerm }   {        forceTerm          }
 *
 * @param autoPasContainer The container for which to update the positions.
 * @param particlePropertiesLibrary The particle properties library for the particles in the container.
 * @param deltaT The time step width.
 * @param globalForce Base force value to which every particle is reset.
 * @param fastParticlesThrow When true throws an exception if particles moved too far for verlet technique
 * (>skin/2/rebuildFrequency).
 */
void calculatePositionsAndUpdateForces(autopas::AutoPas<ParticleType> &autoPasContainer,
                                       const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                                       const std::array<double, 3> &globalForce, bool fastParticlesThrow)

/**
 * Calculate and update the quaternion for every particle. Uses the rotational velocity-verlet algorithm as described by
 * Rozmanov, 2010, Robust rotational-velocity-Verlet integration methods (method A); with slight adaptations to account
 * for md-flexible primarily using (angular) velocities rather than (angular) momentums. Code lines are commented with
 * references to corresponding equations within the paper.
 *
 * @note Throws error if md-flexible is compiled without multi-site support.
 *
 * @param autoPasContainer
 * @param particlePropertiesLibrary
 * @param deltaT
 * @param globalForce
 */
inline void calculateQuaternions(autopas::AutoPas<ParticleType> &autoPasContainer,
                          const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                          const std::array<double, 3> &globalForce) {}

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
template <class ParticleClass>
void calculateVelocities(autopas::AutoPas<ParticleClass> &autoPasContainer,
                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  // helper declarations for operations with vector
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto molecularMass = particlePropertiesLibrary.getMolMass(iter->getTypeId());
    const auto force = iter->getF();
    const auto oldForce = iter->getOldF();
    const auto changeInVel = mulScalar((add(force, oldForce)), deltaT / (2 * molecularMass));
    iter->addV(changeInVel);
  }
}

/**
 * Calculate and update the angular velocity for every particle. Throws error unless ParticleClass is specialised to a
 * rotational molecule, i.e. MultisiteMoleculeLJ.
 * @tparam ParticleClass
 * @param autoPasContainer
 * @param particlePropertiesLibrary
 * @param deltaT
 * @param globalForce
 */
template <class ParticleClass>
inline void calculateAngularVelocities(autopas::AutoPas<ParticleClass> &autoPasContainer,
                                       const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  autopas::utils::ExceptionHandler::exception("calculateAngularVelocities should not be run with a non-rotational molecule type!");
}

template<> inline void calculateAngularVelocities<autopas::MultisiteMoleculeLJ>(autopas::AutoPas<autopas::MultisiteMoleculeLJ> &autoPasContainer,
                                                                     const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::quaternion::rotatePosition;
  using autopas::utils::quaternion::rotatePositionBackwards;

  //#ifdef AUTOPAS_OPENMP
  //#pragma omp parallel
  //#endif
  // for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
  //   const auto torqueW = iter->getTorque();
  //   const auto q = iter->getQ();
  //   const auto I = particlePropertiesLibrary.getMomentOfInertia(iter->getTypeId()); // moment of inertia
  //
  //   // convert torque to molecular-frame
  //   const auto torqueM = rotatePositionBackwards(q, torqueW);
  //
  //   // get I^-1 T in molecular-frame
  //   const auto torqueDivMoIM = div(torqueM, I);
  //
  //   // convert to world-frame
  //   const auto torqueDivMoIW = rotatePosition(q, torqueDivMoIM);
  //
  //   iter->addAngularVel(mulScalar(torqueDivMoIW, 0.5*deltaT)); // (28)
  // }
}


}  // namespace TimeDiscretization
