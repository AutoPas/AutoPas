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
                        const std::array<double, 3> &globalForce) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    auto v = iter->getV();
    auto m = particlePropertiesLibrary.getMass(iter->getTypeId());
    auto f = iter->getF();
    iter->setOldF(f);
    iter->setF(globalForce);
    v = mulScalar(v, deltaT);
    f = mulScalar(f, (deltaT * deltaT / (2 * m)));
    auto newR = add(v, f);
    iter->addR(newR);
  }
}

/**
 * Calculate and update the quaternion for every particle. Throws error unless ParticleClass is specialised to a
 * rotational molecule, i.e. MulticenteredMoleculeLJ.
 * @tparam ParticleClass
 * @param autoPasContainer
 * @param particlePropertiesLibrary
 * @param deltaT
 * @param globalForce
 */
template <class ParticleClass>
void calculateQuaternions(autopas::AutoPas<ParticleClass> &autoPasContainer,
                          const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                          const std::array<double, 3> &globalForce) {
  autopas::utils::ExceptionHandler::exception("calculateQuaternion should not be run with a non-rotational molecule type!");
}

/**
 * Calculate and update the quaternion for every particle. Uses the rotational velocity-verlet algorithm as described by
 * Rozmanov, 2010, Robust rotational-velocity-Verlet integration methods (method A); with slight adaptations to account
 * for md-flexible primarily using (angular) velocities rather than (angular) momentums. Code lines are commented with
 * references to corresponding equations within the paper.
 * @param autoPasContainer
 * @param particlePropertiesLibrary
 * @param deltaT
 * @param globalForce
 */
template<> void calculateQuaternions<MulticenteredMoleculeLJ>(autopas::AutoPas<MulticenteredMoleculeLJ> &autoPasContainer,
                                                   const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                                                   const std::array<double, 3> &globalForce) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::addScalar;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayMath::mul;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::cross;
  using autopas::utils::ArrayMath::L2Norm;
  using autopas::utils::ArrayMath::normalize;
  using autopas::utils::quaternion::qMul;
  using autopas::utils::quaternion::rotatePositionBackwards;

  const auto halfDeltaT = 0.5 * deltaT;

  const double tol = 1e-13; // tolerance given in paper

  // todo sort out how to handle global forces.

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto q = iter->getQ();
    const auto angVelW = iter->getAngularVel(); // angular velocity in world frame
    const auto angVelM = rotatePositionBackwards(q,angVelW); // angular velocity in molecular frame  (equivalent to (17))
    const auto torqueW = iter->getTorque();
    const auto torqueM = rotatePositionBackwards(q,torqueW); // (18)

    const auto I = particlePropertiesLibrary.getMomentOfInertia(iter->getTypeId()); // moment of inertia

    const auto angMomentumM = mul(I,angVelM); // equivalent to (19)
    const auto derivativeAngMomentumM = sub(torqueM, cross(angVelM,angMomentumM)); // (20)
    const auto angMomentumMHalfStep = add(angMomentumM, mulScalar(derivativeAngMomentumM, halfDeltaT)); // (21)

    auto derivativeQHalfStep = mulScalar(qMul(q, div(angMomentumMHalfStep, I)), 0.5); // (22)

    auto qHalfStep = normalize(add(q, mulScalar(derivativeQHalfStep,halfDeltaT))); // (23)

    const auto angVelWHalfStep = add(angVelW, mulScalar(div(torqueW, I), halfDeltaT)); // equivalent to (24)

    // (25) start
    // initialise qHalfStepOld to be outside tolerable distance from qHalfStep to satisfy while statement
    auto qHalfStepOld = qHalfStep;
    qHalfStepOld[0] += 2 * tol;

    while (L2Norm(sub(qHalfStep,qHalfStepOld))>tol) {
      qHalfStepOld = qHalfStep;
      auto angVelMHalfStep = rotatePositionBackwards(qHalfStepOld,angVelWHalfStep); // equivalent to first two lines of (25)
      derivativeQHalfStep = mulScalar(qMul(qHalfStepOld,angVelWHalfStep),0.5);
      qHalfStep = normalize(add(q, mulScalar(derivativeQHalfStep, halfDeltaT)));
    }
    // (25) end

    const auto qFullStep = normalize(add(q, mulScalar(derivativeQHalfStep, deltaT)));

    iter->setQ(qFullStep);
    iter->setAngularVel(angVelWHalfStep); // save angular velocity half step, to be used by calculateAngularVelocities
  }
}

/**
 * Calculate and update the velocity for every particle using the the Störmer-Verlet Algorithm.
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
    auto m = particlePropertiesLibrary.getMass(iter->getTypeId());
    auto force = iter->getF();
    auto oldForce = iter->getOldF();
    auto newV = mulScalar((add(force, oldForce)), deltaT / (2 * m));
    iter->addV(newV);
  }
}

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
                                const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  autopas::utils::ExceptionHandler::exception("calculateAngularVelocities should not be run with a non-rotational molecule type!");
}

/**
 * Calculate and update the angular velocity for every particle. Uses the rotational velocity-verlet algorithm as
 * described by Rozmanov, 2010, Robust rotational-velocity-Verlet integration methods (method A); with slight adaptations to account
 * for md-flexible primarily using (angular) velocities rather than (angular) momentums. Code lines are commented with
 * references to corresponding equations within the paper.
 * @param autoPasContainer
 * @param particlePropertiesLibrary
 * @param deltaT
 * @param globalForce
 */
template<> void calculateAngularVelocities<MulticenteredMoleculeLJ>(autopas::AutoPas<MulticenteredMoleculeLJ> &autoPasContainer,
                                                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::div;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto torque = iter->getTorque();
    const auto I = particlePropertiesLibrary.getMomentOfInertia(iter->getTypeId()); // moment of inertia
    iter->addAngularVel(mulScalar(div(torque, I), 0.5*deltaT)); // (28)
  }
}

}  // namespace TimeDiscretization
