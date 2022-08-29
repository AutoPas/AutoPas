/**
* @file TimeDiscretization.cpp
* @author N. Fottner & S. Newcome
* @date 13/05/19
*/

#include "TimeDiscretization.h"

namespace TimeDiscretization {

template <> void calculatePositions<autopas::MoleculeLJ>(autopas::AutoPas<autopas::MoleculeLJ> &autoPasContainer,
                       const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                       const std::array<double, 3> &globalForce);

template <> void calculatePositions<autopas::MulticenteredMoleculeLJ>(autopas::AutoPas<autopas::MulticenteredMoleculeLJ> &autoPasContainer,
                                             const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                                             const std::array<double, 3> &globalForce);

template <class ParticleClass>
void calculateQuaternions(autopas::AutoPas<ParticleClass> &autoPasContainer,
                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                         const std::array<double, 3> &globalForce) {
 autopas::utils::ExceptionHandler::exception("calculateQuaternion should not be run with a non-rotational molecule type!");
}


template<> void calculateQuaternions<autopas::MulticenteredMoleculeLJ>(autopas::AutoPas<autopas::MulticenteredMoleculeLJ> &autoPasContainer,
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
 using autopas::utils::quaternion::rotatePosition;
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

   const auto angVelWHalfStep = add(angVelW, mulScalar(rotatePosition(q,div(torqueM, I)), halfDeltaT)); // equivalent to (24)

   // (25) start
   // initialise qHalfStepOld to be outside tolerable distance from qHalfStep to satisfy while statement
   auto qHalfStepOld = qHalfStep;
   qHalfStepOld[0] += 2 * tol;

   while (L2Norm(sub(qHalfStep,qHalfStepOld))>tol) {
     qHalfStepOld = qHalfStep;
     auto angVelMHalfStep = rotatePositionBackwards(qHalfStepOld,angVelWHalfStep); // equivalent to first two lines of (25)
     derivativeQHalfStep = mulScalar(qMul(qHalfStepOld,angVelMHalfStep),0.5);
     qHalfStep = normalize(add(q, mulScalar(derivativeQHalfStep, halfDeltaT)));
   }
   // (25) end

   const auto qFullStep = normalize(add(q, mulScalar(derivativeQHalfStep, deltaT))); // (26)

   iter->setQ(qFullStep);
   iter->setAngularVel(angVelWHalfStep); // save angular velocity half step, to be used by calculateAngularVelocities
 }
}

template <> void calculateVelocities<autopas::MoleculeLJ>(autopas::AutoPas<autopas::MoleculeLJ> &autoPasContainer,
                        const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT);

template <> void calculateVelocities<autopas::MulticenteredMoleculeLJ>(autopas::AutoPas<autopas::MulticenteredMoleculeLJ> &autoPasContainer,
                                              const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT);


template <class ParticleClass>
void calculateAngularVelocities(autopas::AutoPas<ParticleClass> &autoPasContainer,
                               const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
 autopas::utils::ExceptionHandler::exception("calculateAngularVelocities should not be run with a non-rotational molecule type!");
}

template<> void calculateAngularVelocities<autopas::MulticenteredMoleculeLJ>(autopas::AutoPas<autopas::MulticenteredMoleculeLJ> &autoPasContainer,
                                                                 const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
 using autopas::utils::ArrayMath::mulScalar;
 using autopas::utils::ArrayMath::div;
 using autopas::utils::quaternion::rotatePosition;
 using autopas::utils::quaternion::rotatePositionBackwards;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
 for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
   const auto torqueW = iter->getTorque();
   const auto q = iter->getQ();
   const auto I = particlePropertiesLibrary.getMomentOfInertia(iter->getTypeId()); // moment of inertia

   // convert torque to molecular-frame
   const auto torqueM = rotatePositionBackwards(q, torqueW);

   // get I^-1 T in molecular-frame
   const auto torqueDivMoIM = div(torqueM, I);

   // convert to world-frame
   const auto torqueDivMoIW = rotatePosition(q, torqueDivMoIM);

   iter->addAngularVel(mulScalar(torqueDivMoIW, 0.5*deltaT)); // (28)
 }
}

}  // namespace TimeDiscretization