/**
* @file TimeDiscretization.cpp
* @author N. Fottner & S. Newcome
* @date 13/05/19
*/

#include "TimeDiscretization.h"

namespace TimeDiscretization {

void calculatePositionsAndUpdateForces(autopas::AutoPas<ParticleType> &autoPasContainer,
                       const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                       const std::array<double, 3> &globalForce, bool fastParticlesThrow) {
  using autopas::utils::ArrayUtils::operator<<;
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayMath::mulScalar;

  const auto maxAllowedDistanceMoved =
      autoPasContainer.getVerletSkin() / (2 * autoPasContainer.getVerletRebuildFrequency());

  bool throwException = false;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(|| : throwException)
#endif
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto m = particlePropertiesLibrary.getMolMass(iter->getTypeId());
    auto v = iter->getV();
    auto f = iter->getF();
    iter->setOldF(f);
    iter->setF(globalForce);
    v = mulScalar(v, deltaT);
    f = mulScalar(f, (deltaT * deltaT / (2 * m)));
    const auto displacement = add(v, f);
    // sanity check that particles are not too fast for the Verlet skin technique.
    // If this condition is violated once this is not necessarily an error. Only if the total distance traveled over
    // the whole rebuild frequency is farther than the skin we lose interactions.
    const auto distanceMovedSquared = dot(displacement, displacement);
    if (distanceMovedSquared > maxAllowedDistanceMoved) {
#pragma omp critical
      std::cerr << "A particle moved farther than verletSkinPerTimestep/2: " << std::sqrt(distanceMovedSquared) << " > "
                << autoPasContainer.getVerletSkinPerTimestep() << "/2 = " << maxAllowedDistanceMoved << "\n"
                << *iter << "\nNew Position: " << add(iter->getR(), displacement) << std::endl;
      if (fastParticlesThrow) {
        throwException = true;
      }
    }
    iter->addR(displacement);
  }

  if (throwException) {
    throw std::runtime_error("At least one particle was too fast!");
  }
}

void calculateQuaternions(autopas::AutoPas<ParticleType> &autoPasContainer,
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

#if defined(MD_FLEXIBLE_USE_MULTI_SITE)

 const auto halfDeltaT = 0.5 * deltaT;

 const double tol = 1e-13; // tolerance given in paper

 // todo sort out how to handle global forces.

 std::cout << "here0";
 bool flag = true;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
 for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
   if (flag == true) {
     std::cout << "here1";
   }
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

   if (flag == true) {
     std::cout << "here3";
   }

   int i = 0;
   while (L2Norm(sub(qHalfStep,qHalfStepOld))>tol) {
     qHalfStepOld = qHalfStep;
     auto angVelMHalfStep = rotatePositionBackwards(qHalfStepOld,angVelWHalfStep); // equivalent to first two lines of (25)
     derivativeQHalfStep = mulScalar(qMul(qHalfStepOld,angVelMHalfStep),0.5);
     qHalfStep = normalize(add(q, mulScalar(derivativeQHalfStep, halfDeltaT)));
     if (i > 30) {
       std::cout << i << ":";
     }
   }
   // (25) end

   const auto qFullStep = normalize(add(q, mulScalar(derivativeQHalfStep, deltaT))); // (26)

   if (flag == true) {
     std::cout << "here4";
   }

   iter->setQ(qFullStep);
   iter->setAngularVel(angVelWHalfStep); // save angular velocity half step, to be used by calculateAngularVelocities

   flag = false;
 }

#else
  autopas::utils::ExceptionHandler::exception("Attempting to perform rotational integrations when md-flexible has not been compiled with multi-site support!");
#endif
}


void calculateVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
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


void calculateAngularVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
                                                                 const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
 using autopas::utils::ArrayMath::mulScalar;
 using autopas::utils::ArrayMath::div;
 using autopas::utils::quaternion::rotatePosition;
 using autopas::utils::quaternion::rotatePositionBackwards;

#if defined(MD_FLEXIBLE_USE_MULTI_SITE)

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

#else
 autopas::utils::ExceptionHandler::exception("Attempting to perform rotational integrations when md-flexible has not been compiled with multi-site support!");
#endif
}

}  // namespace TimeDiscretization