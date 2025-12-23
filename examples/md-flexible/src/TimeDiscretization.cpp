/**
 * @file TimeDiscretization.cpp
 * @author N. Fottner
 * @date 13/05/19
 */
#include "TimeDiscretization.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <autopas/utils/ExceptionHandler.h>
#include <autopas/utils/WrapOpenMP.h>
#include <autopas/utils/WrapKokkos.h>

#include <autopas/utils/KokkosStorage.h>

namespace TimeDiscretization {

void calculatePositionsAndResetForces(autopas::AutoPas<ParticleType> &autoPasContainer,
                                      const ParticlePropertiesLibraryType &particlePropertiesLibrary,
                                      const double &deltaT, const std::array<double, 3> &globalForce,
                                      bool fastParticlesThrow) {
  using autopas::utils::ArrayUtils::operator<<;
  using autopas::utils::ArrayMath::dot;
  using namespace autopas::utils::ArrayMath::literals;

  bool kokkosForEach = autoPasContainer.containerAllowsKokkos();

  if (!kokkosForEach) {
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    AUTOPAS_OPENMP(parallel)
  #else
    const auto maxAllowedDistanceMoved =
        autoPasContainer.getVerletSkin() / autoPasContainer.getVerletRebuildFrequency() / 2.;
    const auto maxAllowedDistanceMovedSquared = maxAllowedDistanceMoved * maxAllowedDistanceMoved;

    bool throwException = false;

    AUTOPAS_OPENMP(parallel reduction(|| : throwException))
  #endif

    for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      const auto m = particlePropertiesLibrary.getMolMass(iter->getTypeId());
      auto v = iter->getV();
      auto f = iter->getF();
      iter->setOldF(f);
      iter->setF(globalForce);
      v *= deltaT;
      f *= (deltaT * deltaT / (2 * m));
      const auto displacement = v + f;
#ifdef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
      iter->addR(displacement);
#else
      // Sanity check that particles are not too fast for the Verlet skin technique. Only makes sense if skin > 0.
      if (not iter->addRDistanceCheck(displacement, maxAllowedDistanceMovedSquared) and
          maxAllowedDistanceMovedSquared > 0) {
        const auto distanceMoved = std::sqrt(dot(displacement, displacement));
        // If this condition is violated once this is not necessarily an error. Only if the total distance traveled over
        // the whole rebuild frequency is farther than the skin we lose interactions.
        AUTOPAS_OPENMP(critical)
        std::cerr << "A particle moved farther than verletSkinPerTimestep/2: " << distanceMoved << " > "
                  << autoPasContainer.getVerletSkin() / autoPasContainer.getVerletRebuildFrequency()
                  << "/2 = " << maxAllowedDistanceMoved << "\n"
                  << *iter << "\nNew Position: " << iter->getR() + displacement << std::endl;
        if (fastParticlesThrow) {
          throwException = true;
        }
          }
#endif
    }
  } else {
    // TODO: sth else than HostSpace
   autoPasContainer.forEachKokkos<Kokkos::HostSpace::execution_space>(KOKKOS_LAMBDA(int i, const autopas::utils::KokkosStorage<ParticleType>& storage) {
    //auto m = particlePropertiesLibrary.getMolMass(storage.template get<ParticleType::AttributeNames::typeId, true>(i));
    double m = 1.;
    auto vX = storage.template operator()<ParticleType::AttributeNames::velocityX, true>(i);
    auto vY = storage.template operator()<ParticleType::AttributeNames::velocityY, true>(i);
    auto vZ = storage.template operator()<ParticleType::AttributeNames::velocityZ, true>(i);

    auto fX = storage.template operator()<ParticleType::AttributeNames::forceX, true>(i);
    auto fY = storage.template operator()<ParticleType::AttributeNames::forceY, true>(i);
    auto fZ = storage.template operator()<ParticleType::AttributeNames::forceZ, true>(i);

    storage.template operator()<ParticleType::AttributeNames::oldForceX, true>(i) = fX;
    storage.template operator()<ParticleType::AttributeNames::oldForceY, true>(i) = fY;
    storage.template operator()<ParticleType::AttributeNames::oldForceZ, true>(i) = fZ;

    storage.template operator()<ParticleType::AttributeNames::forceX, true>(i) = globalForce[0];
    storage.template operator()<ParticleType::AttributeNames::forceY, true>(i) = globalForce[1];
    storage.template operator()<ParticleType::AttributeNames::forceZ, true>(i) = globalForce[2];

    vX *= deltaT;
    vY *= deltaT;
    vZ *= deltaT;

    fX *= (deltaT * deltaT / (2 * m));
    fY *= (deltaT * deltaT / (2 * m));
    fZ *= (deltaT * deltaT / (2 * m));

    const auto displacementX = vX + fX;
    const auto displacementY = vY + fY;
    const auto displacementZ = vZ + fZ;

    const auto pX = storage.template operator()<ParticleType::AttributeNames::posX, true>(i);
    const auto pY = storage.template operator()<ParticleType::AttributeNames::posY, true>(i);
    const auto pZ = storage.template operator()<ParticleType::AttributeNames::posZ, true>(i);
    storage.template operator()<ParticleType::AttributeNames::posX, true>(i) = pX + displacementX;
    storage.template operator()<ParticleType::AttributeNames::posY, true>(i) = pY + displacementY;
    storage.template operator()<ParticleType::AttributeNames::posZ, true>(i) = pZ + displacementZ;

  }, autopas::IteratorBehavior::owned);
  }
#ifndef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
  if (throwException) {
    throw std::runtime_error("At least one particle was too fast!");
  }
#endif
}

void calculateQuaternionsAndResetTorques(autopas::AutoPas<ParticleType> &autoPasContainer,
                                         const ParticlePropertiesLibraryType &particlePropertiesLibrary,
                                         const double &deltaT, const std::array<double, 3> &globalForce) {
  using namespace autopas::utils::ArrayMath::literals;
  using autopas::utils::ArrayMath::cross;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayMath::L2Norm;
  using autopas::utils::ArrayMath::normalize;
  using autopas::utils::quaternion::qMul;
  using autopas::utils::quaternion::rotatePosition;
  using autopas::utils::quaternion::rotatePositionBackwards;
  using autopas::utils::quaternion::rotateVectorOfPositions;

#if MD_FLEXIBLE_MODE == MULTISITE

  const auto halfDeltaT = 0.5 * deltaT;

  const double tol = 1e-13;  // tolerance given in paper
  const double tolSquared = tol * tol;

  AUTOPAS_OPENMP(parallel)
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    // Calculate Quaternions
    const auto q = iter->getQuaternion();
    const auto angVelW = iter->getAngularVel();  // angular velocity in world frame
    const auto angVelM =
        rotatePositionBackwards(q, angVelW);  // angular velocity in molecular frame  (equivalent to (17))
    const auto torqueW = iter->getTorque();
    const auto torqueM = rotatePositionBackwards(q, torqueW);  // (18)

    const auto I = particlePropertiesLibrary.getMomentOfInertia(iter->getTypeId());  // moment of inertia

    const auto angMomentumM = I * angVelM;                                                 // equivalent to (19)
    const auto derivativeAngMomentumM = torqueM - cross(angVelM, angMomentumM);            // (20)
    const auto angMomentumMHalfStep = angMomentumM + derivativeAngMomentumM * halfDeltaT;  // (21)

    auto derivativeQHalfStep = qMul(q, div(angMomentumMHalfStep, I)) * 0.5;  // (22)

    auto qHalfStep = normalize(q + derivativeQHalfStep * halfDeltaT);  // (23)

    const auto angVelWHalfStep = angVelW + rotatePosition(q, torqueM / I) * halfDeltaT;  // equivalent to (24)

    // (25) start
    // initialise qHalfStepOld to be outside tolerable distance from qHalfStep to satisfy while statement
    auto qHalfStepOld = qHalfStep;
    qHalfStepOld[0] += 2 * tol;

    while (dot(qHalfStep - qHalfStepOld, qHalfStep - qHalfStepOld) > tolSquared) {
      qHalfStepOld = qHalfStep;
      const auto angVelMHalfStep =
          rotatePositionBackwards(qHalfStepOld, angVelWHalfStep);  // equivalent to first two lines of (25)
      derivativeQHalfStep = qMul(qHalfStepOld, angVelMHalfStep) * 0.5;
      qHalfStep = normalize(q + derivativeQHalfStep * halfDeltaT);
    }
    // (25) end

    const auto qFullStep = normalize(q + derivativeQHalfStep * deltaT);  // (26)

    iter->setQuaternion(qFullStep);
    iter->setAngularVel(angVelWHalfStep);  // save angular velocity half step, to be used by calculateAngularVelocities

    // Reset torque
    iter->setTorque({0., 0., 0.});
    if (std::any_of(globalForce.begin(), globalForce.end(),
                    [](double i) { return std::abs(i) > std::numeric_limits<double>::epsilon(); })) {
      // Get torque from global force
      const auto unrotatedSitePositions = particlePropertiesLibrary.getSitePositions(iter->getTypeId());
      const auto rotatedSitePositions = rotateVectorOfPositions(qFullStep, unrotatedSitePositions);
      for (size_t site = 0; site < particlePropertiesLibrary.getNumSites(iter->getTypeId()); site++) {
        iter->addTorque(cross(rotatedSitePositions[site], globalForce));
      }
    }
  }

#else
  autopas::utils::ExceptionHandler::exception(
      "Attempting to perform rotational integrations when md-flexible has not been compiled with multi-site support!");
#endif
}

void calculateVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  // helper declarations for operations with vector
  using namespace autopas::utils::ArrayMath::literals;

  bool kokkosForEach = autoPasContainer.containerAllowsKokkos();

  if (!kokkosForEach) {
    AUTOPAS_OPENMP(parallel)
    for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
      const auto molecularMass = particlePropertiesLibrary.getMolMass(iter->getTypeId());
      const auto force = iter->getF();
      const auto oldForce = iter->getOldF();
      const auto changeInVel = (force + oldForce) * (deltaT / (2 * molecularMass));
      iter->addV(changeInVel);
    }
  }
  else {

    // TODO: sth else than HostSpace
    autoPasContainer.forEachKokkos<Kokkos::HostSpace::execution_space>(KOKKOS_LAMBDA(int i, const autopas::utils::KokkosStorage<ParticleType>& storage) {
      //const auto mass = particlePropertiesLibrary.getMolMass(storage.template get<ParticleType::AttributeNames::typeId, true>(i));
        double mass = 1.;
        auto vX = storage.template operator()<ParticleType::AttributeNames::velocityX, true>(i);
      auto vY = storage.template operator()<ParticleType::AttributeNames::velocityY, true>(i);
      auto vZ = storage.template operator()<ParticleType::AttributeNames::velocityZ, true>(i);

      const auto fX = storage.template operator()<ParticleType::AttributeNames::forceX, true>(i);
      const auto fY = storage.template operator()<ParticleType::AttributeNames::forceY, true>(i);
      const auto fZ = storage.template operator()<ParticleType::AttributeNames::forceZ, true>(i);

      const auto oldFx = storage.template operator()<ParticleType::AttributeNames::oldForceX, true>(i);
      const auto oldFy = storage.template operator()<ParticleType::AttributeNames::oldForceY, true>(i);
      const auto oldFz = storage.template operator()<ParticleType::AttributeNames::oldForceZ, true>(i);

      const auto vUpdateX = (fX + oldFx) * (deltaT / (2 * mass));
      const auto vUpdateY = (fY + oldFy) * (deltaT / (2 * mass));
      const auto vUpdateZ = (fZ + oldFz) * (deltaT / (2 * mass));

      storage.template operator()<ParticleType::AttributeNames::velocityX, true>(i) = vX + vUpdateX;
      storage.template operator()<ParticleType::AttributeNames::velocityY, true>(i) = vY + vUpdateY;
      storage.template operator()<ParticleType::AttributeNames::velocityZ, true>(i) = vZ + vUpdateZ;
    }, autopas::IteratorBehavior::owned);

  }
}

void calculateAngularVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
                                const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  using namespace autopas::utils::ArrayMath::literals;
  using autopas::utils::quaternion::rotatePosition;
  using autopas::utils::quaternion::rotatePositionBackwards;

#if MD_FLEXIBLE_MODE == MULTISITE

  AUTOPAS_OPENMP(parallel)
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto torqueW = iter->getTorque();
    const auto q = iter->getQuaternion();
    const auto I = particlePropertiesLibrary.getMomentOfInertia(iter->getTypeId());  // moment of inertia

    // convert torque to molecular-frame
    const auto torqueM = rotatePositionBackwards(q, torqueW);

    // get I^-1 T in molecular-frame
    const auto torqueDivMoIM = torqueM / I;

    // convert to world-frame
    const auto torqueDivMoIW = rotatePosition(q, torqueDivMoIM);

    iter->addAngularVel(torqueDivMoIW * 0.5 * deltaT);  // (28)
  }

#else
  autopas::utils::ExceptionHandler::exception(
      "Attempting to perform rotational integrations when md-flexible has not been compiled with multi-site support!");
#endif
}

}  // namespace TimeDiscretization
