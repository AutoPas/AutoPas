/**
 * @file TimeDiscretization.cpp
 * @author N. Fottner
 * @date 13/05/19
 */
#include "TimeDiscretization.h"
#include "../../src/autopas/utils/Quaternion.h"
#include "autopas/utils/ArrayMath.h"

namespace TimeDiscretization {

#if not defined MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH or MD_FLEXIBLE_MODE!=MULTISITE
void calculatePositionsAndResetForces(autopas::AutoPas<ParticleType> &autoPasContainer,
                                      const ParticlePropertiesLibraryType &particlePropertiesLibrary,
                                      const double &deltaT, const std::array<double, 3> &globalForce,
                                      bool fastParticlesThrow) {
  using autopas::utils::ArrayUtils::operator<<;
  using autopas::utils::ArrayMath::dot;
  using namespace autopas::utils::ArrayMath::literals;

  const auto maxAllowedDistanceMoved = autoPasContainer.getVerletSkinPerTimestep() / 2.;
  const auto maxAllowedDistanceMovedSquared = maxAllowedDistanceMoved * maxAllowedDistanceMoved;

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
    v *= deltaT;
    f *= (deltaT * deltaT / (2 * m));
    const auto displacement = v + f;
    // sanity check that particles are not too fast for the Verlet skin technique.
    // If this condition is violated once this is not necessarily an error. Only if the total distance traveled over
    // the whole rebuild frequency is farther than the skin we lose interactions.
    const auto distanceMovedSquared = dot(displacement, displacement);
    if (distanceMovedSquared > maxAllowedDistanceMovedSquared) {
#pragma omp critical
      std::cerr << "A particle moved farther than verletSkinPerTimestep/2: " << std::sqrt(distanceMovedSquared) << " > "
                << autoPasContainer.getVerletSkinPerTimestep() << "/2 = " << maxAllowedDistanceMoved << "\n"
                << *iter << "\nNew Position: " << iter->getR() + displacement << std::endl;
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
#else

//needed for reduction in this calculatePositionsAndResetForces
#pragma omp declare reduction(vec_of_arrays_plus : std::vector<std::array<double, 3>> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(),[](const auto& lhs, const auto& rhs){return std::array<double,3>{lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]};})) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))


void calculatePositionsAndResetForces(autopas::AutoPas<ParticleType> &autoPasContainer, MoleculeContainer& moleculeContainer,
                                      const ParticlePropertiesLibraryType &particlePropertiesLibrary,
                                      const double &deltaT, const std::array<double, 3> &globalForce,
                                      bool fastParticlesThrow){
  using autopas::utils::ArrayUtils::operator<<;
  using autopas::utils::ArrayMath::dot;
  using namespace autopas::utils::ArrayMath::literals;

  const auto maxAllowedDistanceMoved = autoPasContainer.getVerletSkinPerTimestep() / 2.;
  const auto maxAllowedDistanceMovedSquared = maxAllowedDistanceMoved * maxAllowedDistanceMoved;

  bool throwException = false;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for shared(moleculeContainer) default(none)
#endif
  //store oldForce of all molecules
  for(size_t i=0; i < moleculeContainer.size(); i++) {
    auto molecule = moleculeContainer.get(i);
    molecule.setOldF(molecule.getF());
    molecule.setF({0,0,0});
  }

  //accumulate all forces working on molecules rn (i am sure there is a prettier way to this that doesn't involve all this SoA-, AoS-juggling)
  std::vector<std::array<double, 3>> forceAccumulator = std::vector<std::array<double, 3>>(moleculeContainer.size(), {0., 0., 0.});
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(vec_of_arrays_plus : forceAccumulator) shared(autoPasContainer, globalForce) default(none)
#endif
  //accumulate all forces acting on sites of a molecule in that molecule
  for(auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    //MoleculeType molecule = moleculeContainer.get(iter->getMoleculeId());
    //molecule.addF(iter->getF());    //@TODO considering that acting force shouldn't just be translated as a translational force there should be some factor based on the angle between the relative site position and the force vector
    forceAccumulator[iter->getMoleculeId()] = autopas::utils::ArrayMath::add(forceAccumulator[iter->getMoleculeId()], iter->getF());
    iter->setF(globalForce);
  }

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for shared(forceAccumulator, moleculeContainer) default(none)
  for(size_t i=0; i < moleculeContainer.size(); i++) {
    moleculeContainer.get(i).setF(forceAccumulator[i]);
  }
#endif

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for shared(moleculeContainer, particlePropertiesLibrary, deltaT, maxAllowedDistanceMoved, maxAllowedDistanceMovedSquared, autoPasContainer, fastParticlesThrow, throwException) default(none)
#endif
  //compute velocity and new position based on that info
  for(size_t i=0; i < moleculeContainer.size(); i++) {
    auto molecule = moleculeContainer.get(i);
    auto m = particlePropertiesLibrary.getMolMass(molecule.getTypeId());
    auto v = molecule.getV();
    auto f = molecule.getF();
    v *= deltaT;
    f*= (deltaT * deltaT / (2 * m));
    const auto displacement = v + f;

    // sanity check that particles are not too fast for the Verlet skin technique.
    // If this condition is violated once this is not necessarily an error. Only if the total distance traveled over
    // the whole rebuild frequency is farther than the skin we lose interactions.
    const auto distanceMovedSquared = dot(displacement, displacement);
    if (distanceMovedSquared > maxAllowedDistanceMovedSquared) {
#pragma omp critical
      std::cerr << "A molecule moved farther than verletSkinPerTimestep/2: " << std::sqrt(distanceMovedSquared) << " > "
                << autoPasContainer.getVerletSkinPerTimestep() << "/2 = " << maxAllowedDistanceMoved << "\n"
                << molecule << "\nNew Position: " << molecule.getR() + displacement << std::endl;
      if (fastParticlesThrow) {
        throwException = true;
      }
    }
    molecule.addR(displacement);

    if (throwException) {
      throw std::runtime_error("At least one molecule was too fast!");
    }
  }


  /*
   //old implementation (still kept around because i will definitely need it when debugging)
#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(|| : throwException) shared(autoPasContainer, particlePropertiesLibrary, globalForce, deltaT, fastParticlesThrow, maxAllowedDistanceMoved, maxAllowedDistanceMovedSquared) default(none)
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
    // sanity check that particles are not too fast for the Verlet skin technique.
    // If this condition is violated once this is not necessarily an error. Only if the total distance traveled over
    // the whole rebuild frequency is farther than the skin we lose interactions.
    const auto distanceMovedSquared = dot(displacement, displacement);
    if (distanceMovedSquared > maxAllowedDistanceMovedSquared) {
#pragma omp critical
      std::cerr << "A particle moved farther than verletSkinPerTimestep/2: " << std::sqrt(distanceMovedSquared) << " > "
                << autoPasContainer.getVerletSkinPerTimestep() << "/2 = " << maxAllowedDistanceMoved << "\n"
                << *iter << "\nNew Position: " << iter->getR() + displacement << std::endl;
      if (fastParticlesThrow) {
        throwException = true;
      }
    }
    iter->addR(displacement);
  }

  if (throwException) {
    throw std::runtime_error("At least one particle was too fast!");
  }*/
}
#endif

#if not defined MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH or MD_FLEXIBLE_MODE!=MULTISITE
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
  using autopas::utils::quaternion::getRotationBetweenVectors;

#if MD_FLEXIBLE_MODE == MULTISITE

  const auto halfDeltaT = 0.5 * deltaT;

  const double tol = 1e-13;  // tolerance given in paper
  const double tolSquared = tol * tol;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
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

#if not defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
    iter->setQuaternion(qFullStep);
#else
    iter->setQuaternion(qFullStep, particlePropertiesLibrary);
    //iter->setQuaternion(qFullStep, particlePropertiesLibrary.getSitePositions(iter->getTypeId()));
#endif

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
#else
void calculateQuaternionsAndResetTorques(autopas::AutoPas<ParticleType> &autoPasContainer, MoleculeContainer& moleculeContainer,
                                         const ParticlePropertiesLibraryType &particlePropertiesLibrary,
                                         const double &deltaT, const std::array<double, 3> &globalForce){
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
  using autopas::utils::quaternion::getRotationBetweenVectors;

#if MD_FLEXIBLE_MODE == MULTISITE

  const auto halfDeltaT = 0.5 * deltaT;

  const double tol = 1e-13;  // tolerance given in paper
  const double tolSquared = tol * tol;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for shared(autoPasContainer) default(none)
  for(auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
  }
#endif
  //gather total torque acting on molecules

  //do actual quaternion-calculations

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for
#endif
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

#if not defined(MD_FLEXIBLE_FUNCTOR_ABSOLUTE_POS)
    iter->setQuaternion(qFullStep);
#else
    iter->setQuaternion(qFullStep, particlePropertiesLibrary);
    //iter->setQuaternion(qFullStep, particlePropertiesLibrary.getSitePositions(iter->getTypeId()));
#endif

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
#endif

void calculateVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  // helper declarations for operations with vector
  using namespace autopas::utils::ArrayMath::literals;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto molecularMass = particlePropertiesLibrary.getMolMass(iter->getTypeId());
    const auto force = iter->getF();
    const auto oldForce = iter->getOldF();
    const auto changeInVel = (force + oldForce) * (deltaT / (2 * molecularMass));
    iter->addV(changeInVel);
  }
}

void calculateAngularVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
                                const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  using namespace autopas::utils::ArrayMath::literals;
  using autopas::utils::quaternion::rotatePosition;
  using autopas::utils::quaternion::rotatePositionBackwards;
  using autopas::utils::quaternion::getRotationBetweenVectors;
  using autopas::utils::ArrayMath::normalize;

#if MD_FLEXIBLE_MODE == MULTISITE

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
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
