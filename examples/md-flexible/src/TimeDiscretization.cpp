/**
 * @file TimeDiscretization.cpp
 * @author N. Fottner
 * @date 13/05/19
 */
#include "TimeDiscretization.h"
#include "../../src/autopas/utils/Quaternion.h"
#include "autopas/utils/ArrayMath.h"
#include <stdexcept>

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
#pragma omp parallel reduction(|| : throwException) shared(autoPasContainer, particlePropertiesLibrary, globalForce, deltaT, maxAllowedDistanceMoved, maxAllowedDistanceMovedSquared, fastParticlesThrow, std::cerr) default(none)
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
#pragma omp parallel for reduction(|| : throwException) shared(autoPasContainer, moleculeContainer, particlePropertiesLibrary, globalForce, deltaT, maxAllowedDistanceMovedSquared, maxAllowedDistanceMoved, fastParticlesThrow, std::cerr) default(none)
#endif
  for(size_t i = 0; i < moleculeContainer.size(); i++) {
  auto& molecule = moleculeContainer.get(i);
  const auto m = particlePropertiesLibrary.getMolMass(molecule.getTypeId());
  auto v = molecule.getV();
  auto f = molecule.getF();
  molecule.setOldF(f);
  molecule.setF(globalForce);
  v *=deltaT;
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
              << molecule << "\nNew Position: " << molecule.getR() + displacement << std::endl;
    if (fastParticlesThrow) {
      throwException = true;
    }
  }
  molecule.addR(displacement);
  }

  if (throwException) {
    throw std::runtime_error("At least one particle was too fast!");
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

#if defined MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH and MD_FLEXIBLE_MODE==MULTISITE
void accumulateSiteForcesInMol(autopas::AutoPas<ParticleType> &autoPasContainer, MoleculeContainer& moleculeContainer){
#ifdef AUTOPAS_OPENMP
#pragma omp parallel shared(moleculeContainer, autoPasContainer) default(none)
#endif
  //accumulate all forces acting on sites of a molecule in that molecule
  for(auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    MoleculeType& molecule = moleculeContainer.get(iter->getMoleculeId());
    molecule.addF(iter->getF());    //@TODO considering that acting force shouldn't just be translated as a translational force there should be some factor based on the angle between the relative site position and the force vector
    //iter->setF(globalForce);           //don't reset site force since we need that later for torque calculation
  }
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

  //do actual quaternion-calculations
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for shared(moleculeContainer, particlePropertiesLibrary, halfDeltaT, tol, tolSquared, globalForce, deltaT) default(none)
#endif
  for(size_t i=0; i < moleculeContainer.size(); i++){
    MoleculeType& molecule = moleculeContainer.get(i);
    const auto q = molecule.getQuaternion();
    const auto angVelW = molecule.getAngularVel();// angular velocity in world frame
    const auto angVelM =
        rotatePositionBackwards(q, angVelW);  // angular velocity in molecular frame  (equivalent to (17))
    const auto torqueW = molecule.getTorque();
    const auto torqueM = rotatePositionBackwards(q, torqueW);  // (18)
    const auto I = particlePropertiesLibrary.getMomentOfInertia(molecule.getTypeId());  // moment of inertia

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

   molecule.setQuaternion(qFullStep);

   molecule.setAngularVel(angVelWHalfStep);  // save angular velocity half step, to be used by calculateAngularVelocities

   // Reset torque
   molecule.setTorque({0., 0., 0.});

   if (std::any_of(globalForce.begin(), globalForce.end(),
                   [](double i) { return std::abs(i) > std::numeric_limits<double>::epsilon(); })) {
     // Get torque from global force (this part would be unnecessary if every site simply took the global force)
     const auto unrotatedSitePositions = particlePropertiesLibrary.getSitePositions(molecule.getTypeId());
     const auto rotatedSitePositions = rotateVectorOfPositions(qFullStep, unrotatedSitePositions);
     for (size_t site = 0; site < particlePropertiesLibrary.getNumSites(molecule.getTypeId()); site++) {
        molecule.addTorque(cross(rotatedSitePositions[site], globalForce));
     }
   }
  }


  //update the site position based on the position and rotation of the molecule
  //@TODO: multithread
  for(auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
   MoleculeType& molecule = moleculeContainer.get(iter->getMoleculeId());
   const auto r = molecule.getR();
   const auto q = molecule.getQuaternion();
   const auto rotatedSitePosition = rotatePosition(q,particlePropertiesLibrary.getSitePositions(iter->getTypeId())[iter->getIndexInsideMolecule()]);
   iter->setR(autopas::utils::ArrayMath::add(r, rotatedSitePosition));
  }


#else
  autopas::utils::ExceptionHandler::exception(
      "Attempting to perform rotational integrations when md-flexible has not been compiled with multi-site support!");
#endif
}
#endif

#if defined MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH and MD_FLEXIBLE_MODE==MULTISITE
void gatherTorquesFromForces(autopas::AutoPas<ParticleType> &autoPasContainer, MoleculeContainer& moleculeContainer,
                                         const ParticlePropertiesLibraryType &particlePropertiesLibrary) {
  using autopas::utils::quaternion::rotatePosition;
  //gather total torque acting on molecules
  for(auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    MoleculeType& molecule = moleculeContainer.get(iter->getMoleculeId());
    const auto q = molecule.getQuaternion();
    const auto rotatedSitePosition = rotatePosition(q,particlePropertiesLibrary.getSitePositions(iter->getTypeId())[iter->getIndexInsideMolecule()]);
    const auto torqueOnSite = autopas::utils::ArrayMath::cross(rotatedSitePosition, iter->getF());
    molecule.addTorque(torqueOnSite);
    iter->setF({0.,0.,0.}); //global force instantly gets stored in molecule
  }
}
#endif

#if not defined MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH or MD_FLEXIBLE_MODE!=MULTISITE
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
#else
void calculateVelocities(autopas::AutoPas<ParticleType> &autoPasContainer, MoleculeContainer& moleculeContainer,
                                const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  // helper declarations for operations with vector
  using namespace autopas::utils::ArrayMath::literals;

  //this loop is not thread save for some reason that still boggles my mind
#ifdef AUTOPAS_OPENMP
#pragma omp parallel for shared(moleculeContainer, particlePropertiesLibrary, deltaT) default(none)
#endif
  for(size_t i = 0; i < moleculeContainer.size(); i++) {
   auto& molecule = moleculeContainer.get(i);
   const auto molecularMass = particlePropertiesLibrary.getMolMass(molecule.getTypeId());
   const auto oldForce = molecule.getOldF();
   const auto force = molecule.getF();
   const auto changeInVel = (force + oldForce) * (deltaT / (2 * molecularMass));
   molecule.addV(changeInVel);
  }
}
#endif

#if not defined MD_FLEXIBLE_USE_BUNDLING_MULTISITE_APPROACH or MD_FLEXIBLE_MODE!=MULTISITE
void calculateAngularVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
                                const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  using namespace autopas::utils::ArrayMath::literals;
  using autopas::utils::quaternion::rotatePosition;
  using autopas::utils::quaternion::rotatePositionBackwards;
  using autopas::utils::quaternion::getRotationBetweenVectors;
  using autopas::utils::ArrayMath::normalize;

#if MD_FLEXIBLE_MODE == MULTISITE

  //you cannot parallelize loops using iters that easily with openmp. #pragma omp parallel (without the "for") will only
  // lead to the whole loop being executed by all threads instead of the workload being distributed
  //adding #pragma omp parallel for won't work since openmp doesn't know at the start of the iteration whether there will
  // be a next iteration
//#ifdef AUTOPAS_OPENMP
//#pragma omp parallel
//#endif
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
#else
void calculateAngularVelocities(autopas::AutoPas<ParticleType> &autoPasContainer, MoleculeContainer& moleculeContainer,
                                const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  using namespace autopas::utils::ArrayMath::literals;
  using autopas::utils::quaternion::rotatePosition;
  using autopas::utils::quaternion::rotatePositionBackwards;
  using autopas::utils::quaternion::getRotationBetweenVectors;
  using autopas::utils::ArrayMath::normalize;

#if MD_FLEXIBLE_MODE == MULTISITE

#ifdef AUTOPAS_OPENMP
#pragma omp parallel for shared(moleculeContainer, particlePropertiesLibrary, deltaT) default(none)
#endif
  for(size_t i = 0; i < moleculeContainer.size(); i++) {
    MoleculeType& molecule = moleculeContainer.get(i);
    const auto torqueW = molecule.getTorque();
    const auto q = molecule.getQuaternion();
    const auto I = particlePropertiesLibrary.getMomentOfInertia(molecule.getTypeId());

    // convert torque to molecular-frame
    const auto torqueM = rotatePositionBackwards(q, torqueW);
    // get I^-1 T in molecular-frame
    const auto torqueDivMoIM = torqueM / I;

    // convert to world-frame
    const auto torqueDivMoIW = rotatePosition(q, torqueDivMoIM);

    molecule.addAngularVel(torqueDivMoIW * 0.5 * deltaT);  // (28)
  }

#else
  autopas::utils::ExceptionHandler::exception(
      "Attempting to perform rotational integrations when md-flexible has not been compiled with multi-site support!");
#endif
}
#endif
}  // namespace TimeDiscretization
