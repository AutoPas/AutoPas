/**
 * @file TimeDiscretization.cpp
 * @author N. Fottner
 * @date 13/05/19
 */
#include "TimeDiscretization.h"

#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/Quaternion.h"

/**
 * Functions for updating velocities and positions as simulation time progresses.
 */
namespace TimeDiscretization {
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

template <class ParticleClass>
void calculateQuaternions(autopas::AutoPas<ParticleClass> &autoPasContainer,
                          const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                          const std::array<double, 3> &globalForce) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::sub;
  using autopas::utils::ArrayMath::mul;
  using autopas::utils::ArrayMath::mulScalar;
  using autopas::utils::ArrayMath::div;
  using autopas::utils::ArrayMath::cross;
  using autopas::utils::quaternion::rotatePosition;

  const auto halfDeltaT = 0.5 * deltaT;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
#endif
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    auto q = iter->getQ();
    auto angVelW = iter->getAngularVel(); // angular velocity in world frame
    auto angVelM = rotatePosition(q,angVelW); // angular velocity in molecular frame  (equivalent to (17))
    auto torqueM = iter->getT();
    auto torqueW = rotatePosition(q,torqueM); // (18)

    auto I = particlePropertiesLibrary.getMomentOfInteria(iter->getTypeId()); // moment of inertia

    auto angMomentumM = mul(I,angVelM); // equivalent to (19)
    auto derivativeAngMomentumM = sub(torqueM, cross(angVelM,angMomentumM)); // (20)
    auto angMomentumMHalfStep = add(angMomentumM, mulScalar(derivativeAngMomentumM, halfDeltaT)); // (21)

    auto derivativeQHalfStep = mulScalar(mul(q, div(angMomentumMHalfStep, I)), 0.5); // (22)

    auto qHalfStep = add(q, mulScalar(derivativeQHalfStep,halfDeltaT)); // (23)

    auto angVelWHalfStep = an
  }
}

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
}  // namespace TimeDiscretization
