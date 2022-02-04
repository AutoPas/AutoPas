/**
 * @file TimeDiscretization.cpp
 * @author N. Fottner
 * @date 13/05/19
 */
#include "TimeDiscretization.h"

#include "autopas/molecularDynamics/LJFunctor.h"
#include "autopas/molecularDynamics/LJFunctorAVX.h"
#include "autopas/utils/ArrayMath.h"

/**
 * Functions for updating velocities and positions as simulation time progresses.
 */
namespace TimeDiscretization {
void calculatePositions(autopas::AutoPas<ParticleType> &autoPasContainer,
                        const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

  autoPasContainer.forEachParallel(
      [&](auto &p) {
        auto v = p.getV();
        auto m = particlePropertiesLibrary.getMass(p.getTypeId());
        auto f = p.getF();
        p.setOldF(f);
        p.setF({0., 0., 0.});
        v = mulScalar(v, deltaT);
        f = mulScalar(f, (deltaT * deltaT / (2 * m)));
        auto newR = add(v, f);
        p.addR(newR);
      },
      autopas::IteratorBehavior::owned);
}

void calculateVelocities(autopas::AutoPas<ParticleType> &autoPasContainer,
                         const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
  // helper declarations for operations with vector
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::mulScalar;

  autoPasContainer.forEachParallel(
      [&](auto &p) {
        auto m = particlePropertiesLibrary.getMass(p.getTypeId());
        auto force = p.getF();
        auto oldForce = p.getOldF();
        auto newV = mulScalar((add(force, oldForce)), deltaT / (2 * m));
        p.addV(newV);
      },
      autopas::IteratorBehavior::owned);
}
}  // namespace TimeDiscretization
