/**
 * @file TimeDiscretization.cpp
 * @author N. Fottner
 * @date 13/05/19
 */
#include "TimeDiscretization.h"

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/ArrayUtils.h"

/**
 * Functions for updating velocities and positions as simulation time progresses.
 */
namespace TimeDiscretization {
void calculatePositions(autopas::AutoPas<ParticleType> &autoPasContainer,
                        const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                        const std::array<double, 3> &globalForce, bool fastParticlesThrow, bool fastParticleWarn) {
  using autopas::utils::ArrayUtils::operator<<;
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayMath::mulScalar;

  const auto maxAllowedDistanceMoved = autoPasContainer.getVerletSkinPerTimestep() / 2.;
  const auto maxAllowedDistanceMovedSquared = maxAllowedDistanceMoved * maxAllowedDistanceMoved;

  bool throwException = false;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel reduction(|| : throwException)
#endif
  for (auto iter = autoPasContainer.begin(autopas::IteratorBehavior::owned); iter.isValid(); ++iter) {
    const auto m = particlePropertiesLibrary.getMass(iter->getTypeId());
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
    if (fastParticleWarn) {
      const auto distanceMovedSquared = dot(displacement, displacement);
      if (distanceMovedSquared > maxAllowedDistanceMovedSquared) {
#pragma omp critical
        std::cerr << "A particle moved farther than verletSkinPerTimestep/2: " << std::sqrt(distanceMovedSquared)
                  << " > " << autoPasContainer.getVerletSkinPerTimestep() << "/2 = " << maxAllowedDistanceMoved << "\n"
                  << *iter << "\nNew Position: " << add(iter->getR(), displacement) << std::endl;
        if (fastParticlesThrow) {
          throwException = true;
        }
      }
    }
    iter->addR(displacement);
  }

  if (throwException) {
    throw std::runtime_error("At least one particle was too fast!");
  }
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
    auto m = particlePropertiesLibrary.getMass(iter->getTypeId());
    auto force = iter->getF();
    auto oldForce = iter->getOldF();
    auto newV = mulScalar((add(force, oldForce)), deltaT / (2 * m));
    iter->addV(newV);
  }
}
}  // namespace TimeDiscretization
