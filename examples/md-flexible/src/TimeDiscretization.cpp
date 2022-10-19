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
                        const std::array<double, 3> &globalForce) {
  using autopas::utils::ArrayUtils::operator<<;
  using autopas::utils::ArrayMath::add;
  using autopas::utils::ArrayMath::dot;
  using autopas::utils::ArrayMath::mulScalar;

  const auto maxAllowedDisplacement =
      autoPasContainer.getVerletSkin() / (2 * autoPasContainer.getVerletRebuildFrequency());
  const auto maxAllowedDisplacementSquared = maxAllowedDisplacement * maxAllowedDisplacement;

#ifdef AUTOPAS_OPENMP
#pragma omp parallel
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
    // If this condition is violated once this is not necessarily an error. Only if the total displacement over
    // the whole rebuild frequency is farther than the skin we lose interactions.
    const auto displacementDistSquared = dot(displacement, displacement);
    if (displacementDistSquared > maxAllowedDisplacementSquared) {
#pragma omp critical
      std::cerr << "A particle moved farther than skin/2/rebuildFrequency: " << std::sqrt(displacementDistSquared)
                << " > " << autoPasContainer.getVerletSkin() << "/2/" << autoPasContainer.getVerletRebuildFrequency()
                << " = " << maxAllowedDisplacement << "\n"
                << *iter << "\nNew Position: " << add(iter->getR(), displacement) << std::endl;
    }
    iter->addR(displacement);
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
