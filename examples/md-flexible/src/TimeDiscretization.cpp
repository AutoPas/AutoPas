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
/**
 * Calculate and update the position for every particle using the Störmer-Verlet Algorithm.
 * @param autopas
 * @param particlePropertiesLibrary
 * @param deltaT time step width
 */
void calculatePositions(autopas::AutoPas<ParticleType> &autoPasContainer,
                        const ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT) {
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
    iter->setF({0., 0., 0.});
    v = mulScalar(v, deltaT);
    f = mulScalar(f, (deltaT * deltaT / (2 * m)));
    auto newR = add(v, f);
    iter->addR(newR);
  }
}

/**
 * Calculate and update the velocity for every particle using the the Störmer-Verlet Algorithm.
 * @param autopas
 * @param particlePropertiesLibrary
 * @param deltaT time step width
 */
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

void calculatePairwiseForces(autopas::AutoPas<ParticleType> &autoPasContainer,
                             ParticlePropertiesLibraryType &particlePropertiesLibrary, const double &deltaT,
                             MDFlexConfig::FunctorOption functorOption, bool &wasTuningIteration) {
  switch (functorOption) {
    case MDFlexConfig::FunctorOption::lj12_6: {
      autopas::LJFunctor<ParticleType, true, true> functor{autoPasContainer.getCutoff(), particlePropertiesLibrary};
      wasTuningIteration = autoPasContainer.iteratePairwise(&functor);
      break;
    }
    case MDFlexConfig::FunctorOption::lj12_6_Globals: {
      autopas::LJFunctor<ParticleType, true, true, autopas::FunctorN3Modes::Both, true> functor{
          autoPasContainer.getCutoff(), particlePropertiesLibrary};
      wasTuningIteration = autoPasContainer.iteratePairwise(&functor);
      break;
    }
    case MDFlexConfig::FunctorOption::lj12_6_AVX: {
      autopas::LJFunctorAVX<ParticleType, true, true> functor{autoPasContainer.getCutoff(), particlePropertiesLibrary};
      wasTuningIteration = autoPasContainer.iteratePairwise(&functor);
      break;
    }
  }
}

void calculateGlobalForces(autopas::AutoPas<ParticleType> &autoPasContainer, const std::array<double, 3> &globalForce) {
#ifdef AUTOPAS_OPENMP
#pragma omp parallel default(none) shared(_autoPasContainer)
#endif
  for (auto particle = autoPasContainer.begin(autopas::IteratorBehavior::owned); particle.isValid(); ++particle) {
    particle->addF(globalForce);
  }
}
}  // namespace TimeDiscretization
