/**
 * @file SPHCalcDensityFunctor.h
 * @author seckler
 * @date 19.01.18
 */

#pragma once

#include "SPHKernels.h"
#include "SPHParticle.h"
#include "autopasIncludes.h"

namespace autopas {
namespace sph {
/**
 * Class that defines the density functor.
 * It is used to calculate the density based on the given SPH kernel.
 */
class SPHCalcDensityFunctor : public Functor<SPHParticle, FullParticleCell<SPHParticle>> {
 public:
  typedef SPHParticle Particle;
  typedef SPHParticle::SoAArraysType SoAArraysType;
  typedef FullParticleCell<Particle> ParticleCell;

  /**
   * Calculates the density contribution of the interaction of particle i and j.
   * It is not symmetric, because the smoothing lenghts of the two particles can
   * be different.
   * @param i first particle of the interaction
   * @param j second particle of the interaction
   * @param newton3 defines whether or whether not to use newton 3
   */
  inline void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {
    const std::array<double, 3> dr = ArrayMath::sub(j.getR(), i.getR());  // ep_j[j].pos - ep_i[i].pos;
    const double density =
        j.getMass() * SPHKernels::W(dr, i.getSmoothingLength());  // ep_j[j].mass * W(dr, ep_i[i].smth)
    i.addDensity(density);
    if (newton3) {
      // Newton 3:
      // W is symmetric in dr, so no -dr needed, i.e. we can reuse dr
      const double density2 = i.getMass() * SPHKernels::W(dr, j.getSmoothingLength());
      j.addDensity(density2);
    }
  }

  /**
   * Get the number of floating point operations used in one full kernel call
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall() {
    unsigned long flops = 0;
    flops += 3;                            // calculating dr
    flops += 2 * SPHKernels::getFlopsW();  // flops for calling W
    flops += 2 * 1;                        // calculating density
    flops += 2 * 1;                        // adding density
    return flops;
  }

  /**
   * SoALoader for SPHCalcDensityFunctor.
   * Loads mass, position, smoothing length and density.
   * @param cell
   * @param soa
   * @param offset
   */
  AUTOPAS_FUNCTOR_SOALOADER(
      cell, soa, offset,
      // todo it is probably better to resize the soa only once, before calling
      // SoALoader (verlet-list only)
      soa.resizeArrays(offset + cell.numParticles());

      if (cell.numParticles() == 0) return;

      double *const __restrict__ massptr = soa.template begin<Particle::AttributeNames::mass>();
      double *const __restrict__ xptr = soa.template begin<Particle::AttributeNames::posX>();
      double *const __restrict__ yptr = soa.template begin<Particle::AttributeNames::posY>();
      double *const __restrict__ zptr = soa.template begin<Particle::AttributeNames::posZ>();
      double *const __restrict__ smthptr = soa.template begin<Particle::AttributeNames::smth>();
      double *const __restrict__ densptr = soa.template begin<Particle::AttributeNames::density>();

      auto cellIter = cell.begin();
      // load particles in SoAs
      for (size_t i = offset; cellIter.isValid(); ++cellIter, ++i) {
        massptr[i] = cellIter->getMass();
        xptr[i] = cellIter->getR()[0];
        yptr[i] = cellIter->getR()[1];
        zptr[i] = cellIter->getR()[2];
        smthptr[i] = cellIter->getSmoothingLength();
        densptr[i] = cellIter->getDensity();
      })

  /**
   * SoAExtractor for SPHCalcDensityFunctor.
   * Extracts density.
   * @param cell
   * @param soa
   * @param offset
   */
  AUTOPAS_FUNCTOR_SOAEXTRACTOR(
      cell, soa, offset,
      // function body
      if (cell.numParticles() == 0) return;
      double *const __restrict__ densptr = soa.template begin<Particle::AttributeNames::density>();

      auto cellIter = cell.begin();
      // load particles in SoAs
      for (size_t i = offset; cellIter.isValid(); ++cellIter, ++i) { cellIter->setDensity(densptr[i]); })
};
}  // namespace sph
}  // namespace autopas
