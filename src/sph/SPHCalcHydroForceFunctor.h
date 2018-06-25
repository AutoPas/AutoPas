/**
 * @file SPHCalcHydroForceFunctor.h
 * @author seckler
 * @date 22.01.18
 */

#pragma once

#include "SPHParticle.h"
#include "autopasIncludes.h"

namespace autopas {
namespace sph {
/**
 * Class that defines the hydrodynamic force functor.
 * It is used to calculate the force based on the given SPH kernels.
 */
class SPHCalcHydroForceFunctor
    : public autopas::Functor<SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> {
 public:
  /// particle type
  typedef SPHParticle Particle;
  /// soa arrays type
  typedef SPHParticle::SoAArraysType SoAArraysType;
  /// particle cell type
  typedef FullParticleCell<Particle> ParticleCell;

  /**
   * Calculates the contribution of the interaction of particle i and j to the
   * hydrodynamic force.
   * It is not symmetric, because the smoothing lenghts of the two particles can
   * be different.
   * @param i first particle of the interaction
   * @param j second particle of the interaction
   * @param newton3 defines whether or whether not to use newton 3
   */
  void AoSFunctor(SPHParticle &i, SPHParticle &j, bool newton3 = true) override {
    const std::array<double, 3> dr = ArrayMath::sub(i.getR(), j.getR());
    // const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;

    double cutoff = i.getSmoothingLength() * autopas::sph::SPHKernels::getKernelSupportRadius();

    if (autopas::ArrayMath::dot(dr, dr) >= cutoff * cutoff) {
      return;
    }

    const std::array<double, 3> dv = ArrayMath::sub(i.getV(), j.getV());
    // const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;

    double dvdr = ArrayMath::dot(dv, dr);
    const double w_ij = (dvdr < 0) ? dvdr / sqrt(ArrayMath::dot(dr, dr)) : 0;
    // const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;

    const double v_sig = i.getSoundSpeed() + j.getSoundSpeed() - 3.0 * w_ij;
    // const PS::F64 v_sig = ep_i[i].snds + ep_j[j].snds - 3.0 * w_ij;

    i.checkAndSetVSigMax(v_sig);
    if (newton3) {
      j.checkAndSetVSigMax(v_sig);  // Newton 3
      // v_sig_max = std::max(v_sig_max, v_sig);
    }
    const double AV = -0.5 * v_sig * w_ij / (0.5 * (i.getDensity() + j.getDensity()));
    // const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens +
    // ep_j[j].dens));

    const std::array<double, 3> gradW_ij = ArrayMath::mulScalar(
        ArrayMath::add(SPHKernels::gradW(dr, i.getSmoothingLength()), SPHKernels::gradW(dr, j.getSmoothingLength())),
        0.5);
    // const PS::F64vec gradW_ij = 0.5 * (gradW(dr, ep_i[i].smth) + gradW(dr,
    // ep_j[j].smth));

    double scale =
        i.getPressure() / (i.getDensity() * i.getDensity()) + j.getPressure() / (j.getDensity() * j.getDensity()) + AV;
    i.subAcceleration(ArrayMath::mulScalar(gradW_ij, scale * j.getMass()));
    // hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
    // ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) *
    // gradW_ij;
    if (newton3) {
      j.addAcceleration(ArrayMath::mulScalar(gradW_ij, scale * i.getMass()));
      // Newton3, gradW_ij = -gradW_ji
    }
    double scale2i = j.getMass() * (i.getPressure() / (i.getDensity() * i.getDensity()) + 0.5 * AV);
    i.addEngDot(ArrayMath::dot(gradW_ij, dv) * scale2i);
    // hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
    // ep_i[i].dens) + 0.5 * AV) * dv * gradW_ij;

    if (newton3) {
      double scale2j = i.getMass() * (j.getPressure() / (j.getDensity() * j.getDensity()) + 0.5 * AV);
      j.addEngDot(ArrayMath::dot(gradW_ij, dv) * scale2j);
      // Newton 3
    }
  }

  /**
   * Get the number of floating point operations used in one full kernel call
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall() {
    ///@todo return correct flopcount
    return 1ul;
  }

  /**
   * SoALoader for SPHCalcDensityFunctor.
   * Loads mass, position, smoothing length, density, velocity, speed of sound, pressure, vsigmax, acceleration and
   * engdot.
   * @param cell
   * @param soa
   * @param offset
   */
  AUTOPAS_FUNCTOR_SOALOADER(cell, soa, offset, {
    // todo it is probably better to resize the soa only once, before calling
    // SoALoader (verlet-list only)
    soa.resizeArrays(offset + cell.numParticles());

    if (cell.numParticles() == 0) return;

    auto massptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::mass>();
    auto densityptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::density>();
    auto smthlngthptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::smth>();
    auto soundSpeedptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::soundSpeed>();
    auto pressureptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::pressure>();
    auto vsigmaxptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::vsigmax>();
    auto engDotptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::engDot>();
    auto xptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posX>();
    auto yptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posY>();
    auto zptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posZ>();
    auto velXptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velX>();
    auto velYptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velY>();
    auto velZptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velZ>();
    auto accXptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accX>();
    auto accYptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accY>();
    auto accZptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accZ>();

    auto cellIter = cell.begin();
    // load particles in SoAs
    for (size_t i = offset; cellIter.isValid(); ++cellIter, ++i) {
      massptr[i] = cellIter->getMass();
      densityptr[i] = cellIter->getDensity();
      smthlngthptr[i] = cellIter->getSmoothingLength();
      soundSpeedptr[i] = cellIter->getSoundSpeed();
      pressureptr[i] = cellIter->getPressure();
      vsigmaxptr[i] = cellIter->getVSigMax();
      engDotptr[i] = cellIter->getEngDot();
      xptr[i] = cellIter->getR()[0];
      yptr[i] = cellIter->getR()[1];
      zptr[i] = cellIter->getR()[2];
      velXptr[i] = cellIter->getV()[0];
      velYptr[i] = cellIter->getV()[1];
      velZptr[i] = cellIter->getV()[2];
      accXptr[i] = cellIter->getAcceleration()[0];
      accYptr[i] = cellIter->getAcceleration()[1];
      accZptr[i] = cellIter->getAcceleration()[2];
    }
  })

  /**
   * SoAExtractor for SPHCalcDensityFunctor.
   * Extracts vsigmax, acceleration and engdot.
   * @param cell
   * @param soa
   * @param offset
   */
  AUTOPAS_FUNCTOR_SOAEXTRACTOR(cell, soa, offset, {
    // function body
    if (cell.numParticles() == 0) return;

    double *const __restrict__ vsigmaxPtr = soa.begin<Particle::AttributeNames::vsigmax>();
    double *const __restrict__ engDotPtr = soa.begin<Particle::AttributeNames::engDot>();
    double *const __restrict__ accXPtr = soa.begin<Particle::AttributeNames::accX>();
    double *const __restrict__ accYPtr = soa.begin<Particle::AttributeNames::accY>();
    double *const __restrict__ accZPtr = soa.begin<Particle::AttributeNames::accZ>();

    auto cellIter = cell.begin();
    // load particles in SoAs
    for (size_t i = offset; cellIter.isValid(); ++cellIter, ++i) {
      cellIter->setVSigMax(vsigmaxPtr[i]);
      cellIter->setEngDot(engDotPtr[i]);
      cellIter->setAcceleration({accXPtr[i],accYPtr[i],accZPtr[i]});
    }
  })
};

}  // namespace sph
}  // namespace autopas
