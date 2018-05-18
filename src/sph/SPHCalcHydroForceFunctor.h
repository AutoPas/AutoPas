//
// Created by seckler on 22.01.18.
//
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
    : public autopas::Functor<
          SPHParticle, autopas::FullParticleCell<autopas::sph::SPHParticle>> {
 public:
  /**
   * Calculates the contribution of the interaction of particle i and j to the
   * hydrodynamic force.
   * It is not symmetric, because the smoothing lenghts of the two particles can
   * be different.
   * @param i first particle of the interaction
   * @param j second particle of the interaction
   * @param newton3 defines whether or whether not to use newton 3
   */
  void AoSFunctor(SPHParticle &i, SPHParticle &j, bool newton3 = true) override{
    const std::array<double, 3> dr = arrayMath::sub(i.getR(), j.getR());
    // const PS::F64vec dr = ep_i[i].pos - ep_j[j].pos;

    double cutoff = i.getSmoothingLength() *
        autopas::sph::SPHKernels::getKernelSupportRadius();

    if (autopas::arrayMath::dot(dr, dr) >= cutoff * cutoff) {
      return;
    }

    const std::array<double, 3> dv = arrayMath::sub(i.getV(), j.getV());
    // const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;

    double dvdr = arrayMath::dot(dv, dr);
    const double w_ij = (dvdr < 0) ? dvdr / sqrt(arrayMath::dot(dr, dr)) : 0;
    // const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;

    const double v_sig = i.getSoundSpeed() + j.getSoundSpeed() - 3.0 * w_ij;
    // const PS::F64 v_sig = ep_i[i].snds + ep_j[j].snds - 3.0 * w_ij;

    i.checkAndSetVSigMax(v_sig);
    if (newton3) {
      j.checkAndSetVSigMax(v_sig);  // Newton 3
      // v_sig_max = std::max(v_sig_max, v_sig);
    }
    const double AV =
        -0.5 * v_sig * w_ij / (0.5 * (i.getDensity() + j.getDensity()));
    // const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens +
    // ep_j[j].dens));

    const std::array<double, 3> gradW_ij = arrayMath::mulScalar(
        arrayMath::add(SPHKernels::gradW(dr, i.getSmoothingLength()),
                       SPHKernels::gradW(dr, j.getSmoothingLength())),
        0.5);
    // const PS::F64vec gradW_ij = 0.5 * (gradW(dr, ep_i[i].smth) + gradW(dr,
    // ep_j[j].smth));

    double scale = i.getPressure() / (i.getDensity() * i.getDensity()) +
        j.getPressure() / (j.getDensity() * j.getDensity()) + AV;
    i.subAcceleration(arrayMath::mulScalar(gradW_ij, scale * j.getMass()));
    // hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
    // ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) *
    // gradW_ij;
    if (newton3) {
      j.addAcceleration(arrayMath::mulScalar(gradW_ij, scale * i.getMass()));
      // Newton3, gradW_ij = -gradW_ji
    }
    double scale2i =
        j.getMass() *
            (i.getPressure() / (i.getDensity() * i.getDensity()) + 0.5 * AV);
    i.addEngDot(arrayMath::dot(gradW_ij, dv) * scale2i);
    // hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
    // ep_i[i].dens) + 0.5 * AV) * dv * gradW_ij;

    if (newton3) {
      double scale2j =
          i.getMass() *
              (j.getPressure() / (j.getDensity() * j.getDensity()) + 0.5 * AV);
      j.addEngDot(arrayMath::dot(gradW_ij, dv) * scale2j);
      // Newton 3
    }
  }

  /**
   * Get the number of floating point operations used in one full kernel call
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall(){
    ///@todo return correct flopcount
    return 1ul;
  }
};

}  // namespace sph
}  // namespace autopas
