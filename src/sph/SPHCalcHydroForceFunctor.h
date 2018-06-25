/**
 * @file SPHCalcHydroForceFunctor.h
 * @author seckler
 * @date 22.01.18
 */

#pragma once

#include "SPHKernels.h"
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
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType>&, bool)
   */
  void SoAFunctor(SoA<SoAArraysType> &soa, bool newton3) override {
    if (soa.getNumParticles() == 0) return;

    double *const __restrict__ massptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::mass>();
    double *const __restrict__ densityptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::density>();
    double *const __restrict__ smthptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::smth>();
    double *const __restrict__ soundSpeedptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::soundSpeed>();
    double *const __restrict__ pressureptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::pressure>();
    double *const __restrict__ vsigmaxptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::vsigmax>();
    double *const __restrict__ engDotptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::engDot>();

    double *const __restrict__ xptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posZ>();
    double *const __restrict__ velXptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velX>();
    double *const __restrict__ velYptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velY>();
    double *const __restrict__ velZptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velZ>();
    double *const __restrict__ accXptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accX>();
    double *const __restrict__ accYptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accY>();
    double *const __restrict__ accZptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accZ>();

    for (unsigned int i = 0; i < soa.getNumParticles(); ++i) {
      double localvsigmax = 0.;
      double localengdotsum = 0.;
      double localAccX = 0.;
      double localAccY = 0.;
      double localAccZ = 0.;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : localengdotsum, localAccX, localAccY, localAccZ), reduction(max : localvsigmax)
      for (unsigned int j = i + 1; j < soa.getNumParticles(); ++j) {
        const double drx = xptr[i] - xptr[j];
        const double dry = yptr[i] - yptr[j];
        const double drz = zptr[i] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;
        double cutoff = smthptr[i] * autopas::sph::SPHKernels::getKernelSupportRadius();
        if (dr2 >= cutoff * cutoff) continue;

        const double dvX = velXptr[i] - velXptr[j];
        const double dvY = velYptr[i] - velYptr[j];
        const double dvZ = velZptr[i] - velZptr[j];
        // const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;

        double dvdr = dvX * drx + dvY * dry + dvZ * drz;
        const double w_ij = (dvdr < 0) ? dvdr / sqrt(dr2) : 0;
        // const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;

        const double v_sig = soundSpeedptr[i] + soundSpeedptr[j] - 3.0 * w_ij;
        // const PS::F64 v_sig = ep_i[i].snds + ep_j[j].snds - 3.0 * w_ij;

        localvsigmax = std::max(localvsigmax, v_sig);
        if (newton3) {
          vsigmaxptr[j] = std::max(vsigmaxptr[j], v_sig);  // Newton 3
          // v_sig_max = std::max(v_sig_max, v_sig);
        }
        const double AV = -0.5 * v_sig * w_ij / (0.5 * (densityptr[i] + densityptr[j]));
        // const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens +
        // ep_j[j].dens));

        const std::array<double, 3> gradW_ij =
            ArrayMath::mulScalar(ArrayMath::add(SPHKernels::gradW({drx, dry, drz}, smthptr[i]),
                                                SPHKernels::gradW({drx, dry, drz}, smthptr[j])),
                                 0.5);
        // const PS::F64vec gradW_ij = 0.5 * (gradW(dr, ep_i[i].smth) + gradW(dr,
        // ep_j[j].smth));

        double scale =
            pressureptr[i] / (densityptr[i] * densityptr[i]) + pressureptr[j] / (densityptr[j] * densityptr[j]) + AV;
        const double massscale = scale * massptr[j];
        localAccX -= gradW_ij[0] * massscale;
        localAccY -= gradW_ij[1] * massscale;
        localAccZ -= gradW_ij[2] * massscale;
        // hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
        // ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) *
        // gradW_ij;
        if (newton3) {
          const double massscale = scale * massptr[i];
          accXptr[j] += gradW_ij[0] * massscale;
          accYptr[j] += gradW_ij[1] * massscale;
          accZptr[j] += gradW_ij[2] * massscale;
          // Newton3, gradW_ij = -gradW_ji
        }
        double scale2i = massptr[j] * (pressureptr[i] / (densityptr[i] * densityptr[i]) + 0.5 * AV);
        localengdotsum += (gradW_ij[0] * dvX + gradW_ij[1] * dvY + gradW_ij[2] * dvZ) * scale2i;
        // hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
        // ep_i[i].dens) + 0.5 * AV) * dv * gradW_ij;

        if (newton3) {
          double scale2j = massptr[i] * (pressureptr[j] / (densityptr[j] * densityptr[j]) + 0.5 * AV);
          engDotptr[j] += (gradW_ij[0] * dvX + gradW_ij[1] * dvY + gradW_ij[2] * dvZ) * scale2j;
          // Newton 3
        }
      }

      engDotptr[i] += localengdotsum;
      accXptr[i] += localAccX;
      accYptr[i] += localAccY;
      accZptr[i] += localAccZ;
      vsigmaxptr[i] = std::max(localvsigmax, vsigmaxptr[i]);
    }
  }

  /**
   * @copydoc Functor::SoAFunctor(SoA<SoAArraysType>&, SoA<SoAArraysType>&, bool)
   */
  void SoAFunctor(SoA<SoAArraysType> &soa1, SoA<SoAArraysType> &soa2, bool newton3) override {
    if (soa1.getNumParticles() == 0 || soa2.getNumParticles() == 0) return;

    double *const __restrict__ massptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::mass>();
    double *const __restrict__ densityptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::density>();
    double *const __restrict__ smthptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::smth>();
    double *const __restrict__ soundSpeedptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::soundSpeed>();
    double *const __restrict__ pressureptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::pressure>();
    double *const __restrict__ vsigmaxptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::vsigmax>();
    double *const __restrict__ engDotptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::engDot>();

    double *const __restrict__ xptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::posX>();
    double *const __restrict__ yptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::posY>();
    double *const __restrict__ zptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::posZ>();
    double *const __restrict__ velXptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::velX>();
    double *const __restrict__ velYptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::velY>();
    double *const __restrict__ velZptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::velZ>();
    double *const __restrict__ accXptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::accX>();
    double *const __restrict__ accYptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::accY>();
    double *const __restrict__ accZptr1 = soa1.begin<autopas::sph::SPHParticle::AttributeNames::accZ>();

    double *const __restrict__ massptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::mass>();
    double *const __restrict__ densityptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::density>();
    double *const __restrict__ smthptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::smth>();
    double *const __restrict__ soundSpeedptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::soundSpeed>();
    double *const __restrict__ pressureptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::pressure>();
    double *const __restrict__ vsigmaxptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::vsigmax>();
    double *const __restrict__ engDotptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::engDot>();

    double *const __restrict__ xptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::posX>();
    double *const __restrict__ yptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::posY>();
    double *const __restrict__ zptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::posZ>();
    double *const __restrict__ velXptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::velX>();
    double *const __restrict__ velYptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::velY>();
    double *const __restrict__ velZptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::velZ>();
    double *const __restrict__ accXptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::accX>();
    double *const __restrict__ accYptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::accY>();
    double *const __restrict__ accZptr2 = soa2.begin<autopas::sph::SPHParticle::AttributeNames::accZ>();

    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      double localvsigmax = 0.;
      double localengdotsum = 0.;
      double localAccX = 0.;
      double localAccY = 0.;
      double localAccZ = 0.;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : localengdotsum, localAccX, localAccY, localAccZ), reduction(max : localvsigmax)
      for (unsigned int j = 0; j < soa2.getNumParticles(); ++j) {
        const double drx = xptr1[i] - xptr2[j];
        const double dry = yptr1[i] - yptr2[j];
        const double drz = zptr1[i] - zptr2[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;
        double cutoff = smthptr1[i] * autopas::sph::SPHKernels::getKernelSupportRadius();
        if (dr2 >= cutoff * cutoff) continue;

        const double dvX = velXptr1[i] - velXptr2[j];
        const double dvY = velYptr1[i] - velYptr2[j];
        const double dvZ = velZptr1[i] - velZptr2[j];
        // const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;

        double dvdr = dvX * drx + dvY * dry + dvZ * drz;
        const double w_ij = (dvdr < 0) ? dvdr / sqrt(dr2) : 0;
        // const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;

        const double v_sig = soundSpeedptr1[i] + soundSpeedptr2[j] - 3.0 * w_ij;
        // const PS::F64 v_sig = ep_i[i].snds + ep_j[j].snds - 3.0 * w_ij;

        localvsigmax = std::max(localvsigmax, v_sig);
        if (newton3) {
          vsigmaxptr2[j] = std::max(vsigmaxptr2[j], v_sig);  // Newton 3
          // v_sig_max = std::max(v_sig_max, v_sig);
        }
        const double AV = -0.5 * v_sig * w_ij / (0.5 * (densityptr1[i] + densityptr2[j]));
        // const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens +
        // ep_j[j].dens));

        const std::array<double, 3> gradW_ij =
            ArrayMath::mulScalar(ArrayMath::add(SPHKernels::gradW({drx, dry, drz}, smthptr1[i]),
                                                SPHKernels::gradW({drx, dry, drz}, smthptr2[j])),
                                 0.5);
        // const PS::F64vec gradW_ij = 0.5 * (gradW(dr, ep_i[i].smth) + gradW(dr,
        // ep_j[j].smth));

        double scale = pressureptr1[i] / (densityptr1[i] * densityptr1[i]) +
                       pressureptr2[j] / (densityptr2[j] * densityptr2[j]) + AV;
        const double massscale = scale * massptr2[j];
        localAccX -= gradW_ij[0] * massscale;
        localAccY -= gradW_ij[1] * massscale;
        localAccZ -= gradW_ij[2] * massscale;
        // hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
        // ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) *
        // gradW_ij;
        if (newton3) {
          const double massscale = scale * massptr1[i];
          accXptr2[j] += gradW_ij[0] * massscale;
          accYptr2[j] += gradW_ij[1] * massscale;
          accZptr2[j] += gradW_ij[2] * massscale;
          // Newton3, gradW_ij = -gradW_ji
        }
        double scale2i = massptr2[j] * (pressureptr1[i] / (densityptr1[i] * densityptr1[i]) + 0.5 * AV);
        localengdotsum += (gradW_ij[0] * dvX + gradW_ij[1] * dvY + gradW_ij[2] * dvZ) * scale2i;
        // hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
        // ep_i[i].dens) + 0.5 * AV) * dv * gradW_ij;

        if (newton3) {
          double scale2j = massptr1[i] * (pressureptr2[j] / (densityptr2[j] * densityptr2[j]) + 0.5 * AV);
          engDotptr2[j] += (gradW_ij[0] * dvX + gradW_ij[1] * dvY + gradW_ij[2] * dvZ) * scale2j;
          // Newton 3
        }
      }

      engDotptr1[i] += localengdotsum;
      accXptr1[i] += localAccX;
      accYptr1[i] += localAccY;
      accZptr1[i] += localAccZ;
      vsigmaxptr1[i] = std::max(localvsigmax, vsigmaxptr1[i]);
    }
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

    double *const __restrict__ massptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::mass>();
    double *const __restrict__ densityptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::density>();
    double *const __restrict__ smthlngthptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::smth>();
    double *const __restrict__ soundSpeedptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::soundSpeed>();
    double *const __restrict__ pressureptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::pressure>();
    double *const __restrict__ vsigmaxptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::vsigmax>();
    double *const __restrict__ engDotptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::engDot>();
    double *const __restrict__ xptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posX>();
    double *const __restrict__ yptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posY>();
    double *const __restrict__ zptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::posZ>();
    double *const __restrict__ velXptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velX>();
    double *const __restrict__ velYptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velY>();
    double *const __restrict__ velZptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::velZ>();
    double *const __restrict__ accXptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accX>();
    double *const __restrict__ accYptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accY>();
    double *const __restrict__ accZptr = soa.begin<autopas::sph::SPHParticle::AttributeNames::accZ>();

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
      cellIter->setAcceleration({accXPtr[i], accYPtr[i], accZPtr[i]});
    }
  })

  /**
   * Get the number of floating point operations used in one full kernel call
   * @return the number of floating point operations
   */
  static unsigned long getNumFlopsPerKernelCall() {
    ///@todo return correct flopcount
    return 1ul;
  }
};

}  // namespace sph
}  // namespace autopas
