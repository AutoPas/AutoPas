/**
 * @file SPHCalcHydroForceFunctor.h
 * @author seckler
 * @date 22.01.18
 */

#pragma once

#include "autopas/autopasIncludes.h"
#include "autopas/sph/SPHKernels.h"
#include "autopas/sph/SPHParticle.h"

namespace autopas {
namespace sph {
/**
 * Class that defines the hydrodynamic force functor.
 * It is used to calculate the force based on the given SPH kernels.
 */
class SPHCalcHydroForceFunctor
    : public Functor<SPHParticle, FullParticleCell<SPHParticle>, SPHParticle::SoAArraysType, SPHCalcHydroForceFunctor> {
 public:
  /// particle type
  typedef SPHParticle Particle;
  /// soa arrays type
  typedef SPHParticle::SoAArraysType SoAArraysType;
  /// particle cell type
  typedef FullParticleCell<Particle> ParticleCell;

  SPHCalcHydroForceFunctor()
      // the actual cutoff used is dynamic. 0 is used to pass the sanity check.
      : autopas::Functor<Particle, ParticleCell, SoAArraysType, SPHCalcHydroForceFunctor>(0.){};

  bool isRelevantForTuning() override { return true; }

  bool allowsNewton3() override { return true; }

  bool allowsNonNewton3() override { return true; }

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
    const double w_ij = (dvdr < 0) ? dvdr / ArrayMath::L2Norm(dr) : 0;
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
   * @copydoc Functor::SoAFunctor(SoAView<SoAArraysType>, bool)
   * This functor ignores the newton3 value, as we do not expect any benefit from disabling newton3.
   */
  void SoAFunctor(SoAView<SoAArraysType> soa, bool newton3) override {
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

    for (unsigned int indexFirst = 0; indexFirst < soa.getNumParticles(); ++indexFirst) {
      double localvsigmax = 0.;
      double localengdotsum = 0.;
      double localAccX = 0.;
      double localAccY = 0.;
      double localAccZ = 0.;

      // icpc vectorizes this.
      // g++ only with -ffast-math or -funsafe-math-optimizations
      //#pragma omp simd reduction(+ : localengdotsum, localAccX, localAccY, localAccZ), reduction(max : localvsigmax)
      for (unsigned int j = indexFirst + 1; j < soa.getNumParticles(); ++j) {
        const double drx = xptr[indexFirst] - xptr[j];
        const double dry = yptr[indexFirst] - yptr[j];
        const double drz = zptr[indexFirst] - zptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;
        double cutoff = smthptr[indexFirst] * autopas::sph::SPHKernels::getKernelSupportRadius();
        if (dr2 >= cutoff * cutoff) continue;

        const double dvX = velXptr[indexFirst] - velXptr[j];
        const double dvY = velYptr[indexFirst] - velYptr[j];
        const double dvZ = velZptr[indexFirst] - velZptr[j];
        // const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;

        double dvdr = dvX * drx + dvY * dry + dvZ * drz;
        const double w_ij = (dvdr < 0) ? dvdr / sqrt(dr2) : 0;
        // const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;

        const double v_sig = soundSpeedptr[indexFirst] + soundSpeedptr[j] - 3.0 * w_ij;
        // const PS::F64 v_sig = ep_i[i].snds + ep_j[j].snds - 3.0 * w_ij;

        localvsigmax = std::max(localvsigmax, v_sig);
        // vsigmaxptr[j] = std::max(vsigmaxptr[j], v_sig);  // Newton 3
        vsigmaxptr[j] = vsigmaxptr[j] > v_sig ? vsigmaxptr[j] : v_sig;  // Newton 3
        // v_sig_max = std::max(v_sig_max, v_sig);

        const double AV = -0.5 * v_sig * w_ij / (0.5 * (densityptr[indexFirst] + densityptr[j]));
        // const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens +
        // ep_j[j].dens));

        const std::array<double, 3> gradW_ij =
            ArrayMath::mulScalar(ArrayMath::add(SPHKernels::gradW({drx, dry, drz}, smthptr[indexFirst]),
                                                SPHKernels::gradW({drx, dry, drz}, smthptr[j])),
                                 0.5);
        // const PS::F64vec gradW_ij = 0.5 * (gradW(dr, ep_i[i].smth) + gradW(dr,
        // ep_j[j].smth));

        double scale = pressureptr[indexFirst] / (densityptr[indexFirst] * densityptr[indexFirst]) +
                       pressureptr[j] / (densityptr[j] * densityptr[j]) + AV;
        const double massscale = scale * massptr[j];
        localAccX -= gradW_ij[0] * massscale;
        localAccY -= gradW_ij[1] * massscale;
        localAccZ -= gradW_ij[2] * massscale;
        // hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
        // ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) *
        // gradW_ij;

        const double massscale2 = scale * massptr[indexFirst];
        accXptr[j] += gradW_ij[0] * massscale2;
        accYptr[j] += gradW_ij[1] * massscale2;
        accZptr[j] += gradW_ij[2] * massscale2;
        // Newton3, gradW_ij = -gradW_ji

        double scale2i =
            massptr[j] * (pressureptr[indexFirst] / (densityptr[indexFirst] * densityptr[indexFirst]) + 0.5 * AV);
        localengdotsum += (gradW_ij[0] * dvX + gradW_ij[1] * dvY + gradW_ij[2] * dvZ) * scale2i;
        // hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
        // ep_i[i].dens) + 0.5 * AV) * dv * gradW_ij;

        double scale2j = massptr[indexFirst] * (pressureptr[j] / (densityptr[j] * densityptr[j]) + 0.5 * AV);
        engDotptr[j] += (gradW_ij[0] * dvX + gradW_ij[1] * dvY + gradW_ij[2] * dvZ) * scale2j;
        // Newton 3
      }

      engDotptr[indexFirst] += localengdotsum;
      accXptr[indexFirst] += localAccX;
      accYptr[indexFirst] += localAccY;
      accZptr[indexFirst] += localAccZ;
      vsigmaxptr[indexFirst] = std::max(localvsigmax, vsigmaxptr[indexFirst]);
    }
  }

  /**
   * @copydoc Functor::SoAFunctor(SoAView<SoAArraysType>, SoAView<SoAArraysType>, bool)
   */
  void SoAFunctor(SoAView<SoAArraysType> soa1, SoAView<SoAArraysType> soa2, bool newton3) override {
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

    for (unsigned int indexFirst = 0; indexFirst < soa1.getNumParticles(); ++indexFirst) {
      double localvsigmax = 0.;
      double localengdotsum = 0.;
      double localAccX = 0.;
      double localAccY = 0.;
      double localAccZ = 0.;

      // icpc vectorizes this.
      // g++ only with -ffast-math or -funsafe-math-optimizations
      //#pragma omp simd reduction(+ : localengdotsum, localAccX, localAccY, localAccZ), reduction(max : localvsigmax)
      for (unsigned int j = 0; j < soa2.getNumParticles(); ++j) {
        const double drx = xptr1[indexFirst] - xptr2[j];
        const double dry = yptr1[indexFirst] - yptr2[j];
        const double drz = zptr1[indexFirst] - zptr2[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;
        double cutoff = smthptr1[indexFirst] * autopas::sph::SPHKernels::getKernelSupportRadius();
        if (dr2 >= cutoff * cutoff) continue;

        const double dvX = velXptr1[indexFirst] - velXptr2[j];
        const double dvY = velYptr1[indexFirst] - velYptr2[j];
        const double dvZ = velZptr1[indexFirst] - velZptr2[j];
        // const PS::F64vec dv = ep_i[i].vel - ep_j[j].vel;

        double dvdr = dvX * drx + dvY * dry + dvZ * drz;
        const double w_ij = (dvdr < 0) ? dvdr / sqrt(dr2) : 0;
        // const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;

        const double v_sig = soundSpeedptr1[indexFirst] + soundSpeedptr2[j] - 3.0 * w_ij;
        // const PS::F64 v_sig = ep_i[i].snds + ep_j[j].snds - 3.0 * w_ij;

        localvsigmax = std::max(localvsigmax, v_sig);
        if (newton3) {
          // vsigmaxptr2[j] = std::max(vsigmaxptr2[j], v_sig);  // Newton 3
          vsigmaxptr2[j] = vsigmaxptr2[j] > v_sig ? vsigmaxptr2[j] : v_sig;  // Newton 3
          // v_sig_max = std::max(v_sig_max, v_sig);
        }
        const double AV = -0.5 * v_sig * w_ij / (0.5 * (densityptr1[indexFirst] + densityptr2[j]));
        // const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens +
        // ep_j[j].dens));

        const std::array<double, 3> gradW_ij =
            ArrayMath::mulScalar(ArrayMath::add(SPHKernels::gradW({drx, dry, drz}, smthptr1[indexFirst]),
                                                SPHKernels::gradW({drx, dry, drz}, smthptr2[j])),
                                 0.5);
        // const PS::F64vec gradW_ij = 0.5 * (gradW(dr, ep_i[i].smth) + gradW(dr,
        // ep_j[j].smth));

        double scale = pressureptr1[indexFirst] / (densityptr1[indexFirst] * densityptr1[indexFirst]) +
                       pressureptr2[j] / (densityptr2[j] * densityptr2[j]) + AV;
        const double massscale = scale * massptr2[j];
        localAccX -= gradW_ij[0] * massscale;
        localAccY -= gradW_ij[1] * massscale;
        localAccZ -= gradW_ij[2] * massscale;
        // hydro[i].acc     -= ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
        // ep_i[i].dens) + ep_j[j].pres / (ep_j[j].dens * ep_j[j].dens) + AV) *
        // gradW_ij;
        if (newton3) {
          const double massscale = scale * massptr1[indexFirst];
          accXptr2[j] += gradW_ij[0] * massscale;
          accYptr2[j] += gradW_ij[1] * massscale;
          accZptr2[j] += gradW_ij[2] * massscale;
          // Newton3, gradW_ij = -gradW_ji
        }
        double scale2i =
            massptr2[j] * (pressureptr1[indexFirst] / (densityptr1[indexFirst] * densityptr1[indexFirst]) + 0.5 * AV);
        localengdotsum += (gradW_ij[0] * dvX + gradW_ij[1] * dvY + gradW_ij[2] * dvZ) * scale2i;
        // hydro[i].eng_dot += ep_j[j].mass * (ep_i[i].pres / (ep_i[i].dens *
        // ep_i[i].dens) + 0.5 * AV) * dv * gradW_ij;

        if (newton3) {
          double scale2j = massptr1[indexFirst] * (pressureptr2[j] / (densityptr2[j] * densityptr2[j]) + 0.5 * AV);
          engDotptr2[j] += (gradW_ij[0] * dvX + gradW_ij[1] * dvY + gradW_ij[2] * dvZ) * scale2j;
          // Newton 3
        }
      }

      engDotptr1[indexFirst] += localengdotsum;
      accXptr1[indexFirst] += localAccX;
      accYptr1[indexFirst] += localAccY;
      accZptr1[indexFirst] += localAccZ;
      vsigmaxptr1[indexFirst] = std::max(localvsigmax, vsigmaxptr1[indexFirst]);
    }
  }
  // clang-format off
  /**
   * @copydoc Functor::SoAFunctor(SoAView<SoAArraysType> soa, const size_t indexFirst, const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3)
   */
  // clang-format on
  void SoAFunctor(SoAView<SoAArraysType> soa, const size_t indexFirst,
                  const std::vector<size_t, autopas::AlignedAllocator<size_t>> &neighborList, bool newton3) override {
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

    double localvsigmax = 0.;
    double localengdotsum = 0.;
    double localAccX = 0.;
    double localAccY = 0.;
    double localAccZ = 0.;

    auto &currentList = neighborList;
    size_t listSize = currentList.size();

    // icpc vectorizes this.
    // g++ only with -ffast-math or -funsafe-math-optimizations
    //#pragma omp simd reduction(+ : localengdotsum, localAccX, localAccY, localAccZ), reduction(max : localvsigmax)
    for (unsigned int j = 0; j < listSize; ++j) {
      const double drx = xptr[indexFirst] - xptr[currentList[j]];
      const double dry = yptr[indexFirst] - yptr[currentList[j]];
      const double drz = zptr[indexFirst] - zptr[currentList[j]];

      const double drx2 = drx * drx;
      const double dry2 = dry * dry;
      const double drz2 = drz * drz;

      const double dr2 = drx2 + dry2 + drz2;
      double cutoff = smthptr[indexFirst] * autopas::sph::SPHKernels::getKernelSupportRadius();
      if (dr2 >= cutoff * cutoff) continue;

      const double dvX = velXptr[indexFirst] - velXptr[currentList[j]];
      const double dvY = velYptr[indexFirst] - velYptr[currentList[j]];
      const double dvZ = velZptr[indexFirst] - velZptr[currentList[j]];
      // const PS::F64vec dv = ep_i[i].vel - ep_j[currentList[j]].vel;

      double dvdr = dvX * drx + dvY * dry + dvZ * drz;
      const double w_ij = (dvdr < 0) ? dvdr / sqrt(dr2) : 0;
      // const PS::F64 w_ij = (dv * dr < 0) ? dv * dr / sqrt(dr * dr) : 0;

      const double v_sig = soundSpeedptr[indexFirst] + soundSpeedptr[currentList[j]] - 3.0 * w_ij;
      // const PS::F64 v_sig = ep_i[i].snds + ep_j[currentList[j]].snds - 3.0 * w_ij;

      localvsigmax = std::max(localvsigmax, v_sig);
      if (newton3) {
        // vsigmaxptr[currentList[j]] = std::max(vsigmaxptr[currentList[j]], v_sig);  // Newton 3
        vsigmaxptr[currentList[j]] =
            vsigmaxptr[currentList[j]] > v_sig ? vsigmaxptr[currentList[j]] : v_sig;  // Newton 3
        // v_sig_max = std::max(v_sig_max, v_sig);
      }
      const double AV = -0.5 * v_sig * w_ij / (0.5 * (densityptr[indexFirst] + densityptr[currentList[j]]));
      // const PS::F64 AV = - 0.5 * v_sig * w_ij / (0.5 * (ep_i[i].dens +
      // ep_j[currentList[j]].dens));

      const std::array<double, 3> gradW_ij =
          ArrayMath::mulScalar(ArrayMath::add(SPHKernels::gradW({drx, dry, drz}, smthptr[indexFirst]),
                                              SPHKernels::gradW({drx, dry, drz}, smthptr[currentList[j]])),
                               0.5);
      // const PS::F64vec gradW_ij = 0.5 * (gradW(dr, ep_i[i].smth) + gradW(dr,
      // ep_j[currentList[j]].smth));

      double scale = pressureptr[indexFirst] / (densityptr[indexFirst] * densityptr[indexFirst]) +
                     pressureptr[currentList[j]] / (densityptr[currentList[j]] * densityptr[currentList[j]]) + AV;
      const double massscale = scale * massptr[currentList[j]];
      localAccX -= gradW_ij[0] * massscale;
      localAccY -= gradW_ij[1] * massscale;
      localAccZ -= gradW_ij[2] * massscale;
      // hydro[i].acc     -= ep_j[currentList[j]].mass * (ep_i[i].pres / (ep_i[i].dens *
      // ep_i[i].dens) + ep_j[currentList[j]].pres / (ep_j[currentList[j]].dens * ep_j[currentList[j]].dens) + AV) *
      // gradW_ij;
      if (newton3) {
        const double massscale = scale * massptr[indexFirst];
        accXptr[currentList[j]] += gradW_ij[0] * massscale;
        accYptr[currentList[j]] += gradW_ij[1] * massscale;
        accZptr[currentList[j]] += gradW_ij[2] * massscale;
        // Newton3, gradW_ij = -gradW_ji
      }
      double scale2i = massptr[currentList[j]] *
                       (pressureptr[indexFirst] / (densityptr[indexFirst] * densityptr[indexFirst]) + 0.5 * AV);
      localengdotsum += (gradW_ij[0] * dvX + gradW_ij[1] * dvY + gradW_ij[2] * dvZ) * scale2i;
      // hydro[i].eng_dot += ep_j[currentList[j]].mass * (ep_i[i].pres / (ep_i[i].dens *
      // ep_i[i].dens) + 0.5 * AV) * dv * gradW_ij;

      if (newton3) {
        double scale2j =
            massptr[indexFirst] *
            (pressureptr[currentList[j]] / (densityptr[currentList[j]] * densityptr[currentList[j]]) + 0.5 * AV);
        engDotptr[currentList[j]] += (gradW_ij[0] * dvX + gradW_ij[1] * dvY + gradW_ij[2] * dvZ) * scale2j;
        // Newton 3
      }
    }

    engDotptr[indexFirst] += localengdotsum;
    accXptr[indexFirst] += localAccX;
    accYptr[indexFirst] += localAccY;
    accZptr[indexFirst] += localAccZ;
    vsigmaxptr[indexFirst] = std::max(localvsigmax, vsigmaxptr[indexFirst]);
  }

  /**
   * @copydoc Functor::getNeededAttr()
   */
  constexpr static const std::array<typename SPHParticle::AttributeNames, 16> getNeededAttr() {
    ///@todo distinguish between N3 and notN3
    return std::array<typename SPHParticle::AttributeNames, 16>{
        SPHParticle::AttributeNames::mass,     SPHParticle::AttributeNames::density,
        SPHParticle::AttributeNames::smth,     SPHParticle::AttributeNames::soundSpeed,
        SPHParticle::AttributeNames::pressure, SPHParticle::AttributeNames::vsigmax,
        SPHParticle::AttributeNames::engDot,   SPHParticle::AttributeNames::posX,
        SPHParticle::AttributeNames::posY,     SPHParticle::AttributeNames::posZ,
        SPHParticle::AttributeNames::velX,     SPHParticle::AttributeNames::velY,
        SPHParticle::AttributeNames::velZ,     SPHParticle::AttributeNames::accX,
        SPHParticle::AttributeNames::accY,     SPHParticle::AttributeNames::accZ};
  }

  /**
   * @copydoc Functor::getNeededAttr(std::false_type)
   */
  constexpr static const std::array<typename SPHParticle::AttributeNames, 11> getNeededAttr(std::false_type) {
    ///@todo distinguish between N3 and notN3
    return std::array<typename SPHParticle::AttributeNames, 11>{
        SPHParticle::AttributeNames::mass,     SPHParticle::AttributeNames::density,
        SPHParticle::AttributeNames::smth,     SPHParticle::AttributeNames::soundSpeed,
        SPHParticle::AttributeNames::pressure, SPHParticle::AttributeNames::posX,
        SPHParticle::AttributeNames::posY,     SPHParticle::AttributeNames::posZ,
        SPHParticle::AttributeNames::velX,     SPHParticle::AttributeNames::velY,
        SPHParticle::AttributeNames::velZ};
  }

  /**
   * @copydoc Functor::getComputedAttr()
   */
  constexpr static const std::array<typename sph::SPHParticle::AttributeNames, 5> getComputedAttr() {
    return std::array<typename SPHParticle::AttributeNames, 5>{
        SPHParticle::AttributeNames::vsigmax, SPHParticle::AttributeNames::engDot, SPHParticle::AttributeNames::accX,
        SPHParticle::AttributeNames::accY, SPHParticle::AttributeNames::accZ};
  }

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
