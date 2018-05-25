/*
 * FlopCounterFunctor.h
 *
 *  Created on: 22 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_PAIRWISEFUNCTORS_FLOPCOUNTERFUNCTOR_H_
#define SRC_PAIRWISEFUNCTORS_FLOPCOUNTERFUNCTOR_H_

#include "Functor.h"
#include "utils/arrayMath.h"

namespace autopas {
/**
 * This class helps in getting the number of performed floating point
 * operations. It is a functor that only calculated the amount of floating point
 * operations.
 * @todo this class currently is limited to the following case:
 *  - constant cutoff radius
 *  - constant amount of floating point operations for one kernel call (distance
 * < cutoff)
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle, class ParticleCell>
class FlopCounterFunctor : public Functor<Particle, ParticleCell> {
 public:
  /**
   * constructor of FlopCounterFunctor
   * @param c the cutoff radius
   */
  explicit FlopCounterFunctor<Particle, ParticleCell>(double c)
      : autopas::Functor<Particle, ParticleCell>(),
        _cutoffSquare(c * c),
        _distanceCalculations(0ul),
        _kernelCalls(0ul) {}

  void AoSFunctor(Particle &i, Particle &j, bool newton3 = true) override {
    auto dr = arrayMath::sub(i.getR(), j.getR());
    double dr2 = arrayMath::dot(dr, dr);

#pragma omp critical
    {
    ++_distanceCalculations;

    if (dr2 <= _cutoffSquare) ++_kernelCalls;
    };
  }

  void SoAFunctor(SoA &soa, bool newton3 = true) override {
    if (soa.getNumParticles() == 0) return;

    double *const __restrict__ x1ptr =
        soa.begin(Particle::AttributeNames::posX);
    double *const __restrict__ y1ptr =
        soa.begin(Particle::AttributeNames::posY);
    double *const __restrict__ z1ptr =
        soa.begin(Particle::AttributeNames::posZ);

    double *const __restrict__ id1ptr = soa.begin(Particle::AttributeNames::id);

    for (unsigned int i = 0; i < soa.getNumParticles(); ++i) {
      unsigned long distanceCalculationsAcc = 0;
      unsigned long kernelCallsAcc = 0;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : kernelCallsAcc, distanceCalculationsAcc)
      for (unsigned int j = i + 1; j < soa.getNumParticles(); ++j) {
        if (id1ptr[i] == id1ptr[j]) continue;

        ++distanceCalculationsAcc;

        const double drx = x1ptr[i] - x1ptr[j];
        const double dry = y1ptr[i] - y1ptr[j];
        const double drz = z1ptr[i] - z1ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 <= _cutoffSquare) ++kernelCallsAcc;
      }

      _distanceCalculations += distanceCalculationsAcc;
      _kernelCalls += kernelCallsAcc;
    }
  }

  void SoAFunctor(SoA &soa1, SoA &soa2, bool newton3 = true) override {
    double *const __restrict__ x1ptr = soa1.begin(posX);
    double *const __restrict__ y1ptr = soa1.begin(posY);
    double *const __restrict__ z1ptr = soa1.begin(posZ);
    double *const __restrict__ x2ptr = soa2.begin(posX);
    double *const __restrict__ y2ptr = soa2.begin(posY);
    double *const __restrict__ z2ptr = soa2.begin(posZ);

    double *const __restrict__ id1ptr = soa1.begin(id);
    double *const __restrict__ id2ptr = soa2.begin(id);

    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      unsigned long distanceCalculationsAcc = 0;
      unsigned long kernelCallsAcc = 0;

// icpc vectorizes this.
// g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : kernelCallsAcc, distanceCalculationsAcc)
      for (unsigned int j = 0; j < soa2.getNumParticles(); ++j) {
        if (*(id1ptr + i) == *(id2ptr + j)) {
          continue;
        }

        ++distanceCalculationsAcc;

        const double drx = x1ptr[i] - x2ptr[j];
        const double dry = y1ptr[i] - y2ptr[j];
        const double drz = z1ptr[i] - z2ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 <= _cutoffSquare) {
          ++kernelCallsAcc;
        }
      }
      _distanceCalculations += distanceCalculationsAcc;
      _kernelCalls += kernelCallsAcc;
    }
  }

  /**
   * get the hit rate of the pair-wise interaction, i.e. the ratio of the number
   * of kernel calls compared to the number of distance calculations
   * @return the hit rate
   */
  double getHitRate() {
    return static_cast<double>(_kernelCalls) /
           static_cast<double>(_distanceCalculations);
  }

  /**
   * get the total number of flops
   * @param numFlopsPerKernelCall
   * @return
   */
  double getFlops(unsigned long numFlopsPerKernelCall) const {
    // 3 sub + 3 square + 2 add
    const double numFlopsPerDistanceCalculation = 8;
    const double distFlops = numFlopsPerDistanceCalculation *
                             static_cast<double>(_distanceCalculations);
    const double kernFlops =
        numFlopsPerKernelCall * static_cast<double>(_kernelCalls);
    return distFlops + kernFlops;
  }

  /**
   * get the number of calculated distance operations
   * @return
   */
  unsigned long getDistanceCalculations() const {
    return _distanceCalculations;
  }

  /**
   * get the number of kernel calls, i.e. the number of pairs of particles with
   * a distance not larger than the cutoff
   * @return
   */
  unsigned long getKernelCalls() const { return _kernelCalls; }

 private:
  enum SoAAttributes { id, posX, posY, posZ };

  double _cutoffSquare;
  unsigned long _distanceCalculations, _kernelCalls;
};

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_FLOPCOUNTERFUNCTOR_H_ */
