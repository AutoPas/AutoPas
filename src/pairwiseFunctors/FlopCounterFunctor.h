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

template <class Particle, class ParticleCell>
class FlopCounterFunctor : public Functor<Particle, ParticleCell> {
 public:
  explicit FlopCounterFunctor<Particle, ParticleCell>(double c)
      : autopas::Functor<Particle, ParticleCell>(),
        _cutoffSquare(c * c),
        _distanceCalculations(0ul),
        _kernelCalls(0ul) {}

  void AoSFunctor(Particle &i, Particle &j) override {
    auto dr = arrayMath::sub(i.getR(), j.getR());
    double dr2 = arrayMath::dot(dr, dr);

    ++_distanceCalculations;

    if (dr2 <= _cutoffSquare) ++_kernelCalls;
  }

  enum SoAAttributes { id, posX, posY, posZ };

  void SoAFunctor(SoA &soa) override {
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

  void SoAFunctor(SoA &soa1, SoA &soa2) override {
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

  double getHitRate() {
    return static_cast<double>(_kernelCalls) /
           static_cast<double>(_distanceCalculations);
  }

  double getFlops(unsigned long numFlopsPerKernelCall) const {
    // 3 sub + 3 square + 2 add
    const double numFlopsPerDistanceCalculation = 8;
    const double distFlops = numFlopsPerDistanceCalculation *
                             static_cast<double>(_distanceCalculations);
    const double kernFlops =
        numFlopsPerKernelCall * static_cast<double>(_kernelCalls);
    return distFlops + kernFlops;
  }

  unsigned long getDistanceCalculations() const {
    return _distanceCalculations;
  }

  unsigned long getKernelCalls() const { return _kernelCalls; }

 private:
  double _cutoffSquare;
  unsigned long _distanceCalculations, _kernelCalls;
};

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_FLOPCOUNTERFUNCTOR_H_ */
