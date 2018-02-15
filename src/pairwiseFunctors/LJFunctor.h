/*
 * LJFunctor.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_
#define SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_

#include <array>
#include "Functor.h"
#include "utils/arrayMath.h"

namespace autopas {

// can we do this without a template? Maybe. But we want to inline it anyway :)
template <class Particle>
class LJFunctor : public Functor<Particle> {
 public:
  // todo: add macroscopic quantities
  void AoSFunctor(Particle &i, Particle &j) override {
    std::array<double, 3> dr = arrayMath::sub(i.getR(), j.getR());
    double dr2 = arrayMath::dot(dr, dr);

    if (dr2 > CUTOFFSQUARE) return;

    double invdr2 = 1. / dr2;
    double lj6 = SIGMASQUARE * invdr2;
    lj6 = lj6 * lj6 * lj6;
    double lj12 = lj6 * lj6;
    double lj12m6 = lj12 - lj6;
    double fac = EPSILON24 * (lj12 + lj12m6) * invdr2;
    std::array<double, 3> f = arrayMath::mulScalar(dr, fac);
    i.addF(f);
    j.subF(f);
  }

  enum SoAAttributes { id, posX, posY, posZ, forceX, forceY, forceZ };

  void SoAFunctor(SoA &soa1, SoA &soa2) override {
    double *const __restrict__ x1ptr = soa1.begin(posX);
    double *const __restrict__ y1ptr = soa1.begin(posY);
    double *const __restrict__ z1ptr = soa1.begin(posZ);
    double *const __restrict__ x2ptr = soa2.begin(posX);
    double *const __restrict__ y2ptr = soa2.begin(posY);
    double *const __restrict__ z2ptr = soa2.begin(posZ);

    double *const __restrict__ fx1ptr = soa1.begin(forceX);
    double *const __restrict__ fy1ptr = soa1.begin(forceY);
    double *const __restrict__ fz1ptr = soa1.begin(forceZ);
    double *const __restrict__ fx2ptr = soa2.begin(forceX);
    double *const __restrict__ fy2ptr = soa2.begin(forceY);
    double *const __restrict__ fz2ptr = soa2.begin(forceZ);

    double *const __restrict__ id1ptr = soa1.begin(id);
    double *const __restrict__ id2ptr = soa2.begin(id);

    const auto numParticles = soa2.getNumParticles();

    for (unsigned int i = 0; i < numParticles; ++i) {
      double fxacc = 0;
      double fyacc = 0;
      double fzacc = 0;

      // icpc vectorizes this.
      // g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
      for (unsigned int j = 0; j < numParticles; ++j) {
        if (*(id1ptr + i) == *(id2ptr + j)) {
          continue;
        }

        const double drx = x1ptr[i] - x2ptr[j];
        const double dry = y1ptr[i] - y2ptr[j];
        const double drz = z1ptr[i] - z2ptr[j];

        const double drx2 = drx * drx;
        const double dry2 = dry * dry;
        const double drz2 = drz * drz;

        const double dr2 = drx2 + dry2 + drz2;

        if (dr2 > CUTOFFSQUARE) {
          continue;
        }

        const double invdr2 = 1. / dr2;
        const double lj2 = SIGMASQUARE * invdr2;
        const double lj6 = lj2 * lj2 * lj2;
        const double lj12 = lj6 * lj6;
        const double lj12m6 = lj12 - lj6;
        const double fac = EPSILON24 * (lj12 + lj12m6) * invdr2;

        const double fx = drx * fac;
        const double fy = dry * fac;
        const double fz = drz * fac;

        fxacc += fx;
        fyacc += fy;
        fzacc += fz;

        fx2ptr[j] -= fx;
        fy2ptr[j] -= fy;
        fz2ptr[j] -= fz;
        //        }
      }

      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
  }

  void SoALoader(std::vector<Particle> &particles, SoA *soa) override {
    soa->initArrays<7>({id, posX, posY, posZ, forceX, forceY, forceZ});

    // load particles in SoAs
    for (auto &&p : particles) {
      soa->push(id, p.getID());

      soa->push(posX, p.getR()[0]);
      soa->push(posY, p.getR()[1]);
      soa->push(posZ, p.getR()[2]);

      soa->push(forceX, p.getF()[0]);
      soa->push(forceY, p.getF()[1]);
      soa->push(forceZ, p.getF()[2]);
    }
  }

  void SoAExtractor(std::vector<Particle> *particles, SoA *soa) override {
    for (unsigned int i = 0; i < soa->getNumParticles(); ++i) {
      Particle p;
      p.setID((soa->read<1>({id}, i))[0]);
      p.setR(soa->read<3>({posX, posY, posZ}, i));
      p.setF(soa->read<3>({forceX, forceY, forceZ}, i));
      particles->push_back(p);
    }
  };

  static void setGlobals(double cutoff, double epsilon, double sigma,
                         double shift) {
    CUTOFFSQUARE = cutoff * cutoff;
    EPSILON24 = epsilon * 24.0;
    SIGMASQUARE = sigma * sigma;
    SHIFT6 = shift * 6.0;
  }

  static double CUTOFFSQUARE, EPSILON24, SIGMASQUARE, SHIFT6;

  static unsigned long getNumFlopsPerKernelCall() {
    // Kernel: 12 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply
    // scale) sum Forces: 6 (forces) kernel total = 12 + 6 = 18
    return 18ul;
  }
};  // namespace autopas

template <class T>
double LJFunctor<T>::CUTOFFSQUARE;

template <class T>
double LJFunctor<T>::EPSILON24;

template <class T>
double LJFunctor<T>::SIGMASQUARE;

template <class T>
double LJFunctor<T>::SHIFT6;

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_ */
