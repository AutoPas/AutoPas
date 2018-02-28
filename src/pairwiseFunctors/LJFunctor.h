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
#include "iterators/SingleCellIterator.h"

namespace autopas {

// can we do this without a template? Maybe. But we want to inline it anyway :)
template<class Particle, class ParticleCell>
class LJFunctor : public Functor<Particle, ParticleCell> {
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

  void SoAFunctor(SoA &soa) override {
    double *const __restrict__ x1ptr = soa.begin(Particle::AttributeNames::posX);
    double *const __restrict__ y1ptr = soa.begin(Particle::AttributeNames::posY);
    double *const __restrict__ z1ptr = soa.begin(Particle::AttributeNames::posZ);
    double *const __restrict__ x2ptr = soa.begin(Particle::AttributeNames::posX);
    double *const __restrict__ y2ptr = soa.begin(Particle::AttributeNames::posY);
    double *const __restrict__ z2ptr = soa.begin(Particle::AttributeNames::posZ);

    double *const __restrict__ fx1ptr = soa.begin(Particle::AttributeNames::forceX);
    double *const __restrict__ fy1ptr = soa.begin(Particle::AttributeNames::forceY);
    double *const __restrict__ fz1ptr = soa.begin(Particle::AttributeNames::forceZ);
    double *const __restrict__ fx2ptr = soa.begin(Particle::AttributeNames::forceX);
    double *const __restrict__ fy2ptr = soa.begin(Particle::AttributeNames::forceY);
    double *const __restrict__ fz2ptr = soa.begin(Particle::AttributeNames::forceZ);

    double *const __restrict__ id1ptr = soa.begin(Particle::AttributeNames::id);
    double *const __restrict__ id2ptr = soa.begin(Particle::AttributeNames::id);

    for (unsigned int i = 0; i < soa.getNumParticles(); ++i) {
      double fxacc = 0;
      double fyacc = 0;
      double fzacc = 0;

      // icpc vectorizes this.
      // g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
      for (unsigned int j = i + 1; j < soa.getNumParticles(); ++j) {
        if (id1ptr[i] == id2ptr[j]) {
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
      }

      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
  }

  void SoAFunctor(SoA &soa1, SoA &soa2) override {
    double *const __restrict__ x1ptr = soa1.begin(Particle::AttributeNames::posX);
    double *const __restrict__ y1ptr = soa1.begin(Particle::AttributeNames::posY);
    double *const __restrict__ z1ptr = soa1.begin(Particle::AttributeNames::posZ);
    double *const __restrict__ x2ptr = soa2.begin(Particle::AttributeNames::posX);
    double *const __restrict__ y2ptr = soa2.begin(Particle::AttributeNames::posY);
    double *const __restrict__ z2ptr = soa2.begin(Particle::AttributeNames::posZ);

    double *const __restrict__ fx1ptr = soa1.begin(Particle::AttributeNames::forceX);
    double *const __restrict__ fy1ptr = soa1.begin(Particle::AttributeNames::forceY);
    double *const __restrict__ fz1ptr = soa1.begin(Particle::AttributeNames::forceZ);
    double *const __restrict__ fx2ptr = soa2.begin(Particle::AttributeNames::forceX);
    double *const __restrict__ fy2ptr = soa2.begin(Particle::AttributeNames::forceY);
    double *const __restrict__ fz2ptr = soa2.begin(Particle::AttributeNames::forceZ);

    double *const __restrict__ id1ptr = soa1.begin(Particle::AttributeNames::id);
    double *const __restrict__ id2ptr = soa2.begin(Particle::AttributeNames::id);

    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      double fxacc = 0;
      double fyacc = 0;
      double fzacc = 0;

      // icpc vectorizes this.
      // g++ only with -ffast-math or -funsafe-math-optimizations
#pragma omp simd reduction(+ : fxacc, fyacc, fzacc)
      for (unsigned int j = 0; j < soa2.getNumParticles(); ++j) {
        if (id1ptr[i] == id2ptr[j]) {
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
      }

      fx1ptr[i] += fxacc;
      fy1ptr[i] += fyacc;
      fz1ptr[i] += fzacc;
    }
  }

  void SoALoader(ParticleCell &cell, SoA *soa) override {
    soa->resizeArrays(cell.numParticles());

    double *const __restrict__ idptr = soa->begin(Particle::AttributeNames::id);
    double *const __restrict__ xptr = soa->begin(Particle::AttributeNames::posX);
    double *const __restrict__ yptr = soa->begin(Particle::AttributeNames::posY);
    double *const __restrict__ zptr = soa->begin(Particle::AttributeNames::posZ);
    double *const __restrict__ fxptr = soa->begin(Particle::AttributeNames::forceX);
    double *const __restrict__ fyptr = soa->begin(Particle::AttributeNames::forceY);
    double *const __restrict__ fzptr = soa->begin(Particle::AttributeNames::forceZ);

    SingleCellIterator<Particle, ParticleCell> cellIter(&cell);
    // load particles in SoAs
    for (size_t i = 0; cellIter.isValid(); ++cellIter, ++i) {
      idptr[i] = (*cellIter).getID();
      xptr[i] = (*cellIter).getR()[0];
      yptr[i] = (*cellIter).getR()[1];
      zptr[i] = (*cellIter).getR()[2];
      fxptr[i] = (*cellIter).getF()[0];
      fyptr[i] = (*cellIter).getF()[1];
      fzptr[i] = (*cellIter).getF()[2];
    }
  }

  void SoAExtractor(ParticleCell *cell, SoA *soa) override {

    SingleCellIterator<Particle, ParticleCell> cellIter(cell);

#ifndef NDEBUG
    double *const __restrict__ idptr = soa->begin(Particle::AttributeNames::id);
#endif

    double *const __restrict__ fxptr = soa->begin(Particle::AttributeNames::forceX);
    double *const __restrict__ fyptr = soa->begin(Particle::AttributeNames::forceY);
    double *const __restrict__ fzptr = soa->begin(Particle::AttributeNames::forceZ);

    for (unsigned int i = 0; i < soa->getNumParticles(); ++i, ++cellIter) {
      assert(idptr[i] == (*cellIter).getID());
      (*cellIter).setF({fxptr[i], fyptr[i], fzptr[i]});
    }
  }

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

template<class T, class U>
double LJFunctor<T, U>::CUTOFFSQUARE;

template<class T, class U>
double LJFunctor<T, U>::EPSILON24;

template<class T, class U>
double LJFunctor<T, U>::SIGMASQUARE;

template<class T, class U>
double LJFunctor<T, U>::SHIFT6;

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_ */
