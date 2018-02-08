/*
 * LJFunctor.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_
#define SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_

#include "Functor.h"
#include "utils/arrayMath.h"
#include <array>

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
//#pragma omp simd collapse 2 // TODO: moep?
    for (unsigned int i = 0; i < soa1.getNumParticles(); ++i) {
      for (unsigned int j = 0; j < soa2.getNumParticles(); ++j) {
        if (soa1.read<1>({id}, i)[0] == soa2.read<1>({id}, j)[0]) continue;

        std::array<double, 3> dr =
            arrayMath::sub(soa1.read<3>({posX, posY, posZ}, i),
                           soa2.read<3>({posX, posY, posZ}, j));
        double dr2 = arrayMath::dot(dr, dr);

        if (dr2 > CUTOFFSQUARE) continue;
        double invdr2 = 1. / dr2;
        double lj6 = SIGMASQUARE * invdr2;
        lj6 = lj6 * lj6 * lj6;
        double lj12 = lj6 * lj6;
        double lj12m6 = lj12 - lj6;
        double fac = EPSILON24 * (lj12 + lj12m6) * invdr2;

        std::array<double, 3> f = arrayMath::mulScalar(dr, fac);

        std::array<double, 3> newf1 =
            arrayMath::add(soa1.read<3>({forceX, forceY, forceZ}, i), f);
        std::array<double, 3> newf2 =
            arrayMath::sub(soa2.read<3>({forceX, forceY, forceZ}, j), f);

        soa1.write<3>({forceX, forceY, forceZ}, i, newf1);
        soa2.write<3>({forceX, forceY, forceZ}, j, newf2);
      }
    }
  }

  void SoALoader(std::vector<Particle> &particles, SoA *soa) override {
    soa->initArrays({id, posX, posY, posZ, forceX, forceY, forceZ});

    // load particles in SoAs
    for (auto p : particles) {
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
    // Kernel: 12 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply scale)
    // sum Forces: 6 (forces)
    // kernel total = 12 + 6 = 18
    return 18ul;
  }
};

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
