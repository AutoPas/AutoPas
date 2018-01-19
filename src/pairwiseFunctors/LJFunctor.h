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
//#include "particles/MoleculeLJ.h"
#include <array>

namespace autopas {

// can we do this without a template? Maybe. But we want to inline it anyway :)
template<class Particle>
class LJFunctor: public Functor<Particle> {
public:
	// todo: add macroscopic quantities
	void AoSFunctor(Particle & i, Particle & j) override {
		std::array<double, 3> dr = arrayMath::sub(i.getR(), j.getR());
		double dr2 = arrayMath::dot(dr, dr);

		if (dr2 > CUTOFFSQUARE)
			return;

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

	void AoSFlopFunctor(Particle & i, Particle & j) override {
		std::array<double, 3> dr = arrayMath::sub(i.getR(), j.getR());
		double dr2 = arrayMath::dot(dr, dr);

		++DISTANCECALCULATIONS;

		if (dr2 <= CUTOFFSQUARE) ++KERNELCALLS;

	}

	void SoAFunctor() override {

	}

	static void setGlobals(double cutoff, double epsilon, double sigma, double shift) {
		CUTOFFSQUARE = cutoff * cutoff;
		EPSILON24 = epsilon * 24.0;
		SIGMASQUARE = sigma * sigma;
		SHIFT6 = shift * 6.0;
	}

	static double getFlops() {
		// 3 sub + 3 square + 2 add
		double distanceFlops = 8;

		// Kernel: 12 = 1 (inverse R squared) + 8 (compute scale) + 3 (apply scale)
		// sum Forces: 6 (forces)
		// kernel total = 12 + 6 = 18
		double kernelFlops = 18;

		return DISTANCECALCULATIONS * distanceFlops + KERNELCALLS * kernelFlops;
	}

	static double CUTOFFSQUARE, EPSILON24, SIGMASQUARE, SHIFT6;
	static unsigned long DISTANCECALCULATIONS, KERNELCALLS;

};


template<class T> double LJFunctor<T>::CUTOFFSQUARE;

template<class T> double LJFunctor<T>::EPSILON24;

template<class T> double LJFunctor<T>::SIGMASQUARE;

template<class T> double LJFunctor<T>::SHIFT6;

template<class T> unsigned long LJFunctor<T>::DISTANCECALCULATIONS = 0ul;

template<class T> unsigned long LJFunctor<T>::KERNELCALLS = 0ul;

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_ */
