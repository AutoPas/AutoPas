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
	void AoSFunctor(Particle & i, Particle & j) {
		std::array<double, 3> dr = arrayMath::sub(i.getR(), j.getR());
		double dr2 = arrayMath::dot(dr, dr);
		if (dr2 > CUTOFFSQUARE)
			return;

		double invdr2 = 1. / dr2;
		double lj6 = SIGMASQUARE * invdr2;
		lj6 = lj6 * lj6 * lj6;
		double lj12 = lj6 * lj6;
		double lj12m6 = lj12 - lj6;
//		u6 = EPSILON24 * lj12m6 + SHIFT6;
		double fac = EPSILON24 * (lj12 + lj12m6) * invdr2;
		std::array<double, 3> f = arrayMath::mulScalar(dr, fac);
		i.addF(f);
		j.subF(f);
	}
	void SoAFunctor() {

	}

	static void setGlobals(double cutoff, double epsilon, double sigma, double shift) {
		CUTOFFSQUARE = cutoff * cutoff;
		EPSILON24 = epsilon * 24.0;
		SIGMASQUARE = sigma * sigma;
		SHIFT6 = shift * 6.0;
	}

	static double CUTOFFSQUARE, EPSILON24, SIGMASQUARE, SHIFT6;

};


template<class T> double LJFunctor<T>::CUTOFFSQUARE;

template<class T> double LJFunctor<T>::EPSILON24;

template<class T> double LJFunctor<T>::SIGMASQUARE;

template<class T> double LJFunctor<T>::SHIFT6;

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_ */
