/*
 * LJFunctor.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_
#define SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_

#include "Functor.h"
#include "particles/MoleculeLJ.h"

namespace autopas {

// can we do this without a template? Maybe. But we want to inline it anyway :)
class LJFunctor: public Functor<MoleculeLJ> {
	void AoSFunctor(MoleculeLJ & i, MoleculeLJ & j) {
		std::array<double, 3> dr, iR, jR;

		double f[3];
		double u;
		double drs[3], dr2; // site distance vector & length^2

		iR = i.getR();
		jR = j.getR();

		double r2 = 0.0;

		for (int d = 0; d < 3 ; ++d ) {
			dr[d] = iR - jR;
			r2 += dr[d] * dr[d];
		}

		if (r2 > CUTOFFSQUARE)
			return;

		double invdr2 = 1. / dr2;
		double lj6 = SIGMA2 * invdr2; lj6 = lj6 * lj6 * lj6;
		double lj12 = lj6 * lj6;
		double lj12m6 = lj12 - lj6;
//		u6 = EPSILON24 * lj12m6;
		double fac = EPSILON24 * (lj12 + lj12m6) * invdr2;
		for (unsigned short d = 0; d < 3; ++d)
			f[d] = fac * dr[d];

		TODO: continue here

	}
	void SoAFunctor() {

	}

	static double CUTOFFSQUARE, EPSILON24, SIGMA2, SHIFT6;

};

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_LJFUNCTOR_H_ */
