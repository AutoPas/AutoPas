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

template<class Particle>
class FlopCounterFunctor: public Functor<Particle> {
public:
	FlopCounterFunctor<Particle>(double c) :
			Functor<Particle>(), _cutoffSquare(c * c),
			_distanceCalculations(0ul), _kernelCalls(0ul) {
	}

	void AoSFunctor(Particle & i, Particle & j) override {
		std::array<double, 3> dr = arrayMath::sub(i.getR(), j.getR());
		double dr2 = arrayMath::dot(dr, dr);

		++_distanceCalculations;

		if (dr2 <= _cutoffSquare) ++_kernelCalls;

	}

	double getHitRate() {
		return static_cast<double>(_kernelCalls) / static_cast<double>(_distanceCalculations);
	}

	double getFlops(unsigned numFlopsPerKernelCall) const {
		// 3 sub + 3 square + 2 add
		const double numFlopsPerDistanceCalculation = 8;
		const double distFlops = numFlopsPerDistanceCalculation * static_cast<double>(_distanceCalculations);
		const double kernFlops = numFlopsPerKernelCall * static_cast<double>(_kernelCalls);
		return distFlops + kernFlops;
	}

	unsigned long getDistanceCalculations() const {
		return _distanceCalculations;
	}

	unsigned long getKernelCalls() const {
		return _kernelCalls;
	}

private:
	double _cutoffSquare;
	unsigned long _distanceCalculations, _kernelCalls;
};

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_FLOPCOUNTERFUNCTOR_H_ */
