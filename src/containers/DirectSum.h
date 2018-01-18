/*
 * Direct.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_DIRECTSUM_H_
#define SRC_CONTAINERS_DIRECTSUM_H_

#include "ParticleContainer.h"

namespace autopas {

template<class Particle, class ParticleCell>
class DirectSum : public ParticleContainer<Particle, ParticleCell> {
public:
	void init() {
		this->_data.resize(1);
	}

	void addParticle(Particle& p) {
		getCell()->addParticle(p);
	}

	void iteratePairwise(Functor<Particle>* f, bool countFlops = false) {
		// AoS version

		if (not countFlops) {
			for (auto outer = getIt(); outer.isValid(); ++outer) {
				Particle & p1 = *outer;

				int ind = outer.getIndex() + 1;

				for (auto inner = getIt(ind); inner.isValid(); ++inner) {
					Particle & p2 = *inner;

					f->AoSFunctor(p1, p2);
				}
			}
		} else {
			for (auto outer = getIt(); outer.isValid(); ++outer) {
				Particle & p1 = *outer;

				int ind = outer.getIndex() + 1;

				for (auto inner = getIt(ind); inner.isValid(); ++inner) {
					Particle & p2 = *inner;

					f->AoSFlopFunctor(p1, p2);
				}
			}
		}

	}

private:
	// for convenience

	typedef SingleCellIterator<Particle, ParticleCell> singIterator;

	singIterator getIt(int index = 0) {
		return singIterator(getCell(), index);
	}

	ParticleCell * getCell() {return &(this->_data.at(0));};
};

} /* namespace autopas */

#endif /* SRC_CONTAINERS_DIRECTSUM_H_ */
