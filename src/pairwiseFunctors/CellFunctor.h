/*
 * CellFunctor.h
 *
 *  Created on: 22 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_PAIRWISEFUNCTORS_CELLFUNCTOR_H_
#define SRC_PAIRWISEFUNCTORS_CELLFUNCTOR_H_

#include "iterators/SingleCellIterator.h"

namespace autopas {

template<class Particle, class ParticleCell, class ParticleFunctor>
class CellFunctor {
public:
	explicit CellFunctor(ParticleFunctor* f) : _functor(f) {
	}
	void processCell(ParticleCell & cell) {
		processCellAoSN3(cell);
	}
	void processCellPair(ParticleCell & cell1, ParticleCell & cell2) {
		processCellPairAoSN3(cell1, cell2);
	}

	void processCellAoSN3(ParticleCell & cell) {
		SingleCellIterator<Particle, ParticleCell> outer(&cell);
		for (; outer.isValid(); ++outer) {
			Particle & p1 = *outer;

			int ind = outer.getIndex() + 1;

			SingleCellIterator<Particle, ParticleCell> inner(&cell, ind);
			for (; inner.isValid(); ++inner) {
				Particle & p2 = *inner;

				_functor->AoSFunctor(p1, p2);
			}
		}
	}

	void processCellPairAoSN3(ParticleCell& cell1, ParticleCell& cell2) {
		SingleCellIterator<Particle, ParticleCell> outer(&cell1);
		for (; outer.isValid(); ++outer) {
			Particle & p1 = *outer;

			SingleCellIterator<Particle, ParticleCell> inner(&cell2);
			for (; inner.isValid(); ++inner) {
				Particle & p2 = *inner;

				_functor->AoSFunctor(p1, p2);
			}
		}
	}

private:
	ParticleFunctor* _functor;
};

} /* namespace autopas */

#endif /* SRC_PAIRWISEFUNCTORS_CELLFUNCTOR_H_ */
