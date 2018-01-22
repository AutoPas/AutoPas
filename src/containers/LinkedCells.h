/*
 * LinkedCells.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_LINKEDCELLS_H_
#define SRC_CONTAINERS_LINKEDCELLS_H_

#include "ParticleContainer.h"
#include "CellBlock3D.h"

namespace autopas {

template<class Particle, class ParticleCell>
class LinkedCells : public ParticleContainer<Particle, ParticleCell> {
public:
	void init() override {
		this->_data.resize(5); // TODO
	}
	void addParticle(Particle& p) override {
		this->_data[0].addParticle(p);
	}
	void addParticle(Particle& p, int i) {
		// TODO only for the tests
		this->_data.at(i).addParticle(p); // at performs an out of bounds check
	}

	void iteratePairwise(Functor<Particle>* f) override {

	}

private:
	CellBlock3D<ParticleCell> _cellBlock;
	// ThreeDimensionalCellHandler
};

} /* namespace autopas */

#endif /* SRC_CONTAINERS_LINKEDCELLS_H_ */
