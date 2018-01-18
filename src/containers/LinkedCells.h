/*
 * LinkedCells.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_LINKEDCELLS_H_
#define SRC_CONTAINERS_LINKEDCELLS_H_

#include "ParticleContainer.h"

namespace autopas {

template<class Particle, class ParticleCell>
class LinkedCells : public ParticleContainer<Particle, ParticleCell> {
public:
	void init() {
		this->_data.resize(5); // TODO
	}
	void addParticle(Particle& p) {
		this->_data[0].addParticle(p);
	}
	void addParticle(Particle& p, int i) {
		// TODO only for the tests
		this->_data.at(i).addParticle(p); // at performs an out of bounds check
	}
private:
	// ThreeDimensionalCellHandler
};

} /* namespace autopas */

#endif /* SRC_CONTAINERS_LINKEDCELLS_H_ */
