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

	void iteratePairwise(Functor<Particle> f) {

	}

private:
	ParticleCell * getCell() {return &(this->_data.at(0));};
};

} /* namespace autopas */

#endif /* SRC_CONTAINERS_DIRECTSUM_H_ */
