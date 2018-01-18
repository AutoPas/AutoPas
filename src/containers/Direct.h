/*
 * Direct.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_DIRECT_H_
#define SRC_CONTAINERS_DIRECT_H_

#include "ParticleContainer.h"

namespace autopas {

template<class Particle, class ParticleCell>
class Direct : public ParticleContainer<Particle, ParticleCell> {
public:
	void init() {
		this->_data.resize(1);
		_singleCell = &(this->_data.front());
	}

	void addParticle(Particle& p) {
		_singleCell->addParticle(p);
	}

private:
	ParticleCell * _singleCell;
};

} /* namespace autopas */

#endif /* SRC_CONTAINERS_DIRECT_H_ */
