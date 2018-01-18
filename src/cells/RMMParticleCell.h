/*
 * RMMParticleCell.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_RMMPARTICLECELL_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_RMMPARTICLECELL_H_

#include "ParticleCell.h"
#include "utils/SoA.h"

namespace autopas {

template<class Particle>
class RMMParticleCell : public ParticleCell<Particle> {
public:
	void moleculesAt(int i, Particle*& rmm_or_not_pointer) {
		buildMoleculeFromSoA(i, rmm_or_not_pointer);
	}
	void buildMoleculeFromSoA(int i, Particle*& rmm_or_not_pointer) {
		rmm_or_not_pointer->setR(_soa.readPos(i));
		rmm_or_not_pointer->setF(_soa.readFor(i));
	}
	void addParticle(Particle& m) {
		_soa.pushPos(m.getR());
		_soa.pushFor(m.getF());
	}
	int numParticles() const {return _soa.getNumParticles();}
	bool isNotEmpty() const { return numParticles() > 0; }
	SoA _soa;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_RMMPARTICLECELL_H_ */
