/*
 * FullParticleCell.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_FULLPARTICLECELL_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_FULLPARTICLECELL_H_

#include "ParticleCell.h"
#include <vector>

namespace autopas {

template<class Particle>
class FullParticleCell : public ParticleCell<Particle> {
public:
	void moleculesAt(int i, Particle *& rmm_or_not_pointer) override {
		rmm_or_not_pointer = &_mols.at(i);
	}
	void addParticle(Particle& m) override {
		_mols.push_back(m);
	}
	unsigned long numParticles() const override {return _mols.size();}
	bool isNotEmpty() const override { return numParticles() > 0 ; }
	std::vector<Particle> _mols;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_FULLPARTICLECELL_H_ */
