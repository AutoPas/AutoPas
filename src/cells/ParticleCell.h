/*
 * ParticleCell.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECELL_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECELL_H_

namespace autopas {

template<class Particle>
class ParticleCell {
public:
	virtual ~ParticleCell() {}
	virtual void addParticle(Particle& p) = 0;
	virtual void moleculesAt(int i, Particle *& rmm_or_not_pointer) = 0; // TODO: consider making return type Particle*
	virtual int numParticles() const = 0;
	virtual bool isNotEmpty() const = 0;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECELL_H_ */
