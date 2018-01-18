/*
 * ParticleContainer.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECONTAINER_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECONTAINER_H_

#include "iterators/ParticleIterator.h"
#include <array>

namespace autopas {

template<class Particle, class ParticleCell>
class ParticleContainer {
public:
	virtual ~ParticleContainer(){}

	typedef ParticleIterator<Particle, ParticleCell> iterator;

	virtual void init() {}

	virtual void addParticle(Particle& p) = 0;

	iterator begin() {return ParticleIterator<Particle, ParticleCell>(&_data);}

	const std::array<double, 3>& getBoxMax() const { return _boxMax; }

	void setBoxMax(const std::array<double, 3>& boxMax) { _boxMax = boxMax; }

	const std::array<double, 3>& getBoxMin() const { return _boxMin; }

	void setBoxMin(const std::array<double, 3>& boxMin) { _boxMin = boxMin; }

	double getCutoff() const { return _cutoff; }

	void setCutoff(double cutoff) { _cutoff = cutoff; }

protected:
	std::vector<ParticleCell> _data;

private:
	std::array<double, 3> _boxMin;
	std::array<double, 3> _boxMax;
	double _cutoff;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECONTAINER_H_ */
