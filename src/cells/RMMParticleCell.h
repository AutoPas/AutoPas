/*
 * RMMParticleCell.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef AUTOPAS_RMMPARTICLECELL_H_
#define AUTOPAS_RMMPARTICLECELL_H_

#include <utils/SoA.h>
#include "ParticleCell.h"

namespace autopas {

template <class Particle>
class RMMParticleCell : public ParticleCell<Particle> {
 public:
  RMMParticleCell() {
    _molsSoABuffer.initArrays(
        {Particle::AttributeNames::posX, Particle::AttributeNames::posY,
         Particle::AttributeNames::posZ, Particle::AttributeNames::forceX,
         Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ});
  }

  void moleculesAt(int i, Particle *&rmm_or_not_pointer) override {
    buildMoleculeFromSoA(i, rmm_or_not_pointer);
  }
  void buildMoleculeFromSoA(unsigned int i, Particle *&rmm_or_not_pointer) {
    rmm_or_not_pointer->setR(_molsSoABuffer.read<3>(
        {Particle::AttributeNames::posX, Particle::AttributeNames::posY,
         Particle::AttributeNames::posZ},
        i));
    rmm_or_not_pointer->setF(_molsSoABuffer.read<3>(
        {Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
         Particle::AttributeNames::forceZ},
        i));
  }
  void addParticle(Particle &m) override {
    _molsSoABuffer.push(Particle::AttributeNames::posX, m.getR()[0]);
    _molsSoABuffer.push(Particle::AttributeNames::posY, m.getR()[1]);
    _molsSoABuffer.push(Particle::AttributeNames::posZ, m.getR()[2]);
    _molsSoABuffer.push(Particle::AttributeNames::forceX, m.getF()[0]);
    _molsSoABuffer.push(Particle::AttributeNames::forceY, m.getF()[1]);
    _molsSoABuffer.push(Particle::AttributeNames::forceZ, m.getF()[2]);
  }

  unsigned long numParticles() const override {
    return _molsSoABuffer.getNumParticles();
  }
  bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override {
    //TODO
  }
  SoA _molsSoABuffer;
};

} /* namespace autopas */

#endif /* AUTOPAS_RMMPARTICLECELL_H_ */
