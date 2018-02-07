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
    _soa.initArrays({posX, posY, posZ, forceX, forceY, forceZ});
  }

  void moleculesAt(int i, Particle *&rmm_or_not_pointer) override {
    buildMoleculeFromSoA(i, rmm_or_not_pointer);
  }
  void buildMoleculeFromSoA(unsigned int i, Particle *&rmm_or_not_pointer) {
    rmm_or_not_pointer->setR(_soa.readParticle<3>({posX, posY, posZ}, i));
    rmm_or_not_pointer->setF(_soa.readParticle<3>({forceX, forceY, forceZ}, i));
  }
  void addParticle(Particle &m) override {
    _soa.push(posX, m.getR()[0]);
    _soa.push(posY, m.getR()[1]);
    _soa.push(posZ, m.getR()[2]);
    _soa.push(forceX, m.getF()[0]);
    _soa.push(forceY, m.getF()[1]);
    _soa.push(forceZ, m.getF()[2]);
  }

  unsigned long numParticles() const override { return _soa.getNumParticles(); }
  bool isNotEmpty() const override { return numParticles() > 0; }
  SoA _soa;
  enum attributes { posX, posY, posZ, forceX, forceY, forceZ };
};

} /* namespace autopas */

#endif /* AUTOPAS_RMMPARTICLECELL_H_ */
