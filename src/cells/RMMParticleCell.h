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

/**
 * Reduced Memory Mode ParticleCell
 * This cell type does not store particles explicitly. Instead, the particles
 * are stored directly in a structure of array.
 * @todo this currently does not store all information of the particles. Only
 * position and forces and thus does not work
 * @tparam Particle type of particle to be stored
 */
template <class Particle>
class RMMParticleCell : public ParticleCell<Particle> {
 public:
  /**
   * Constructor of RMMParticleCell
   */
  RMMParticleCell() {
    _molsSoABuffer.initArrays(
        {Particle::AttributeNames::posX, Particle::AttributeNames::posY,
         Particle::AttributeNames::posZ, Particle::AttributeNames::forceX,
         Particle::AttributeNames::forceY, Particle::AttributeNames::forceZ});
  }

  // TODO: Rethink this function, what is it supposed to do? and what are you
  // supposed to pass?
  void particleAt(int i, Particle *&rmm_or_not_pointer) override {
    buildMoleculeFromSoA(i, rmm_or_not_pointer);
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
    // TODO
  }

  /**
   * the soa buffer of the particle, all information is stored here.
   */
  SoA _molsSoABuffer;

 private:
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
};

} /* namespace autopas */

#endif /* AUTOPAS_RMMPARTICLECELL_H_ */
