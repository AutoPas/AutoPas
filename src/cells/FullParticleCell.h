/*
 * FullParticleCell.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef AUTOPAS_SRC_FULLPARTICLECELL_H_
#define AUTOPAS_SRC_FULLPARTICLECELL_H_

#include <vector>
#include <iterators/SingleCellIterator.h>
#include "ParticleCell.h"

namespace autopas {

/**
 * This class handles the storage of particles in their full form.
 * @tparam Particle
 */
template <class Particle>
class FullParticleCell : public ParticleCell<Particle, SingleCellIterator<Particle, FullParticleCell<Particle>>, FullParticleCell<Particle>> {
 public:
  FullParticleCell() {
    _molsSoABuffer.initArrays({
        Particle::AttributeNames::id, Particle::AttributeNames::posX,
        Particle::AttributeNames::posY, Particle::AttributeNames::posZ,
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ,
    });
  }
//  void particleAt(int i, Particle *&rmm_or_not_pointer) override {
//    rmm_or_not_pointer = &_mols.at(i);
//  }
  void addParticle(Particle &m) override { _mols.push_back(m); }
  unsigned long numParticles() const override { return _mols.size(); }
  bool isNotEmpty() const override { return numParticles() > 0; }
  void clear() override { _mols.clear(); }

  /**
   * storage of the molecules of the cell
   */
  std::vector<Particle> _mols;

  /**
   * the soa buffer of this cell
   */
  SoA _molsSoABuffer;

  typedef SingleCellIterator<Particle, FullParticleCell<Particle>> iterator;

  template<class ParticleType, class ParticleCellType>
  friend class SingleCellIterator;
};

} /* namespace autopas */

#endif /* AUTOPAS_SRC_FULLPARTICLECELL_H_ */
