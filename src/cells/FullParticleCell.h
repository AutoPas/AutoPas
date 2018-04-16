/*
 * FullParticleCell.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef AUTOPAS_SRC_FULLPARTICLECELL_H_
#define AUTOPAS_SRC_FULLPARTICLECELL_H_

#include <iterators/SingleCellIterator.h>
#include <vector>
#include "ParticleCell.h"

namespace autopas {

/**
 * This class handles the storage of particles in their full form.
 * @tparam Particle
 */
template <class Particle>
class FullParticleCell
    : public ParticleCell<
          Particle, SingleCellIterator<Particle, FullParticleCell<Particle>>,
          FullParticleCell<Particle>> {
 public:
  FullParticleCell() {
    _particleSoABuffer.initArrays({
        Particle::AttributeNames::id, Particle::AttributeNames::posX,
        Particle::AttributeNames::posY, Particle::AttributeNames::posZ,
        Particle::AttributeNames::forceX, Particle::AttributeNames::forceY,
        Particle::AttributeNames::forceZ,
    });
  }

  void addParticle(Particle &m) override { _mols.push_back(m); }

  unsigned long numParticles() const override { return _mols.size(); }

  bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override { _mols.clear(); }

  void deleteByIndex(int index) override {
    assert(index >= 0 and index < numParticles());

    if (index < numParticles() - 1) {
      std::swap(_mols[index], _mols[numParticles() - 1]);
    }
    _mols.pop_back();
  }

  /**
   * storage of the molecules of the cell
   */
  std::vector<Particle> _mols;

  /**
   * the soa buffer of this cell
   */
  SoA _particleSoABuffer;

  /**
   * iterator to iterate through ParticleCell
   * If you need to explicitly store this iterator use
   * typename FullParticleCell<ParticleType>::iterator iter;
   */
  typedef SingleCellIterator<Particle, FullParticleCell<Particle>> iterator;

  template <class ParticleType, class ParticleCellType>
  friend class SingleCellIterator;
};

} /* namespace autopas */

#endif /* AUTOPAS_SRC_FULLPARTICLECELL_H_ */
