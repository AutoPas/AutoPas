/**
 * @file FullParticleCell.h
 *
 * @date 17.01.2018
 * @author tchipevn
 */

#pragma once

#include <vector>
#include "ParticleCell.h"
#include "iterators/SingleCellIterator.h"
#include "utils/SoA.h"

namespace autopas {

/**
 * This class handles the storage of particles in their full form.
 * @tparam Particle
 */
template <class Particle, class SoAArraysType = typename Particle::SoAArraysType>
class FullParticleCell : public ParticleCell<Particle> {
 public:
  FullParticleCell() {}

  void addParticle(Particle &m) override { _particles.push_back(m); }

  virtual SingleCellIteratorWrapper<Particle> begin() override {
    return SingleCellIteratorWrapper<Particle>(new iterator_t(this));
  }

  unsigned long numParticles() const override { return _particles.size(); }

  bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override { _particles.clear(); }

  void deleteByIndex(size_t index) override {
    assert(index < numParticles());

    if (index < numParticles() - 1) {
      std::swap(_particles[index], _particles[numParticles() - 1]);
    }
    _particles.pop_back();
  }

  /**
   * storage of the molecules of the cell
   */
  std::vector<Particle> _particles;

  /**
   * the soa buffer of this cell
   */
  SoA<SoAArraysType> _particleSoABuffer;

  /**
   * friend class iterator
   * @tparam ParticleType
   * @tparam ParticleCellType
   */
  template <class ParticleType, class ParticleCellType>
  friend class SingleCellIterator;

  /**
   * type of the internal iterator
   */
  typedef internal::SingleCellIterator<Particle, FullParticleCell<Particle, SoAArraysType>> iterator_t;
};

}  // namespace autopas
