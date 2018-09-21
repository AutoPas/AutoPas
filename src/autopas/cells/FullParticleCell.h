/**
 * @file FullParticleCell.h
 * @date 18.01.2018
 * @author seckler
 */

#pragma once

#include <vector>
#include "autopas/cells/ParticleCell.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class handles the storage of particles in their full form.
 * @tparam Particle
 */
template <class Particle, class SoAArraysType = typename Particle::SoAArraysType>
class FullParticleCell : public ParticleCell<Particle> {
 public:
  FullParticleCell() { autopas_init_lock(&addParticleLock); }

  ~FullParticleCell() { autopas_destroy_lock(&addParticleLock); }

  void addParticle(Particle &m) override {
    autopas_set_lock(&addParticleLock);
    _particles.push_back(m);
    autopas_unset_lock(&addParticleLock);
  }

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
   * type of the internal iterator
   */
  typedef internal::SingleCellIterator<Particle, FullParticleCell<Particle, SoAArraysType>> iterator_t;

 private:
  autopas_lock_t addParticleLock;
};

}  // namespace autopas
