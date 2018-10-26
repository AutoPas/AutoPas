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

namespace autopas {

/**
 * This class handles the storage of particles in their full form.
 * @tparam Particle
 */
template <class Particle, class SoAArraysType = typename Particle::SoAArraysType>
class FullParticleCell : public ParticleCell<Particle> {
 public:
  FullParticleCell() {}

  void addParticle(Particle& m) override { _particles.push_back(m); }

  virtual SingleCellIteratorWrapper<Particle> begin() override {
    return SingleCellIteratorWrapper<Particle>(new iterator_t(this));
  }

  unsigned long numParticles() const override { return _particles.size(); }

  /**
   * Returns a reference to the element at position n in the container
   * @param Position of an element in the container
   * @return Reference to the element
   */
  Particle& operator[](size_t idx) { return _particles[idx]; }

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
   * Resizes the container so that it contains n elements.
   * @param n New container size
   */
  void resize(size_t n) { _particles.resize(n); }

  /**
   * Sort the particles in the cell by their z-dimension
   */
  void sortByZ() {
    std::sort(_particles.begin(), _particles.end(),
              [](const Particle& a, const Particle& b) -> bool { return a.getR()[2] < b.getR()[2]; });
  }

  /**
   * Requests that the vector capacity be at least enough to contain n elements.
   * @param n Minimum capacity for the vector.
   */
  void reserve(size_t n) { _particles.reserve(n); }

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
