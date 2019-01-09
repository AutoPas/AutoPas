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
  FullParticleCell() { autopas_init_lock(&particlesLock); }

  ~FullParticleCell() { autopas_destroy_lock(&particlesLock); }

  /**
   * Move Constructor
   * @param other
   */
  FullParticleCell(FullParticleCell&& other) {
    _particles = std::move(other._particles);
    _particleSoABuffer = std::move(other._particleSoABuffer);
    autopas_init_lock(&particlesLock);
  }

  /**
   * Copy constructor
   * @param other
   */
  FullParticleCell(const FullParticleCell& other) {
    _particles = other._particles;
    _particleSoABuffer = other._particleSoABuffer;
    autopas_init_lock(&particlesLock);
  }

  /**
   * Assignment operator
   * @param other
   * @return reference to this object after copy
   */
  FullParticleCell& operator=(FullParticleCell& other) {
    _particles = other._particles;
    _particleSoABuffer = other._particleSoABuffer;
    autopas_init_lock(&particlesLock);
    return *this;
  }

  void addParticle(Particle& m) override {
    autopas_set_lock(&particlesLock);
    _particles.push_back(m);
    autopas_unset_lock(&particlesLock);
  }

  virtual SingleCellIteratorWrapper<Particle> begin() override {
    return SingleCellIteratorWrapper<Particle>(new iterator_t(this));
  }

  unsigned long numParticles() const override { return _particles.size(); }

  /**
   * Returns a reference to the element at position n in the cell.
   * @param n Position of an element in the container
   * @return Reference to the element
   */
  Particle& operator[](size_t n) { return _particles[n]; }

  bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override { _particles.clear(); }

  void deleteByIndex(size_t index) override {
    autopas_set_lock(&particlesLock);
    assert(index < numParticles());

    if (index < numParticles() - 1) {
      std::swap(_particles[index], _particles[numParticles() - 1]);
    }
    _particles.pop_back();
    autopas_unset_lock(&particlesLock);
  }

  /**
   * Resizes the container so that it contains n elements.
   * @param n New container size
   */
  void resize(size_t n) { _particles.resize(n); }

  /**
   * Sort the particles in the cell by a dimension.
   * @param dim dimension to sort
   */
  void sortByDim(const size_t dim) {
    std::sort(_particles.begin(), _particles.end(),
              [dim](const Particle& a, const Particle& b) -> bool { return a.getR()[dim] < b.getR()[dim]; });
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
   * device AoS Buffer
   */
  double * _particlesDevice;

  /**
   * device particle SoABuffer
   */
  typename Particle::SoADevice _particleSoABufferDevice;

  /**
   * type of the internal iterator
   */
  typedef internal::SingleCellIterator<Particle, FullParticleCell<Particle, SoAArraysType>> iterator_t;

 private:
  autopas_lock_t particlesLock;
};
}  // namespace autopas
