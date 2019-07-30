/**
 * @file FullParticleCell.h
 * @date 18.01.2018
 * @author seckler
 */

#pragma once

#include <vector>
#include "autopas/cells/ParticleCell.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/utils/CudaSoA.h"
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
  /**
   * Constructs a new FullParticleCell.
   */
  FullParticleCell()
      : _cellLength({std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()}) {}

  /**
   * Constructs a new FullParticleCell with the given cell side length.
   * @param cellLength cell side length
   */
  FullParticleCell(std::array<double, 3> &cellLength) : _cellLength(cellLength) {}

  /**
   * @copydoc ParticleCell::addParticle()
   */
  void addParticle(const Particle &p) override {
    particlesLock.lock();
    _particles.push_back(p);
    particlesLock.unlock();
  }

  SingleCellIteratorWrapper<Particle> begin() override {
    return SingleCellIteratorWrapper<Particle>(new iterator_t(this));
  }

  unsigned long numParticles() const override { return _particles.size(); }

  /**
   * Returns a reference to the element at position n in the cell.
   * @param n Position of an element in the container
   * @return Reference to the element
   */
  Particle &operator[](size_t n) { return _particles[n]; }

  /**
   * Returns a const reference to the element at position n in the cell.
   * @param n Position of an element in the container
   * @return Reference to the element
   */
  const Particle &operator[](size_t n) const { return _particles[n]; }

  bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override { _particles.clear(); }

  void deleteByIndex(size_t index) override {
    particlesLock.lock();
    assert(index < numParticles());

    if (index < numParticles() - 1) {
      std::swap(_particles[index], _particles[numParticles() - 1]);
    }
    _particles.pop_back();
    particlesLock.unlock();
  }

  void setCellLength(const std::array<double, 3> &cellLength) override { _cellLength = cellLength; }

  std::array<double, 3> getCellLength() const override { return _cellLength; }

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
              [dim](const Particle &a, const Particle &b) -> bool { return a.getR()[dim] < b.getR()[dim]; });
  }

  /**
   * Requests that the vector capacity be at least enough to contain n elements.
   * @param n Minimum capacity for the vector.
   */
  void reserve(size_t n) { _particles.reserve(n); }

  /**
   * Storage of the molecules of the cell.
   */
  std::vector<Particle> _particles;

  /**
   * SoA buffer of this cell.
   */
  SoA<SoAArraysType> _particleSoABuffer;

  /**
   * Device particle SoABuffer.
   */
  CudaSoA<typename Particle::CudaDeviceArraysType> _particleSoABufferDevice;

  /**
   * Type of the internal iterator.
   */
  typedef internal::SingleCellIterator<Particle, FullParticleCell<Particle, SoAArraysType>> iterator_t;

 private:
  AutoPasLock particlesLock;
  std::array<double, 3> _cellLength;
};
}  // namespace autopas
