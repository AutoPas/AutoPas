/**
 * @file FullParticleCell.h
 * @date 18.01.2018
 * @author seckler
 */

#pragma once

#include <array>
#include <mutex>
#include <vector>

#include "autopas/cells/ParticleCell.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/utils/CudaSoA.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

#include <Kokkos_Core.hpp>

namespace autopas {

/**
 * This class handles the storage of particles in their full form.
 * @tparam Particle
 */
template <class Particle>
class FullParticleCell : public ParticleCell<Particle> {
 public:
  /**
   * The structure of the SoAs is defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Constructs a new FullParticleCell.
   */
  FullParticleCell()
      : _cellLength({std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()}), capacity(128), size(0),
        _particles(Kokkos::View<Particle *> ("particles device", 128)){}

  /**
   * Constructs a new FullParticleCell with the given cell side length.
   * @param cellLength cell side length
   */
  explicit FullParticleCell(const std::array<double, 3> &cellLength) : _cellLength(cellLength),
  capacity(128), size(0), _particles(Kokkos::View<Particle *> ("particles device", 128)){}

  /**
   * @copydoc ParticleCell::addParticle()
   */
  void addParticle(const Particle &p) override {
    particlesLock.lock();
    printf("Add Particle s:%u c:%u x:%f y:%f z:%f \n", numParticles(), capacity, p.getR()[0], p.getR()[1], p.getR()[2]);
    if(size == capacity) {
      capacity *= 2;
      Kokkos::resize(_particles, capacity);
    }
    _particles[size] = p;
    size++;

    particlesLock.unlock();
  }

  SingleCellIteratorWrapper<Particle, true> begin() override {
    return SingleCellIteratorWrapper<Particle, true>(new iterator_t(this));
  }

  SingleCellIteratorWrapper<Particle, false> begin() const override {
    return SingleCellIteratorWrapper<Particle, false>(new const_iterator_t(this));
  }

  [[nodiscard]] unsigned long numParticles() const override { return size; }

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

  /**
   * Returns the particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  Particle &at(size_t index) { return _particles[index]; }

  /**
   * @copydoc ParticleCell::getParticleCellTypeAsEnum()
   */
  CellType getParticleCellTypeAsEnum() override { return CellType::FullParticleCell; }

  /**
   * Returns the const particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  const Particle &at(size_t index) const { return _particles[index]; }

  [[nodiscard]] bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override { size = 0; }

  void deleteDummyParticles() override {
//    _particles.erase(
//        std::remove_if(_particles.begin(), _particles.end(), [](const auto &particle) { return particle.isDummy(); }),
//        _particles.end());
  }

  void deleteByIndex(size_t index) override {
    std::lock_guard<AutoPasLock> lock(particlesLock);
    if (index >= numParticles()) {
      utils::ExceptionHandler::exception("Index out of range (range: [0, {}[, index: {})", numParticles(), index);
    }

    if (index < numParticles() - 1) {
//      std::swap(_particles[index], _particles[numParticles() - 1]);
      _particles[index] = _particles[numParticles() - 1];
    }
    size--;
  }

  void setCellLength(std::array<double, 3> &cellLength) override { _cellLength = cellLength; }

  [[nodiscard]] std::array<double, 3> getCellLength() const override { return _cellLength; }

  /**
   * Resizes the container so that it contains n elements.
   * @param n New container size
   * @param toInsert Particle to insert. This is needed to allow for non-default-constructible particles.
   */
  void resize(size_t n, const Particle &toInsert) { _particles.resize(n, toInsert); }

  /**
   * Sort the particles in the cell by a dimension.
   * @param dim dimension to sort
   */
  void sortByDim(const size_t dim) {
    //noop atm TODO lgaertner
//    std::sort(_particles.begin(), _particles.end(),
//              [dim](const Particle &a, const Particle &b) -> bool { return a.getR()[dim] < b.getR()[dim]; });
  }

  /**
   * Requests that the vector capacity be at least enough to contain n elements.
   * @param n Minimum capacity for the vector.
   */
  void reserve(size_t n) { _particles.reserve(n); }

  /**
   * Storage of the molecules of the cell.
   */
//  std::vector<Particle> _particles;
  Kokkos::View<Particle *> _particles;

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
  using iterator_t = internal::SingleCellIterator<Particle, FullParticleCell<Particle>, true>;

  /**
   * Type of the internal const iterator.
   */
  using const_iterator_t = internal::SingleCellIterator<Particle, FullParticleCell<Particle>, false>;

 private:
  AutoPasLock particlesLock;
  std::array<double, 3> _cellLength;

  size_t capacity;
  size_t size;
};
}  // namespace autopas
