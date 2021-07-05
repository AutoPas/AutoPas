/**
 * @file ReferenceParticleCell.h
 * @date 05.04.2020
 * @author lunaticcoding
 */

#pragma once

#include <mutex>
#include <vector>

#include "autopas/cells/ParticleCell.h"
#include "autopas/iterators/SingleCellIterator.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/CudaSoA.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class handles the storage of particles in their full form.
 * @tparam Particle
 */
template <class Particle>
class ReferenceParticleCell : public ParticleCell<Particle> {
 public:
  /**
   * The structure of the SoAs is defined by the particle.
   */
  using SoAArraysType = typename Particle::SoAArraysType;

  /**
   * Constructs a new ReferenceParticleCell.
   */
  ReferenceParticleCell()
      : _cellLength({std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()}) {}

  /**
   * Constructs a new ReferenceParticleCell with the given cell side length.
   * @param cellLength cell side length
   */
  explicit ReferenceParticleCell(const std::array<double, 3> &cellLength) : _cellLength(cellLength) {}

  void addParticle(const Particle &p) override {
    autopas::utils::ExceptionHandler::exception("Should use addParticleReference instead");
  }

  /**
   * @copydoc ParticleCell::addParticle()
   */
  void addParticleReference(Particle *p) {
    _particlesLock.lock();
    _particles.push_back(p);
    _particlesLock.unlock();
  }

  SingleCellIteratorWrapper<Particle, true> begin() override {
    return SingleCellIteratorWrapper<Particle, true>(new iterator_t(this));
  }

  SingleCellIteratorWrapper<Particle, false> begin() const override {
    return SingleCellIteratorWrapper<Particle, false>(new const_iterator_t(this));
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
    _forEach<false>(forEachLambda, nullptr, nullptr, behavior);
  }

  template <typename Lambda>
  void forEach(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
               const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    _forEach<true>(forEachLambda, lowerCorner, higherCorner, behavior);
  }

  [[nodiscard]] unsigned long numParticles() const override { return _particles.size(); }

  /**
   * Returns a reference to the element at position n in the cell.
   * @param n Position of an element in the container
   * @return Reference to the element
   */
  Particle &operator[](size_t n) { return *(_particles[n]); }

  /**
   * @copydoc ParticleCell::getParticleCellTypeAsEnum()
   */
  CellType getParticleCellTypeAsEnum() override { return CellType::ReferenceParticleCell; }

  /**
   * Returns a const reference to the element at position n in the cell.
   * @param n Position of an element in the container
   * @return Reference to the element
   */
  const Particle &operator[](size_t n) const { return *(_particles[n]); }

  /**
   * Returns the particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  [[nodiscard]] Particle &at(size_t index) { return *(_particles.at(index)); }

  /**
   * Returns the const particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  [[nodiscard]] const Particle &at(size_t index) const { return *(_particles.at(index)); }

  [[nodiscard]] bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override { _particles.clear(); }

  void deleteDummyParticles() override {
    _particles.erase(
        std::remove_if(_particles.begin(), _particles.end(), [](const auto &particle) { return particle->isDummy(); }),
        _particles.end());
  }

  void deleteByIndex(size_t index) override {
    std::lock_guard<AutoPasLock> lock(_particlesLock);
    if (index >= numParticles()) {
      utils::ExceptionHandler::exception("Index out of range (range: [0, {}[, index: {})", numParticles(), index);
    }

    _particles[index]->setOwnershipState(OwnershipState::dummy);

    if (index < numParticles() - 1) {
      std::swap(_particles[index], _particles[numParticles() - 1]);
    }
    _particles.pop_back();
  }

  void setCellLength(std::array<double, 3> &cellLength) override { _cellLength = cellLength; }

  [[nodiscard]] std::array<double, 3> getCellLength() const override { return _cellLength; }

  /**
   * Resizes the container so that it contains n elements.
   * @param n New container size
   * @param toInsert Particle to insert. This is needed to allow for non-default-constructible particles.
   */
  void resize(size_t n, const Particle &toInsert) { _particles.resize(n, std::unique_ptr<Particle>(toInsert)); }

  /**
   * Sort the particles in the cell by a dimension.
   * @param dim dimension to sort
   */
  void sortByDim(const size_t dim) {
    std::sort(_particles.begin(), _particles.end(),
              [dim](const auto *a, const auto *b) -> bool { return a->getR()[dim] < b->getR()[dim]; });
  }

  /**
   * Requests that the vector capacity be at least enough to contain n elements.
   * @param n Minimum capacity for the vector.
   */
  void reserve(size_t n) { _particles.reserve(n); }

  /**
   * Storage of the molecules of the cell.
   */
  std::vector<Particle *> _particles;

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
  using iterator_t = internal::SingleCellIterator<Particle, ReferenceParticleCell<Particle>, true>;

  /**
   * Type of the internal const iterator.
   */
  using const_iterator_t = internal::SingleCellIterator<Particle, ReferenceParticleCell<Particle>, false>;

 private:
  AutoPasLock _particlesLock;
  std::array<double, 3> _cellLength;

  template <bool regionCheck, typename Lambda>
  void _forEach(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                const std::array<double, 3> &higherCorner,
                IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHaloOrDummy) {
    auto isParticleInRegion = [&](Particle &p) -> bool { return utils::inBox(p.getR(), lowerCorner, higherCorner); };

    auto isParticleValid = [&](Particle &p) -> bool {
      switch (behavior) {
        case options::IteratorBehavior::ownedOrHaloOrDummy:
          return true;
        case options::IteratorBehavior::ownedOrHalo:
          return not p.isDummy();
        case options::IteratorBehavior::halo:
          return p.isHalo();
        case options::IteratorBehavior::owned:
          return p.isOwned();
        default:
          utils::ExceptionHandler::exception("unknown iterator behavior");
          return false;
      }
    };

    for (Particle *p : _particles) {
      if (isParticleValid(*p)) {
        if (regionCheck & isParticleInRegion(*p)) {
          forEachLambda(*p);
        }
      }
    }
  }
};
}  // namespace autopas
