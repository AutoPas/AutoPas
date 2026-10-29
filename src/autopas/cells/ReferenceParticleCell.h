/**
 * @file ReferenceParticleCell.h
 * @date 05.04.2020
 * @author lunaticcoding
 */

#pragma once

#include <mutex>
#include <vector>

#include "autopas/cells/ParticleCell.h"
#include "autopas/iterators/ReferenceCellIterator.h"
#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {

/**
 * This class refers to particles that are stored in a vector outside of the cell by their indices in that vector.
 * Adding particles to that vector does not invalidate the indices, removing particles from it does.
 * @tparam Particle_T
 */
template <class Particle_T>
class ReferenceParticleCell : public ParticleCell<Particle_T> {
 public:
  /**
   * The structure of the SoAs is defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

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

  void addParticle(const Particle_T &p) override {
    autopas::utils::ExceptionHandler::exception("Should use addParticleReference instead");
  }

  /**
   * Sets the vector the indices stored in this cell refer to.
   * @param particleStorage Vector holding the particles. It has to outlive this cell.
   */
  void setParticleStorage(std::vector<Particle_T> *particleStorage) { _particleStorage = particleStorage; }

  /**
   * Adds the particle at the given index of the particle storage to this cell.
   * @param index Index of the particle in the particle storage.
   */
  void addParticleReference(size_t index) {
    if (_particleStorage == nullptr) {
      autopas::utils::ExceptionHandler::exception(
          "ReferenceParticleCell::addParticleReference() called before setParticleStorage().");
    }
    const Particle_T &p = (*_particleStorage)[index];
    // sanity check that ensures that only particles of the cells OwnershipState can be added. Note: if a cell is a
    // dummy-cell, only dummies can be added, otherwise dummies can always be added
    if ((not toInt64(p.getOwnershipState() & this->_ownershipState)) and
        p.getOwnershipState() != OwnershipState::dummy) {
      autopas::utils::ExceptionHandler::exception(
          "ReferenceParticleCell::addParticleReference() can not add a particle with OwnershipState {} to a cell with "
          "OwnershipState {}",
          p.getOwnershipState(), this->_ownershipState);
    }

    std::lock_guard<AutoPasLock> guard(this->_cellLock);
    _particles.push_back(index);
  }

  /**
   * @copydoc autopas::FullParticleCell::begin()
   */
  [[nodiscard]] ReferenceCellIterator<Particle_T, true> begin() {
    return ReferenceCellIterator<Particle_T, true>(getParticleStorageBuffer(), _particles.begin());
  }

  /**
   * @copydoc autopas::FullParticleCell::begin()
   * @note const version
   */
  [[nodiscard]] ReferenceCellIterator<Particle_T, false> begin() const {
    return ReferenceCellIterator<Particle_T, false>(getParticleStorageBuffer(), _particles.cbegin());
  }

  /**
   * @copydoc autopas::FullParticleCell::end()
   */
  [[nodiscard]] ReferenceCellIterator<Particle_T, true> end() {
    return ReferenceCellIterator<Particle_T, true>(getParticleStorageBuffer(), _particles.end());
  }

  /**
   * @copydoc autopas::FullParticleCell::end()
   * @note const version
   */
  [[nodiscard]] ReferenceCellIterator<Particle_T, false> end() const {
    return ReferenceCellIterator<Particle_T, false>(getParticleStorageBuffer(), _particles.cend());
  }

  /**
   * Executes code for every particle in this cell as defined by lambda function.
   * @tparam Lambda (Particle_T &p) -> void
   * @param forEachLambda code to be executed on particles
   * @param behavior ownerships of particles that should be in-/excluded
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
    const std::array<double, 3> dummy{};
    forEachImpl<false>(forEachLambda, dummy, dummy, behavior);
  }

  /**
   * Executes code for every particle in this cell as defined by lambda function.
   * @tparam Lambda (Particle_T &p) -> void
   * @param forEachLambda code to be executed on particles
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior ownerships of particles that should be in-/excluded
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
               const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    forEachImpl<true>(forEachLambda, lowerCorner, higherCorner, behavior);
  }
  /**
   * Reduce properties of particles as defined by a lambda function.
   * @tparam Lambda (Particle_T p, A initialValue) -> void
   * @tparam A type of particle attribute to be reduced
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   * @param behavior ownerships of particles that should be in-/excluded
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior) {
    const std::array<double, 3> dummy{};
    reduceImpl<true, false>(reduceLambda, result, dummy, dummy, behavior);
  }

  /**
   * Reduce properties of particles as defined by a lambda function.
   * @tparam Lambda (Particle_T p, A initialValue) -> void
   * @tparam A type of particle attribute to be reduced
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior ownerships of particles that should be in-/excluded
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
              const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    reduceImpl<true, true>(reduceLambda, result, lowerCorner, higherCorner, behavior);
  }

  /**
   * Get the number of all particles stored in this cell (owned, halo and dummy).
   * @return number of particles stored in this cell (owned, halo and dummy).
   */
  [[nodiscard]] size_t size() const override { return _particles.size(); }

  /**
   * @copydoc autopas::ParticleCell::getNumberOfParticles()
   */
  [[nodiscard]] size_t getNumberOfParticles(IteratorBehavior behavior) const override {
    std::lock_guard<AutoPasLock> guard(this->_cellLock);
    if (_particleStorage == nullptr) {
      utils::ExceptionHandler::exception(
          "ReferenceParticleCell::getNumberOfParticles: No particle storage has been set!");
    }

    size_t numParticles{0};
    numParticles = std::count_if(_particles.begin(), _particles.end(),
                                 [&](size_t index) { return behavior.contains((*_particleStorage)[index]); });

    return numParticles;
  }

  /**
   * Returns a reference to the element at position n in the cell.
   * @param n Position of an element in the container
   * @return Reference to the element
   */
  Particle_T &operator[](size_t n) { return (*_particleStorage)[_particles[n]]; }

  /**
   * @copydoc ParticleCell::getParticleCellTypeAsEnum()
   */
  CellType getParticleCellTypeAsEnum() override { return CellType::ReferenceParticleCell; }

  /**
   * Returns a const reference to the element at position n in the cell.
   * @param n Position of an element in the container
   * @return Reference to the element
   */
  const Particle_T &operator[](size_t n) const { return (*_particleStorage)[_particles[n]]; }

  /**
   * Returns the particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  [[nodiscard]] Particle_T &at(size_t index) {
    if (_particleStorage == nullptr) {
      utils::ExceptionHandler::exception("ReferenceParticleCell::at: No particle storage has been set!");
    }
    return (*_particleStorage)[_particles.at(index)];
  }

  /**
   * Returns the const particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  [[nodiscard]] const Particle_T &at(size_t index) const {
    if (_particleStorage == nullptr) {
      utils::ExceptionHandler::exception("ReferenceParticleCell::at: No particle storage has been set!");
    }
    return (*_particleStorage)[_particles.at(index)];
  }

  [[nodiscard]] bool isEmpty() const override { return size() == 0; }

  void clear() override { _particles.clear(); }

  void deleteDummyParticles() override {
    if (_particleStorage == nullptr) {
      utils::ExceptionHandler::exception(
          "ReferenceParticleCell::deleteDummyParticles: No particle storage has been set!");
    }
    _particles.erase(std::remove_if(_particles.begin(), _particles.end(),
                                    [&](size_t index) { return (*_particleStorage)[index].isDummy(); }),
                     _particles.end());
  }

  void deleteByIndex(size_t index) override {
    std::lock_guard<AutoPasLock> lock(this->_cellLock);
    if (index >= size()) {
      utils::ExceptionHandler::exception("Index out of range (range: [0, {}[, index: {})", size(), index);
    }
    if (_particleStorage == nullptr) {
      utils::ExceptionHandler::exception("ReferenceParticleCell::deleteByIndex: No particle storage has been set!");
    }

    (*_particleStorage)[_particles[index]].setOwnershipState(OwnershipState::dummy);

    if (index < size() - 1) {
      std::swap(_particles[index], _particles[size() - 1]);
    }
    _particles.pop_back();
  }

  void setCellLength(std::array<double, 3> &cellLength) override { _cellLength = cellLength; }

  [[nodiscard]] std::array<double, 3> getCellLength() const override { return _cellLength; }

  /**
   * Sort the particles in the cell by a dimension.
   * @param dim dimension to sort
   */
  void sortByDim(const size_t dim) {
    if (_particleStorage == nullptr) {
      utils::ExceptionHandler::exception("ReferenceParticleCell::sortByDim: No particle storage has been set!");
    }
    std::sort(_particles.begin(), _particles.end(), [&](size_t a, size_t b) -> bool {
      return (*_particleStorage)[a].getR()[dim] < (*_particleStorage)[b].getR()[dim];
    });
  }

  /**
   * Requests that the vector capacity be at least enough to contain n elements.
   * @param n Minimum capacity for the vector.
   */
  void reserve(size_t n) { _particles.reserve(n); }

  /**
   * Indices of the particles of this cell in the particle storage.
   */
  std::vector<size_t> _particles;

  /**
   * SoA buffer of this cell.
   */
  SoA<SoAArraysType> _particleSoABuffer;

 private:
  std::array<double, 3> _cellLength;

  /**
   * Vector holding the particles the indices in _particles refer to.
   */
  std::vector<Particle_T> *_particleStorage{nullptr};

  /**
   * Returns the start of the particle storage's buffer, or nullptr if no storage is set.
   * @return Pointer to the first particle in the particle storage.
   */
  Particle_T *getParticleStorageBuffer() const {
    return _particleStorage == nullptr ? nullptr : _particleStorage->data();
  }

  template <bool regionCheck, typename Lambda>
  void forEachImpl(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                   const std::array<double, 3> &higherCorner,
                   IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHaloOrDummy) {
    Particle_T *particles = getParticleStorageBuffer();
    // Check the storage rather than its buffer, which may legitimately be nullptr if the storage is empty.
    if (_particleStorage == nullptr) {
      utils::ExceptionHandler::exception("ReferenceParticleCell::forEachImpl: No particle storage has been set!");
    }
    for (size_t index : _particles) {
      Particle_T &p = particles[index];
      if (behavior.contains(p)) {
        if ((not regionCheck) or utils::inBox(p.getR(), lowerCorner, higherCorner)) {
          forEachLambda(p);
        }
      }
    }
  }

  template <bool ownershipCheck, bool regionCheck, typename Lambda, typename A>
  void reduceImpl(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                  const std::array<double, 3> &higherCorner,
                  IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHaloOrDummy) {
    Particle_T *particles = getParticleStorageBuffer();
    // Check the storage rather than its buffer, which may legitimately be nullptr if the storage is empty.
    if (_particleStorage == nullptr) {
      utils::ExceptionHandler::exception("ReferenceParticleCell::reduceImpl: No particle storage has been set!");
    }
    for (size_t index : _particles) {
      Particle_T &p = particles[index];
      if ((not ownershipCheck) or behavior.contains(p)) {
        if ((not regionCheck) or utils::inBox(p.getR(), lowerCorner, higherCorner)) {
          reduceLambda(p, result);
        }
      }
    }
  }
};
}  // namespace autopas
