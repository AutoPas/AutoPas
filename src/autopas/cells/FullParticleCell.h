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
#include "autopas/iterators/CellIterator.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/inBox.h"

namespace autopas {

/**
 * This class handles the storage of particles in their full form.
 * @tparam Particle_T
 */
template <class Particle_T>
class FullParticleCell : public ParticleCell<Particle_T> {
 public:
  /**
   * The structure of the SoAs is defined by the particle.
   */
  using SoAArraysType = typename Particle_T::SoAArraysType;

  /**
   * Type that holds or refers to the actual particles.
   */
  using StorageType = std::vector<Particle_T>;

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
  explicit FullParticleCell(const std::array<double, 3> &cellLength) : _cellLength(cellLength) {}

  void addParticle(const Particle_T &p) override {
    std::lock_guard<AutoPasLock> guard(this->_cellLock);

    // sanity check that ensures that only particles of the cells OwnershipState can be added. Note: if a cell is a
    // dummy-cell, only dummies can be added, otherwise dummies can always be added
    if ((not toInt64(p.getOwnershipState() & this->_ownershipState)) and
        p.getOwnershipState() != OwnershipState::dummy) {
      autopas::utils::ExceptionHandler::exception(
          "FullParticleCell::addParticle() can not add a particle with OwnershipState {} to a cell with OwnershipState "
          "{}",
          p.getOwnershipState(), this->_ownershipState);
    }

    _particles.push_back(p);
  }

  /**
   * Get an iterator to the start of a ParticleCell.
   * normal use:
   * for(auto iter = cell.begin(); iter.isValid; ++iter){...}
   * @return the iterator
   */
  [[nodiscard]] CellIterator<StorageType, true> begin() { return CellIterator<StorageType, true>(_particles.begin()); }

  /**
   * @copydoc autopas::FullParticleCell::begin()
   * @note const version
   */
  [[nodiscard]] CellIterator<StorageType, false> begin() const {
    return CellIterator<StorageType, false>(_particles.cbegin());
  }

  /**
   * Get an iterator to the end of a ParticleCell.
   * normal use:
   * for(auto &p : cell){...}
   * @return the iterator
   */
  [[nodiscard]] CellIterator<StorageType, true> end() { return CellIterator<StorageType, true>(_particles.end()); }

  /**
   * @copydoc autopas::FullParticleCell::end()
   * @note const version
   */
  [[nodiscard]] CellIterator<StorageType, false> end() const {
    return CellIterator<StorageType, false>(_particles.cend());
  }

  /**
   * Executes code for every particle in this cell as defined by lambda function.
   * @tparam Lambda (Particle_T &p) -> void
   * @param forEachLambda code to be executed on particles
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda) {
    const std::array<double, 3> dummy{};
    forEachImpl<false, false>(forEachLambda, dummy, dummy);
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
    forEachImpl<true, false>(forEachLambda, dummy, dummy, behavior);
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
    forEachImpl<true, true>(forEachLambda, lowerCorner, higherCorner, behavior);
  }

  /**
   * Reduce properties of particles as defined by a lambda function.
   * @tparam Lambda (Particle_T p, A initialValue) -> void
   * @tparam A type of particle attribute to be reduced
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result) {
    const std::array<double, 3> dummy{};
    reduceImpl<false, false>(reduceLambda, result, dummy, dummy);
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
    return std::count_if(_particles.begin(), _particles.end(), [&behavior](auto p) { return behavior.contains(p); });
  }

  /**
   * Returns a reference to the element at position n in the cell.
   * @param n Position of an element in the container
   * @return Reference to the element
   */
  Particle_T &operator[](size_t n) { return _particles[n]; }

  /**
   * Returns a const reference to the element at position n in the cell.
   * @param n Position of an element in the container
   * @return Reference to the element
   */
  const Particle_T &operator[](size_t n) const { return _particles[n]; }

  /**
   * Returns the particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  Particle_T &at(size_t index) { return _particles.at(index); }

  CellType getParticleCellTypeAsEnum() override { return CellType::FullParticleCell; }

  /**
   * Returns the const particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  const Particle_T &at(size_t index) const { return _particles.at(index); }

  [[nodiscard]] bool isEmpty() const override { return size() == 0; }

  void clear() override { _particles.clear(); }

  void deleteDummyParticles() override {
    _particles.erase(
        std::remove_if(_particles.begin(), _particles.end(), [](const auto &particle) { return particle.isDummy(); }),
        _particles.end());
  }

  void deleteByIndex(size_t index) override {
    std::lock_guard<AutoPasLock> lock(this->_cellLock);
    if (index >= size()) {
      utils::ExceptionHandler::exception("Index out of range (range: [0, {}[, index: {})", size(), index);
    }

    if (index < size() - 1) {
      std::swap(_particles[index], _particles[size() - 1]);
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
  void resize(size_t n, const Particle_T &toInsert) { _particles.resize(n, toInsert); }

  /**
   * Sort the particles in the cell by a dimension.
   * @param dim dimension to sort
   */
  void sortByDim(const size_t dim) {
    std::sort(_particles.begin(), _particles.end(),
              [dim](const Particle_T &a, const Particle_T &b) -> bool { return a.getR()[dim] < b.getR()[dim]; });
  }

  /**
   * Requests that the vector capacity be at least enough to contain n elements.
   * @param n Minimum capacity for the vector.
   */
  void reserve(size_t n) { _particles.reserve(n); }

  /**
   * Storage of the molecules of the cell.
   */
  StorageType _particles{};

  /**
   * SoA buffer of this cell.
   */
  SoA<SoAArraysType> _particleSoABuffer{};

 private:
  std::array<double, 3> _cellLength;

  template <bool ownershipCheck, bool regionCheck, typename Lambda>
  void forEachImpl(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                   const std::array<double, 3> &higherCorner,
                   IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHaloOrDummy) {
    for (Particle_T &p : _particles) {
      if ((not ownershipCheck) or behavior.contains(p)) {
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
    for (Particle_T &p : _particles) {
      if ((not ownershipCheck) or behavior.contains(p)) {
        if ((not regionCheck) or utils::inBox(p.getR(), lowerCorner, higherCorner)) {
          reduceLambda(p, result);
        }
      }
    }
  }
};
}  // namespace autopas
