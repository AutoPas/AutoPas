/**
 * @file ReferenceParticleCell.h
 * @date 05.04.2020
 * @author lunaticcoding
 */

#pragma once

#include <mutex>
#include <vector>

#include "autopas/cells/ParticleCell.h"
#include "autopas/particles/OwnershipState.h"
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
   * Type that holds or refers to the actual particles.
   */
  using StorageType = std::vector<Particle *>;

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

  /**
   * @copydoc autopas::FullParticleCell::begin()
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
   * @copydoc autopas::FullParticleCell::end()
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
   * @tparam Lambda (Particle &p) -> void
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
   * @tparam Lambda (Particle &p) -> void
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
   * @tparam Lambda (Particle p, A initialValue) -> void
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
   * @tparam Lambda (Particle p, A initialValue) -> void
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

  [[nodiscard]] bool isEmpty() const override { return numParticles() == 0; }

  void clear() override { _particles.clear(); }

  void deleteDummyParticles() override {
    _particles.erase(
        std::remove_if(_particles.begin(), _particles.end(), [](const auto &particle) { return particle->isDummy(); }),
        _particles.end());
  }

  void deleteHaloParticles() override {
    _particles.erase(
        std::remove_if(_particles.begin(), _particles.end(), [](const auto &particle) { return particle->isHalo(); }),
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
  StorageType _particles;

  /**
   * SoA buffer of this cell.
   */
  SoA<SoAArraysType> _particleSoABuffer;

 private:
  AutoPasLock _particlesLock;
  std::array<double, 3> _cellLength;

  template <bool regionCheck, typename Lambda>
  void forEachImpl(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                   const std::array<double, 3> &higherCorner,
                   IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHaloOrDummy) {
    for (Particle *p : _particles) {
      if (behavior.contains(*p)) {
        if ((not regionCheck) or utils::inBox(p->getR(), lowerCorner, higherCorner)) {
          forEachLambda(*p);
        }
      }
    }
  }

  template <bool ownershipCheck, bool regionCheck, typename Lambda, typename A>
  void reduceImpl(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                  const std::array<double, 3> &higherCorner,
                  IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHaloOrDummy) {
    for (Particle *p : _particles) {
      if ((not ownershipCheck) or behavior.contains(*p)) {
        if ((not regionCheck) or utils::inBox(p->getR(), lowerCorner, higherCorner)) {
          reduceLambda(*p, result);
        }
      }
    }
  }
};
}  // namespace autopas
