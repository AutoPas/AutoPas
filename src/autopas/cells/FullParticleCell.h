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
#include "autopas/options/IteratorBehavior.h"
#include "autopas/utils/SoA.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

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
                     std::numeric_limits<double>::max()}) {}

  /**
   * Constructs a new FullParticleCell with the given cell side length.
   * @param cellLength cell side length
   */
  explicit FullParticleCell(const std::array<double, 3> &cellLength) : _cellLength(cellLength) {}

  /**
   * @copydoc ParticleCell::addParticle()
   */
  void addParticle(const Particle &p) override {
    particlesLock.lock();
    _particles.push_back(p);
    particlesLock.unlock();
  }

  SingleCellIteratorWrapper<Particle, true> begin() override {
    return SingleCellIteratorWrapper<Particle, true>(new iterator_t(this));
  }

  SingleCellIteratorWrapper<Particle, false> begin() const override {
    return SingleCellIteratorWrapper<Particle, false>(new const_iterator_t(this));
  }

  /**
   * executes code for every particle in this cell as defined by lambda function
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on particles
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda) {
    const std::array<double, 3> dummy{};
    _forEach<false, false>(forEachLambda, dummy, dummy);
  }

  /**
   * executes code for every particle in this cell as defined by lambda function
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on particles
   * @param behavior ownerships of particles that should be in-/excluded
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior) {
    const std::array<double, 3> dummy{};
    _forEach<true, false>(forEachLambda, dummy, dummy, behavior);
  }

  /**
   * executes code for every particle in this cell as defined by lambda function
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on particles
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior ownerships of particles that should be in-/excluded
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
               const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    _forEach<true, true>(forEachLambda, lowerCorner, higherCorner, behavior);
  }

  /**
   * reduce properties of particles as defined by a lambda function
   * @tparam Lambda (Particle p, A initialValue) -> void
   * @tparam A type of particle attribute to be reduced
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result) {
    const std::array<double, 3> dummy{};
    _reduce<false, false>(reduceLambda, result, dummy, dummy);
  }

  /**
   * reduce properties of particles as defined by a lambda function
   * @tparam Lambda (Particle p, A initialValue) -> void
   * @tparam A type of particle attribute to be reduced
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   * @param behavior ownerships of particles that should be in-/excluded
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior) {
    const std::array<double, 3> dummy{};
    _reduce<true, false>(reduceLambda, result, dummy, dummy, behavior);
  }

  /**
   * reduce properties of particles as defined by a lambda function
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
    _reduce<true, true>(reduceLambda, result, lowerCorner, higherCorner, behavior);
  }

  [[nodiscard]] unsigned long numParticles() const override { return _particles.size(); }

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
  Particle &at(size_t index) { return _particles.at(index); }

  /**
   * @copydoc ParticleCell::getParticleCellTypeAsEnum()
   */
  CellType getParticleCellTypeAsEnum() override { return CellType::FullParticleCell; }

  /**
   * Returns the const particle at position index. Needed by SingleCellIterator.
   * @param index the position of the particle to return.
   * @return the particle at position index.
   */
  const Particle &at(size_t index) const { return _particles.at(index); }

  [[nodiscard]] bool isNotEmpty() const override { return numParticles() > 0; }

  void clear() override { _particles.clear(); }

  void deleteDummyParticles() override {
    _particles.erase(
        std::remove_if(_particles.begin(), _particles.end(), [](const auto &particle) { return particle.isDummy(); }),
        _particles.end());
  }

  void deleteByIndex(size_t index) override {
    std::lock_guard<AutoPasLock> lock(particlesLock);
    if (index >= numParticles()) {
      utils::ExceptionHandler::exception("Index out of range (range: [0, {}[, index: {})", numParticles(), index);
    }

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
  void resize(size_t n, const Particle &toInsert) { _particles.resize(n, toInsert); }

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

  inline bool isParticleInRegion(Particle &p, const std::array<double, 3> &lowerCorner,
                                 const std::array<double, 3> &higherCorner) {
    return utils::inBox(p.getR(), lowerCorner, higherCorner);
  }

  inline bool isParticleValid(Particle &p, IteratorBehavior behavior) {
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
  }

  template <bool ownershipCheck, bool regionCheck, typename Lambda>
  void _forEach(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                const std::array<double, 3> &higherCorner,
                IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHaloOrDummy) {
    for (Particle &p : _particles) {
      if ((not ownershipCheck) or isParticleValid(p, behavior)) {
        if ((not regionCheck) or isParticleInRegion(p, lowerCorner, higherCorner)) {
          forEachLambda(p);
        }
      }
    }
  }

  template <bool ownershipCheck, bool regionCheck, typename Lambda, typename A>
  void _reduce(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                const std::array<double, 3> &higherCorner,
                IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHaloOrDummy) {
    for (Particle &p : _particles) {
      if ((not ownershipCheck) or isParticleValid(p, behavior)) {
        if ((not regionCheck) or isParticleInRegion(p, lowerCorner, higherCorner)) {
          reduceLambda(p, result);
        }
      }
    }
  }
};
}  // namespace autopas
