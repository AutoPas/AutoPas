/**
 * @file OctreeNodeWrapper.h
 *
 * @author Johannes Spies
 * @date 03.05.2021
 */
#pragma once

#include <memory>

#include "autopas/containers/octree/OctreeNodeInterface.h"

namespace autopas {
// TODO: Why do you exist?
template <typename Particle>
class OctreeNodeWrapper : public ParticleCell<Particle> {
 public:
  using ParticleCell = OctreeNodeWrapper<Particle>;
  using ParticleType = typename ParticleCell::ParticleType;

  OctreeNodeWrapper(std::array<double, 3> boxMin, std::array<double, 3> boxMax) {
    _pointer = std::make_unique<OctreeLeafNode<Particle>>(boxMin, boxMax, nullptr);
  }

  void appendAllParticles(std::vector<Particle> &ps) { _pointer->appendAllParticles(ps); }

  /**
   * Adds a Particle to the cell.
   * @param p the particle to be added
   */
  void addParticle(const Particle &p) override { _pointer->insert(_pointer, p); }

  /**
   * Get an iterator to the start of a ParticleCell.
   * normal use:
   * for(auto iter = cell.begin(); iter.isValid; ++iter){...}
   * @return the iterator
   */
  SingleCellIteratorWrapper<Particle, true> begin() override {
    return SingleCellIteratorWrapper<ParticleType, true>(new iterator_t(this));
  }

  /**
   * @copydoc begin()
   * @note const version
   */
  SingleCellIteratorWrapper<Particle, false> begin() const override {
    return SingleCellIteratorWrapper<ParticleType, false>(new const_iterator_t(this));
  }

  /**
   * Get the number of particles stored in this cell.
   * @return number of particles stored in this cell
   */
  [[nodiscard]] unsigned long numParticles() const { return _pointer->getNumParticles(); }

  /**
   * Check if the cell is not empty.
   * @return true if at least one particle is stored in this cell
   */
  [[nodiscard]] bool isNotEmpty() const { return numParticles() > 0; }

  /**
   * Deletes all particles in this cell.
   */
  void clear() override { _pointer->clearChildren(_pointer); }

  /**
   * Deletes all dummy particles in this cell.
   */
  void deleteDummyParticles() override {}

  /**
   * Get the ParticleCell type as an ParticleCellTypeEnum
   * @return The Cell type as an Enum
   */
  CellType getParticleCellTypeAsEnum() override { return CellType::FullParticleCell; }

  /**
   * Deletes the index-th particle.
   * @param index the index of the particle that shall be deleted
   */
  void deleteByIndex(size_t index) override {}

  /**
   * Set the side lengths of this cell.
   * @param cellLength cell side length
   */
  void setCellLength(std::array<double, 3> &cellLength) override {}

  /**
   * Get the side lengths of this cell.
   * @return cell side length
   */
  [[nodiscard]] std::array<double, 3> getCellLength() const override {
    std::array<double, 3> result = utils::ArrayMath::sub(_pointer->getBoxMin(), _pointer->getBoxMax());
    return result;
  }

  Particle &at(size_t index) {
    std::vector<Particle> ps;
    _pointer->appendAllParticles(ps);
    return ps[index];
  }

  const Particle &at(size_t index) const {
    std::vector<Particle> ps;
    _pointer->appendAllParticles(ps);
    return ps[index];
  }

  /**
   * Type of the internal iterator.
   */
  using iterator_t = internal::SingleCellIterator<Particle, OctreeNodeWrapper<Particle>, true>;

  /**
   * Type of the internal const iterator.
   */
  using const_iterator_t = internal::SingleCellIterator<Particle, OctreeNodeWrapper<Particle>, false>;

 private:
  std::unique_ptr<OctreeNodeInterface<Particle>> _pointer;
};
}  // namespace autopas