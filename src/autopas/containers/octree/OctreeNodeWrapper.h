/**
 * @file OctreeNodeWrapper.h
 *
 * @author Johannes Spies
 * @date 03.05.2021
 */
#pragma once

#include <memory>

#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/iterators/ParticleIteratorWrapper.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {
/**
 * This class wraps the functionality provided by the octree leaves and inner nodes in a structure that adheres to the
 * ParticleCell concept.
 *
 * What this wrapper does is the following: It _hides implementation details_ of the octree nodes that should not be
 * exposed to the outside. (For instance, the `insert` method of the `OctreeNodeInterface<Particle>` has a very specific
 * method signature that requires the caller to change the pointer to a subtree if necessary. Since this functionality
 * should be transparent for the user, it is wrapped in this class inside of the `addParticle` method. This method does
 * not have to be treated special like `insert`.)
 *
 * This class includes some proxy methods: `appendAllParticles`, `appendAllLeaves` and some more. Those methods only
 * invoke similar functions on the pointer to an octree. This indirection could be removed by implementing the
 * `OctreeNodeInterface<Particle>`. However, if this class inherited directly from `OctreeNodeInterface<Particle>`, it
 * would have to implement the _entire interface_ provided by `OctreeNodeInterface<Particle>`. This does not increase
 * the code quality, since one would have to implement other interface methods (like `insert`) that should not be
 * exposed to the outside. Having proxy calls inside this class keeps the interface clean.
 *
 * @tparam Particle The particle class that should be used in the octree cell.
 */
template <typename Particle>
class OctreeNodeWrapper : public ParticleCell<Particle> {
 public:
  /**
   * The contained particle cell.
   */
  using ParticleCell = OctreeNodeWrapper<Particle>;
  /**
   * The contained particle type.
   */
  using ParticleType = typename ParticleCell::ParticleType;

  /**
   * Constructs a new, empty octree and stores the root.
   *
   * @param boxMin The min coordinate of the box containing the octree
   * @param boxMax The max coordinate of the box containing the octree
   * @param treeSplitThreshold Maximum number of particles inside a leaf before it tries to split itself
   * @param interactionLength The minimum distance at which a force is considered nonzero, cutoff+skin.
   * @param cellSizeFactor The cell size factor
   */
  OctreeNodeWrapper(std::array<double, 3> boxMin, std::array<double, 3> boxMax, int unsigned treeSplitThreshold,
                    double interactionLength, double cellSizeFactor)
      : enclosedParticleCount(0) {
    _pointer = std::make_unique<OctreeLeafNode<Particle>>(boxMin, boxMax, nullptr, treeSplitThreshold,
                                                          interactionLength, cellSizeFactor);
  }

  /**
   * Append all particles in the octree to a list using DFS.
   * @param ps The list to which the particles should be appended to
   */
  void appendAllParticles(std::vector<Particle *> &ps) { _pointer->appendAllParticles(ps); }

  /**
   * Append all leaves in the octree to a list.
   * @param leaves The list to which the leaves should be appended to
   */
  void appendAllLeaves(std::vector<OctreeLeafNode<Particle> *> &leaves) { _pointer->appendAllLeaves(leaves); }

  /**
   * Adds a Particle to the cell.
   * @param p the particle to be added
   */
  void addParticle(const Particle &p) override {
    auto ret = _pointer->insert(p);
    if (ret) _pointer = std::move(ret);

    ++enclosedParticleCount;
  }

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
  [[nodiscard]] unsigned long numParticles() const override { return enclosedParticleCount; }

  /**
   * Check if the cell is not empty.
   * @return true if at least one particle is stored in this cell
   */
  [[nodiscard]] bool isNotEmpty() const override { return enclosedParticleCount > 0; }

  /**
   * Deletes all particles in this cell.
   */
  void clear() override {
    _pointer->clearChildren(_pointer);
    enclosedParticleCount = 0;
  }

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
  void deleteByIndex(size_t index) override {
    throw std::runtime_error("[OctreeNodeWrapper.h] Operation not supported");
  }

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

  /**
   * Get a particle from the iterator
   *
   * @param index The index of the particle
   * @return A ref to a particle
   */
  Particle &at(size_t index) { return getFromReloadingIterator(index); }

  /**
   * Get a particle from the iterator
   *
   * @param index The index of the particle
   * @return A const ref to a particle
   */
  const Particle &at(size_t index) const { return getFromReloadingIterator(index); }

  /**
   * Find all leaves below this subtree that are in the given range.
   * @param min The minimum coordinate in 3D space of the query area
   * @param max The maximum coordinate in 3D space of the query area
   * @return A set of all leaf nodes that are in the query region
   */
  std::set<OctreeLeafNode<Particle> *> getLeavesInRange(std::array<double, 3> min, std::array<double, 3> max) {
    return _pointer->getLeavesInRange(min, max);
  }

  /**
   * Get a raw pointer to the enclosed cell. This should only be used for debugging and insight into the internal
   * structure fo the octree.
   * @return A raw C pointer to the enclosed node
   */
  OctreeNodeInterface<Particle> *getRaw() const { return _pointer.get(); }

  /**
   * Type of the internal iterator.
   */
  using iterator_t = internal::SingleCellIterator<Particle, OctreeNodeWrapper<Particle>, true>;

  /**
   * Type of the internal const iterator.
   */
  using const_iterator_t = internal::SingleCellIterator<Particle, OctreeNodeWrapper<Particle>, false>;

 private:
  Particle &getFromReloadingIterator(size_t index) const {
    lock.lock();
    if (ps.empty() || index == 0) {
      // printf("Reloading particles (ps.empty()=%s, index=%ld) %ld times\n", ps.empty() ? "true" : "false", index,
      // count++);
      ps.clear();
      _pointer->appendAllParticles(ps);
    }
    lock.unlock();
    return *ps[index];
  }

  std::unique_ptr<OctreeNodeInterface<Particle>> _pointer;
  mutable std::vector<Particle *> ps;
  mutable AutoPasLock lock;
  mutable long count = 0;

  long enclosedParticleCount;
};
}  // namespace autopas