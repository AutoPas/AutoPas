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
 * ParticleCell concept. This is necessary to use the octree nodes as particle containers.
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
                    double interactionLength, double cellSizeFactor) {
    _pointer = std::make_unique<OctreeLeafNode<Particle>>(boxMin, boxMax, nullptr, treeSplitThreshold,
                                                          interactionLength, cellSizeFactor);
    // omp_init_lock(&_lock);
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
  }

#if 1
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
#else
  /**
   * Get an iterator to the start of a ParticleCell.
   * normal use:
   * for(auto iter = cell.begin(); iter.isValid; ++iter){...}
   * @return the iterator
   */
  SingleCellIteratorWrapper<Particle, true> begin() override {
    return SingleCellIteratorWrapper<ParticleType, true>(getIterator<true>());
  }

  /**
   * @copydoc begin()
   * @note const version
   */
  SingleCellIteratorWrapper<Particle, false> begin() const override {
    return SingleCellIteratorWrapper<ParticleType, false>(getIterator<false>());
  }
#endif

  /**
   * Get the number of particles stored in this cell.
   * @return number of particles stored in this cell
   */
  [[nodiscard]] unsigned long numParticles() const override { return _pointer->getNumParticles(); }

  /**
   * Check if the cell is not empty.
   * @return true if at least one particle is stored in this cell
   */
  [[nodiscard]] bool isNotEmpty() const override { return numParticles() > 0; }

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

#if 1
  /**
   * Type of the internal iterator.
   */
  using iterator_t = internal::SingleCellIterator<Particle, OctreeNodeWrapper<Particle>, true>;

  /**
   * Type of the internal const iterator.
   */
  using const_iterator_t = internal::SingleCellIterator<Particle, OctreeNodeWrapper<Particle>, false>;
#endif

 private:
  Particle &getFromReloadingIterator(size_t index) const {
#if 0
    // Reload the buffer only when the index is zero, this saves a little bit of compute time, since only one
    // "traversal" of the octree is required in order to gather all particles
    omp_set_lock((omp_lock_t *)&_lock);
    printf("Querying: index=%lu of %lu on T%d\n", index, ps.size(), omp_get_thread_num());
    if (index == 0 or ps.empty()) {
      ps.clear();
    }
    omp_unset_lock((omp_lock_t *)&_lock);
#endif

    // TODO(johannes): This needs severe improvement. If we just copy out all particles, the implementation becomes
    //  unsafe for threading. We need a way to iterate the octree using a better traversal idea.
    //  Referenced in https://github.com/AutoPas/AutoPas/issues/625
    std::vector<Particle *> ps;
    _pointer->appendAllParticles(ps);
    return *ps[index];
  }

  // omp_lock_t _lock;
  std::unique_ptr<OctreeNodeInterface<Particle>> _pointer;
};
}  // namespace autopas