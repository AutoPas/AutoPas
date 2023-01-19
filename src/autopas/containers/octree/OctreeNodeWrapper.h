/**
 * @file OctreeNodeWrapper.h
 *
 * @author Johannes Spies
 * @date 03.05.2021
 */
#pragma once

#include <memory>

#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/containers/octree/OctreeStaticNodeSelector.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"

namespace autopas {
/**
 * This class wraps the functionality provided by the octree leaves and inner nodes in a structure that adheres to the
 * ParticleCell concept.
 *
 * What this wrapper does is the following: It _hides implementation details_ of the octree nodes that should not be
 * exposed to the outside. (For instance, the `OctreeNodeInterface::insert()` method has a very
 * specific method signature that requires the caller to change the pointer to a subtree if necessary. Since the user
 * should not care about this, it is wrapped in this class inside of the `addParticle()` method. This method does not
 * have to be treated special like `insert()`.)
 *
 * This class includes some proxy methods: `appendAllParticles()`, `appendAllLeaves()` and some more. Those methods only
 * invoke similar functions on the pointer to an octree. This indirection could be removed by implementing the
 * `OctreeNodeInterface<Particle>`. However, if this class inherited directly from `OctreeNodeInterface<Particle>`, it
 * would have to implement the _entire interface_ provided by `OctreeNodeInterface<Particle>`. This does not increase
 * the code quality, since one would have to implement other interface methods (like `insert()`) that should not be
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
   * Type that holds or refers to the actual particles.
   */
  using StorageType = std::vector<Particle *>;

  /**
   * Constructs a new, empty octree and stores the root.
   *
   * @param boxMin The min coordinate of the box containing the octree
   * @param boxMax The max coordinate of the box containing the octree
   * @param treeSplitThreshold Maximum number of particles inside a leaf before it tries to split itself
   * @param interactionLength The minimum distance at which a force is considered nonzero, cutoff+skin.
   * @param cellSizeFactor The cell size factor
   */
  OctreeNodeWrapper(std::array<double, 3> const &boxMin, std::array<double, 3> const &boxMax,
                    int unsigned const treeSplitThreshold, double const interactionLength,
                    double const cellSizeFactor) {
    _pointer = std::make_unique<OctreeLeafNode<Particle>>(boxMin, boxMax, nullptr, treeSplitThreshold,
                                                          interactionLength, cellSizeFactor);
  }

  /**
   * Append all particles in the octree to a list using DFS.
   * @param ps The list to which the particles should be appended to
   */
  void collectAllParticles(StorageType &ps) {
    std::lock_guard<AutoPasLock> lock(_lock);
    _pointer->collectAllParticles(ps);
  }

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
    std::lock_guard<AutoPasLock> lock(_lock);
    auto ret = _pointer->insert(p);
    if (ret) _pointer = std::move(ret);

    ++_enclosedParticleCount;
  }

  /**
   * Get an iterator to the start of a ParticleCell.
   * normal use:
   * for(auto iter = cell.begin(); iter.isValid; ++iter){...}
   * @return the iterator
   */
  CellIterator<StorageType, true> begin() {
    std::lock_guard<AutoPasLock> lock(_lock);
    _ps.clear();
    _pointer->collectAllParticles(_ps);
    return CellIterator<StorageType, true>(_ps.begin());
  }

  /**
   * @copydoc begin()
   * @note const version
   */
  CellIterator<StorageType, false> begin() const {
    std::lock_guard<AutoPasLock> lock(_lock);
    _ps.clear();
    _pointer->collectAllParticles(_ps);
    return CellIterator<StorageType, false>(_ps.cbegin());
  }

  /**
   * @copydoc autopas::FullParticleCell::end()
   */
  CellIterator<StorageType, true> end() { return CellIterator<StorageType, true>(_ps.end()); }
  /**
   * @copydoc autopas::FullParticleCell::end()
   * @note const version
   */
  CellIterator<StorageType, false> end() const { return CellIterator<StorageType, false>(_ps.end()); }

  /**
   * Get the number of particles stored in this cell.
   * @return number of particles stored in this cell
   */
  [[nodiscard]] unsigned long numParticles() const override {
    std::lock_guard<AutoPasLock> lock(_lock);
    return _enclosedParticleCount;
  }

  /**
   * Check if the cell is not empty.
   * @return true if at least one particle is stored in this cell
   */
  [[nodiscard]] bool isEmpty() const override {
    std::lock_guard<AutoPasLock> lock(_lock);
    return _enclosedParticleCount == 0;
  }

  /**
   * Deletes all particles in this cell.
   */
  void clear() override {
    std::lock_guard<AutoPasLock> lock(_lock);
    _pointer->clearChildren(_pointer);
    _enclosedParticleCount = 0;
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
   * Delete the given particle from the data structure.
   * This function does not change the tree layout if the node is empty after the operation.
   * @param particle
   * @return True if the given pointer still points to a new, valid particle.
   */
  bool deleteParticle(Particle &particle) {
    --_enclosedParticleCount;
    return _pointer->deleteParticle(particle);
  };

  /**
   * Deletes the index-th particle.
   * @param index the index of the particle that shall be deleted
   */
  void deleteByIndex(size_t index) override {
    throw std::runtime_error("[OctreeNodeWrapper::deleteByIndex()] Operation not supported");
  }

  /**
   * Set the side lengths of this cell.
   * @param cellLength cell side length
   */
  void setCellLength(std::array<double, 3> &cellLength) override {
    throw std::runtime_error("[OctreeNodeWrapper::setCellLength()] Operation not supported");
  }

  /**
   * Get the side lengths of this cell.
   * @return cell side length
   */
  [[nodiscard]] std::array<double, 3> getCellLength() const override {
    return utils::ArrayMath::sub(_pointer->getBoxMin(), _pointer->getBoxMax());
  }

  /**
   * Get a particle from the iterator
   *
   * @param index The index of the particle
   * @return A ref to a particle
   */
  Particle &at(size_t index) {
    std::lock_guard<AutoPasLock> lock(_lock);
    return _ps.at(index);
  }

  /**
   * Get a particle from the iterator
   *
   * @param index The index of the particle
   * @return A const ref to a particle
   */
  const Particle &at(size_t index) const {
    std::lock_guard<AutoPasLock> lock(_lock);
    return _ps.at(index);
  }

  /**
   * Get a particle from the iterator
   *
   * @param index The index of the particle
   * @return A ref to a particle
   */
  Particle &operator[](size_t index) {
    std::lock_guard<AutoPasLock> lock(_lock);
    return *_ps[index];
  }

  /**
   * Get a particle from the iterator
   *
   * @param index The index of the particle
   * @return A ref to a particle
   */
  const Particle &operator[](size_t index) const {
    std::lock_guard<AutoPasLock> lock(_lock);
    return *_ps[index];
  }

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
   * Get a raw pointer to the enclosed cell.
   * @note This should only be used for debugging and insight into the internal structure fo the octree.
   * @return A raw C pointer to the enclosed node
   */
  OctreeNodeInterface<Particle> *getRaw() const { return _pointer.get(); }

  /**
   * Apply the forEach lambda to each particle.
   *
   * @tparam Lambda Function type
   * @param forEachLambda Function to apply
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda) {
    withStaticNodeType(_pointer, [&](auto nodePtr) { nodePtr->forEach(forEachLambda); });
  }

  /**
   * Apply the reduce lambda to each particle.
   *
   * @tparam Lambda Function type
   * @tparam A Initial value type
   * @param reduceLambda Function to apply
   * @param result Initial value
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result) {
    withStaticNodeType(_pointer, [&](auto nodePtr) { nodePtr->reduce(reduceLambda, result); });
  }

  /**
   * Apply the forEach lambda to each particle in the region.
   *
   * @tparam Lambda Function type
   * @param forEachLambda Function to apply
   * @param lowerCorner Lower corner of region
   * @param higherCorner Higher corner of region
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner) {
    withStaticNodeType(_pointer, [&](auto nodePtr) {
      // The iterator behavior is set to ownedOrHalo to include all particles inside this subtree. The baseclass
      // (Octree) decides whether this instance of OctreeNodeWrapper should be included in the iteration or not.
      nodePtr->forEach(forEachLambda, lowerCorner, higherCorner, IteratorBehavior::ownedOrHalo);
    });
  }

  /**
   * Apply the reduce lambda to each particle in the region.
   *
   * @tparam Lambda Function type
   * @tparam A Initial value type
   * @param reduceLambda Function to apply
   * @param result Initial value
   * @param lowerCorner Lower corner of region
   * @param higherCorner Higher corner of region
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner) {
    withStaticNodeType(_pointer, [&](auto nodePtr) {
      // The iterator behavior is set to ownedOrHalo to include all particles inside this subtree. The baseclass
      // (Octree) decides whether this instance of OctreeNodeWrapper should be included in the iteration or not.
      nodePtr->reduce(reduceLambda, result, lowerCorner, higherCorner, IteratorBehavior::ownedOrHalo);
    });
  }

 private:
  /**
   * A pointer to the root node of the enclosed octree.
   */
  std::unique_ptr<OctreeNodeInterface<Particle>> _pointer;

  /**
   * A vector containing all particles when iterating. This serves as a cache such that the octree does not need to
   * traversed every time the `at` function is called.
   * The field is marked mutable since it is used inside the `begin()` method, which is marked `const`. This function
   * loads the cache `ps`, which is therefore required to be mutable.
   */
  mutable StorageType _ps;

  /**
   * This lock is required to synchronize loading the cache `ps` across multiple threads. To change the state of the
   * lock from within the `begin()` method (marked `const`), this field has to be mutable.
   */
  mutable AutoPasLock _lock;

  /**
   * To prevent unnecessary tree-traversals of the octree, this field stores the number of enclosed particles within
   * this node.
   */
  long _enclosedParticleCount{0L};
};
}  // namespace autopas