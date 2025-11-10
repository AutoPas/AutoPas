#pragma once
#include "autopas/containers/ParticleContainerInterface.h"

namespace autopas {

template <typename Particle_T>
class IteratorParticleContainer : public ParticleContainerInterface<Particle_T> {
 public:
  virtual ~IteratorParticleContainer() = default;
  IteratorParticleContainer() = default;
  /**
   * Iterate over all particles using
   * for(auto iter = container.begin(); iter.isValid(); ++iter) .
   * @note The default argument for behavior is necessary to enable range based for loops.
   * @param behavior Behavior of the iterator, see IteratorBehavior.
   * @param additionalVectors Vectors that should be included besides the container.
   * @return Iterator to the first particle.
   */
  [[nodiscard]] virtual ContainerIterator<Particle_T, true, false> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, true, false>::ParticleVecType *additionalVectors = nullptr) = 0;

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   * @note const version
   */
  [[nodiscard]] virtual ContainerIterator<Particle_T, false, false> begin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, false, false>::ParticleVecType *additionalVectors = nullptr) const = 0;

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   * @note cbegin will guarantee to return a const_iterator.
   */
  [[nodiscard]] ContainerIterator<Particle_T, false, false> cbegin(
      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo,
      typename ContainerIterator<Particle_T, false, false>::ParticleVecType *additionalVectors = nullptr) const {
    return begin(behavior, additionalVectors);
  }

  /**
   * Iterate over all particles in a specified region
   * for(auto iter = container.getRegionIterator(lowCorner, highCorner);iter.isValid();++iter) .
   * @param lowerCorner Lower corner of the region
   * @param higherCorner Higher corner of the region
   * @param behavior The behavior of the iterator (shall it iterate over halo particles as well?).
   * @param additionalVectors Vectors that should be included besides the container.
   * @return Iterator to iterate over all particles in a specific region.
   */
  [[nodiscard]] virtual ContainerIterator<Particle_T, true, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle_T, true, true>::ParticleVecType *additionalVectors = nullptr) = 0;

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   * @note const version
   */
  [[nodiscard]] virtual ContainerIterator<Particle_T, false, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner, IteratorBehavior behavior,
      typename ContainerIterator<Particle_T, false, true>::ParticleVecType *additionalVectors = nullptr) const = 0;

  /**
   * @copydoc autopas::AutoPas::end()
   */
  [[nodiscard]] constexpr bool end() const { return false; }

  /**
   * Fetch the pointer to a particle, identified via a cell and particle index.
   * These indices are only meaningful in the context of the current container at its current state.
   * The same indices might (and probably will) yield different particles for different container types or might not
   * even exist.
   * The only guarantee is that the indices {0,0} yield the first particle in the container that satisfies the iterator
   * requirements.
   *
   * @note This function should handle any offsets if used in a parallel iterator.
   *
   * @param cellIndex Index of the cell the particle is located in.
   * @param particleIndex Particle index within the cell.
   * @param iteratorBehavior Which ownership states should be considered for the next particle.
   * @return Pointer to the particle and its indices. tuple<Particle_T*, cellIndex, particleIndex>
   * If a index pair is given that does not exist but is also not beyond the last cell, the next fitting particle shall
   * be returned.
   * Example: If [4,2] does not exist, [5,1] shall be returned
   * (or whatever is the next particle that fulfills the iterator requirements).
   * If there is no next fitting particle {nullptr, 0, 0} is returned.
   */
  virtual std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                                     IteratorBehavior iteratorBehavior) const = 0;

  /**
   * @copydoc getParticle(size_t cellIndex, size_t particleIndex, IteratorBehavior iteratorBehavior) const
   *
   * @param boxMin start of region in which the next particle should be. The coordinates are expected to be within the
   * domain.
   * @param boxMax end of region in which the next particle should be. The coordinates are expected to be within the
   * domain.
   */
  virtual std::tuple<const Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                                     IteratorBehavior iteratorBehavior,
                                                                     const std::array<double, 3> &boxMin,
                                                                     const std::array<double, 3> &boxMax) const = 0;

  // clang-format off
  /**
   * @copydoc getParticle(size_t cellIndex, size_t particleIndex, IteratorBehavior iteratorBehavior, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax) const
   *
   * @note non-const region iter version
   */
  // clang-format on
  std::tuple<Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                       IteratorBehavior iteratorBehavior,
                                                       const std::array<double, 3> &boxMin,
                                                       const std::array<double, 3> &boxMax) {
    const Particle_T *ptr{};
    size_t nextCellIndex{}, nextParticleIndex{};
    std::tie(ptr, nextCellIndex, nextParticleIndex) =
        const_cast<const IteratorParticleContainer<Particle_T> *>(this)->getParticle(cellIndex, particleIndex,
                                                                                     iteratorBehavior, boxMin, boxMax);
    return {const_cast<Particle_T *>(ptr), nextCellIndex, nextParticleIndex};
  }

  /**
   * @copydoc getParticle(size_t cellIndex, size_t particleIndex, IteratorBehavior iteratorBehavior) const
   *
   * @note non-const non-region iter version
   */
  std::tuple<Particle_T *, size_t, size_t> getParticle(size_t cellIndex, size_t particleIndex,
                                                       IteratorBehavior iteratorBehavior) {
    const Particle_T *ptr{};
    size_t nextCellIndex{}, nextParticleIndex{};
    std::tie(ptr, nextCellIndex, nextParticleIndex) =
        const_cast<const IteratorParticleContainer<Particle_T> *>(this)->getParticle(cellIndex, particleIndex,
                                                                                     iteratorBehavior);
    return {const_cast<Particle_T *>(ptr), nextCellIndex, nextParticleIndex};
  }

  /**
   * Deletes the particle at the given index positions as long as this does not compromise the validity of the
   * container. If this is not possible the particle is just marked as deleted. If the positions do not exist the
   * behavior is undefined.
   * @param cellIndex
   * @param particleIndex
   * @return True if the given indices still point to a new particle.
   */
  virtual bool deleteParticle(size_t cellIndex, size_t particleIndex) = 0;

};

}  // namespace autopas
