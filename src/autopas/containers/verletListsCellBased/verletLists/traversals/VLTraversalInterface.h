/**
 * @file VLTraversalInterface.h
 * @date 14.06.2019
 * @author C. Menges
 */

#pragma once

#include <unordered_map>

#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"

namespace autopas {

/**
 * This class provides the Traversal Interface for the verlet lists container.
 *
 * The container only accepts traversals in its computeInteractions() method that implements this interface.
 * @tparam LinkedParticleCell the type of cells
 */
template <class LinkedParticleCell>
class VLTraversalInterface {
  using ParticleType = LinkedParticleCell::ParticleType;

 public:
  /**
   * Destructor
   */
  virtual ~VLTraversalInterface() = default;

  /**
   * Sets the information the traversal needs for the iteration.
   * @param cells         The cells of the underlying LinkedCells container.
   * @param neighborList  The flat CRS neighbor list.
   * @param particleToIndex Map from particle pointers to their dense SoA index.
   * @param indexToParticle Vector mapping SoA indices to particle pointers.
   */
  virtual void setCellsAndNeighborLists(std::vector<LinkedParticleCell> &cells,
                                        VerletListHelpers<ParticleType>::NeighborListCRS &neighborList,
                                        const std::unordered_map<const ParticleType *, size_t> &particleToIndex,
                                        const std::vector<ParticleType *> &indexToParticle) {
    _cells = &cells;
    _neighborList = &neighborList;
    _particleToIndex = &particleToIndex;
    _indexToParticle = &indexToParticle;
  }

 protected:
  /**
   * The cells of the underlying linked cells container of the verlet lists container.
   */
  std::vector<LinkedParticleCell> *_cells = nullptr;

  /**
   * The flat CRS neighbor list.
   */
  VerletListHelpers<ParticleType>::NeighborListCRS *_neighborList = nullptr;

  /**
   * Map from particle pointers to their dense SoA index.
   * Used by the rebuild traversal to create index-based neighbor lists.
   */
  const std::unordered_map<const ParticleType *, size_t> *_particleToIndex = nullptr;

  /**
   * Vector mapping SoA indices to particle pointers.
   * Used by the AoS traversal path to resolve indices back to particle references.
   */
  const std::vector<ParticleType *> *_indexToParticle = nullptr;
};

}  // namespace autopas