/**
 * @file VLTraversalInterface.h
 * @date 14.06.2019
 * @author C. Menges
 */

#pragma once

#include <unordered_map>

#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"

namespace autopas {

/**
 * This class provides the Traversal Interface for the verlet lists container.
 *
 * The container only accepts traversals in its computeInteractions() method that implement this interface.
 * @tparam LinkedParticleCell the type of cells
 */
template <class LinkedParticleCell>
class VLTraversalInterface {
  using ParticleType = typename LinkedParticleCell::ParticleType;

 public:
  /**
   * Destructor
   */
  virtual ~VLTraversalInterface() = default;

  /**
   * Sets the information the traversal needs for the iteration.
   * @param cells         The cells of the underlying LinkedCells container.
   * @param neighborList  The flat CRS neighbor list.
   * @param particleToIndex Map from particle pointer to its dense SoA index.
   */
  virtual void setCellsAndNeighborLists(std::vector<LinkedParticleCell> &cells,
                                        typename VerletListHelpers<ParticleType>::NeighborListCRS &neighborList,
                                        const std::unordered_map<const ParticleType *, size_t> &particleToIndex) {
    _cells = &cells;
    _neighborList = &neighborList;
    _particleToIndex = &particleToIndex;
  }

 protected:
  /**
   * The cells of the underlying linked cells container of the verlet lists container.
   */
  std::vector<LinkedParticleCell> *_cells = nullptr;

  /**
   * The flat CRS neighbor list.
   */
  typename VerletListHelpers<ParticleType>::NeighborListCRS *_neighborList = nullptr;

  /**
   * Map from particle pointer to its dense SoA index.
   * Used by the AoS traversal path to resolve indices back to particle references.
   */
  const std::unordered_map<const ParticleType *, size_t> *_particleToIndex = nullptr;
};

}  // namespace autopas