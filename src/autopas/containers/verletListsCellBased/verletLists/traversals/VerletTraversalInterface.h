/**
 * @file VerletTraversalInterface.h
 * @date 14.06.2019
 * @author C. Menges
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * This class provides the Traversal Interface for the verlet lists container.
 *
 * @tparam LinkedParticleCell the type of cells
 */
template <class LinkedParticleCell>
class VerletTraversalInterface {
 public:
  /**
   * Destructor
   */
  virtual ~VerletTraversalInterface() = default;

  /**
   * Sets the information the traversal needs for the iteration.
   * @param cells The cells of the underlying LinkedCells container.
   * @param aosNeighborLists The AoS neighbor list.
   * @param soaNeighborLists The SoA neighbor list.
   */
  virtual void setCellsAndNeighborLists(
      std::vector<LinkedParticleCell> &cells,
      std::unordered_map<typename LinkedParticleCell::ParticleType *,
                         std::vector<typename LinkedParticleCell::ParticleType *>> &aosNeighborLists,
      std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &soaNeighborLists) {
    _cells = &cells;
    _aosNeighborLists = &aosNeighborLists;
    _soaNeighborLists = &soaNeighborLists;
  }

 protected:
  std::vector<LinkedParticleCell> *_cells;
  std::unordered_map<typename LinkedParticleCell::ParticleType *,
                     std::vector<typename LinkedParticleCell::ParticleType *>> *_aosNeighborLists;
  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> *_soaNeighborLists;
};

}  // namespace autopas