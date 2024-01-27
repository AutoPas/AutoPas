/**
 * @file VLTraversalInterface.h
 * @date 14.06.2019
 * @author C. Menges
 */

#pragma once

#include "autopas/containers/cellTraversals/CellTraversal.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletListHelpers.h"
#include "autopas/options/DataLayoutOption.h"

namespace autopas {

/**
 * This class provides the Traversal Interface for the verlet lists container.
 *
 * The container only accepts traversals in its computeInteractions() method that implement this interface.
 * @tparam LinkedParticleCell the type of cells
 */
template <class LinkedParticleCell>
class VLTraversalInterface {
 public:
  /**
   * Destructor
   */
  virtual ~VLTraversalInterface() = default;

  /**
   * Sets the information the traversal needs for the iteration.
   * @param cells The cells of the underlying LinkedCells container.
   * @param aosNeighborLists The AoS neighbor list.
   * @param soaNeighborLists The SoA neighbor list.
   */
  virtual void setCellsAndNeighborLists(
      std::vector<LinkedParticleCell> &cells,
      typename VerletListHelpers<typename LinkedParticleCell::ParticleType>::NeighborListAoSType &aosNeighborLists,
      std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> &soaNeighborLists,
      typename VerletListHelpers<typename LinkedParticleCell::ParticleType>::PairwiseNeighborListAoSType &pairwiseAosNeighborLists) {
    _cells = &cells;
    _aosNeighborLists = &aosNeighborLists;
    _soaNeighborLists = &soaNeighborLists;
    _pairwiseAosNeighborLists = &pairwiseAosNeighborLists;
  }

 protected:
  /**
   * The cells of the underlying linked cells container of the verlet lists container.
   */
  std::vector<LinkedParticleCell> *_cells = nullptr;
  /**
   * The AoS neighbor list of the verlet lists container.
   */
  typename VerletListHelpers<typename LinkedParticleCell::ParticleType>::NeighborListAoSType *_aosNeighborLists =
      nullptr;
  /**
   * The SoA neighbor list of the verlet lists container.
   */
  std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> *_soaNeighborLists = nullptr;
  /**
   * The AoS neighbor list of the verlet lists container.
   */
  typename VerletListHelpers<typename LinkedParticleCell::ParticleType>::PairwiseNeighborListAoSType *_pairwiseAosNeighborLists =
      nullptr;
};

}  // namespace autopas