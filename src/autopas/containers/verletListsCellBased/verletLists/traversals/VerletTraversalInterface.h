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
   * Constructor
   */
  VerletTraversalInterface() = default;

  /**
   * Destructor
   */
  virtual ~VerletTraversalInterface() = default;

  /**
   * Iterates over the Particles as specified in the Neighbor lists
   * @param aosNeighborLists neighbor lists in aos format
   * @param soaNeighborLists neighbor lists as index list for the soa format
   */
  virtual void iterateVerletLists(
      std::unordered_map<typename LinkedParticleCell::ParticleType *,
                         std::vector<typename LinkedParticleCell::ParticleType *>>
          aosNeighborLists,
      std::vector<std::vector<size_t, autopas::AlignedAllocator<size_t>>> soaNeighborLists) = 0;

  /**
   * Initializes Traversal and copies all relevant data
   * @param cells content of the container the Traversal is to be called on
   */
  virtual void initTraversal(std::vector<LinkedParticleCell> &cells) = 0;

  /**
   * Ends Traversal write back data
   * @param cells content of the container the Traversal is to be called on
   */
  virtual void endTraversal(std::vector<LinkedParticleCell> &cells) = 0;

  /**
   * Returns data layout.
   * @return dataLayout
   */
  virtual autopas::DataLayoutOption getDataLayout() = 0;
};

}  // namespace autopas