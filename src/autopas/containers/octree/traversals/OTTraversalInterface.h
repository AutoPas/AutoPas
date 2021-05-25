/**
 * @file OTTraversalInterface.h
 *
 * @author Johannes Spies
 * @date 11.05.2021
 */

#pragma once

#include "autopas/containers/octree/OctreeNodeWrapper.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/DataLayoutConverter.h"

namespace autopas {

/**
 * This interface exists to provide a row interface for octree to add its cells.
 */
template <typename ParticleCell>
class OTTraversalInterface {
 public:
  /**
   * Notify the traversal about the cells that it is able to traverse.
   * @param cells A vector of size 2 containing the owned and the halo octrees.
   */
  virtual void setCells(std::vector<ParticleCell> *cells) = 0;
};

} // namespace autopas
