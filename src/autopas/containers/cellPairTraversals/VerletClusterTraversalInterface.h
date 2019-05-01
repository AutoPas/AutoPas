/**
 * @file VerletClusterTraversalInterface.h
 * @author jspahl
 * @date 25.3.19
 */

#pragma once

#include <tuple>
#include <vector>

namespace autopas {

/**
 * This Traversal is used to interact all clusters in VerletClusterCluster Container
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell>
class VerletClusterTraversalInterface {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * rebuilds the Traversal and creates the neighbor lists
   * @param dims dimension of the container
   * @param cells in the container
   * @param boundingBoxes of the cells
   * @param distance Maximum distance between interacting cells
   */
  virtual void rebuild(const std::array<unsigned long, 3> &dims, unsigned int clusterSize,
                       std::vector<ParticleCell> &cells,
                       std::vector<std::array<typename Particle::ParticleFloatingPointType, 6>> &boundingBoxes,
                       typename Particle::ParticleFloatingPointType distance) = 0;

  /**
   * This function interacts all cells with the other cells with their index in neighborCellIds
   * @param cells containing the particles
   * @param neighborCellIds Stores the neighbor ids for each cell in cells
   */
  virtual void traverseCellPairs(std::vector<ParticleCell> &cells) = 0;

  /**
   * This function returns the Data Layout Option and use of newton3 to identify a Traversal object
   * @return pair with DataLayoutOption and use of newton3
   */
  virtual std::tuple<TraversalOption, DataLayoutOption, bool> getSignature() = 0;
};

}  // namespace autopas
