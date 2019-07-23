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
   * @param interactionCellRadius radius of cells which can be in the interaction length
   * @param distance Maximum distance between interacting cells
   */
  virtual void rebuildVerlet(const std::array<unsigned long, 3> &dims, std::vector<ParticleCell> &cells,
                             std::vector<std::vector<std::array<double, 6>>> &boundingBoxes, int interactionCellRadius,
                             double distance) = 0;

  /**
   * Sets pointer to the verlet lists stored in the container
   * @param clusterSize pointer to size of clusters in container
   * @param neighborCellIds pointer to neighbor lists in container
   * @param neighborMatrixDim pointer to cuda neighbor matrix dimension
   * @param neighborMatrix pointer to cuda neighbor matrix dimension
   */
  virtual void setVerletListPointer(unsigned int *clusterSize,
                                    std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>> *neighborCellIds,
                                    size_t *neighborMatrixDim,
                                    utils::CudaDeviceVector<unsigned int> *neighborMatrix) = 0;

  /**
   * This function interacts all cells with the other cells with their index in neighborCellIds
   * @param cells containing the particles
   */
  virtual void traverseCellPairs(std::vector<ParticleCell> &cells) = 0;

  /**
   * This function returns the Data Layout Option and use of newton3 to identify a Traversal object
   * @return pair with DataLayoutOption and use of newton3
   */
  virtual std::tuple<TraversalOption, DataLayoutOption, bool> getSignature() = 0;
};

}  // namespace autopas
