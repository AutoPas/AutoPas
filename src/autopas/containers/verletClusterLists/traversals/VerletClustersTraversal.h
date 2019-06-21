/**
 * @file VerletClustersTraversal.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

#include "VerletClustersTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"

namespace autopas {

/**
 * Traversal for VerletClusterLists.
 * @tparam ParticleCell
 * @tparam PairwiseFunctor The type of the functor.
 * @tparam dataLayout The data layout to use. Currently, only AoS is supported.
 * @tparam useNewton3 If newton 3 should be used. Currently, only false is supported.
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
class VerletClustersTraversal : public CellPairTraversal<ParticleCell, dataLayout, useNewton3>,
                                public VerletClustersTraversalInterface<typename ParticleCell::ParticleType> {
  using Particle = typename ParticleCell::ParticleType;
  using index_t = typename VerletClusterMaths::index_t;

 public:
  /**
   * Constructor of the VerletClustersTraversal.
   * @param pairwiseFunctor The functor to use for the traveral.
   */
  explicit VerletClustersTraversal(PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell, dataLayout, useNewton3>({0, 0, 0}), _functor(pairwiseFunctor) {}

  TraversalOption getTraversalType() const override { return TraversalOption::verletClusters; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }
  bool getUseNewton3() const override { return useNewton3; }
  bool isApplicable() const override { return (dataLayout == DataLayoutOption::aos) && not useNewton3; }

  void initTraversal(std::vector<ParticleCell> &cells) override {}
  void endTraversal(std::vector<ParticleCell> &cells) override {}

  /**
   * @copydoc VerletClustersTraversalInterface::traverseParticlePairs
   */
  void traverseParticlePairs(std::array<index_t, 3> cellsPerDim, int clusterSize,
                             std::vector<FullParticleCell<Particle>> &clusters,
                             std::vector<std::vector<std::vector<Particle *>>> &neighborLists) override {
    const index_t end_x = cellsPerDim[0];
    const index_t end_y = cellsPerDim[1];

#if defined(AUTOPAS_OPENMP)
    // @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    for (index_t x = 0; x < end_x; x++) {
      for (index_t y = 0; y < end_y; y++) {
        index_t index = VerletClusterMaths::index1D(x, y, cellsPerDim);
        auto &grid = clusters[index];
        auto &gridVerlet = neighborLists[index];

        const index_t gridSize = grid.numParticles() / clusterSize;
        for (index_t z = 0; z < gridSize; z++) {
          Particle *iClusterStart = &grid[z * clusterSize];
          for (auto neighbor : gridVerlet[z]) {
            if (iClusterStart == neighbor) {
              // self pair
              for (int i = 0; i < clusterSize; i++) {
                for (int j = i + 1; j < clusterSize; j++) {
                  Particle *iParticle = iClusterStart + i;
                  Particle *jParticle = neighbor + j;
                  _functor->AoSFunctor(*iParticle, *jParticle, useNewton3);
                  if (not useNewton3) _functor->AoSFunctor(*jParticle, *iParticle, useNewton3);
                }
              }
            } else {
              for (int i = 0; i < clusterSize; i++) {
                for (int j = 0; j < clusterSize; j++) {
                  Particle *iParticle = iClusterStart + i;
                  Particle *jParticle = neighbor + j;
                  _functor->AoSFunctor(*iParticle, *jParticle, useNewton3);
                }
              }
            }
          }
        }
      }
    }
  }

 private:
  PairwiseFunctor *_functor;
};
}  // namespace autopas
