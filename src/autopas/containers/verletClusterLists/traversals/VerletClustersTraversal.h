/**
 * @file VerletClustersTraversal.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClustersTraversalInterface.h"

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
  void traverseParticlePairs(VerletClusterLists<Particle> &verletClusterLists) override {
    const auto _clusterTraverseFunctor = [functor = _functor](Particle *clusterStart, int clusterSize,
                                                              std::vector<Particle *> &clusterNeighborList) {
      for (auto neighbor : clusterNeighborList) {
        if (clusterStart == neighbor) {
          // self pair
          for (int i = 0; i < clusterSize; i++) {
            for (int j = i + 1; j < clusterSize; j++) {
              Particle *iParticle = clusterStart + i;
              Particle *jParticle = neighbor + j;
              functor->AoSFunctor(*iParticle, *jParticle, useNewton3);
              if (not useNewton3) functor->AoSFunctor(*jParticle, *iParticle, useNewton3);
            }
          }
        } else {
          for (int i = 0; i < clusterSize; i++) {
            for (int j = 0; j < clusterSize; j++) {
              Particle *iParticle = clusterStart + i;
              Particle *jParticle = neighbor + j;
              functor->AoSFunctor(*iParticle, *jParticle, useNewton3);
            }
          }
        }
      }
    };

    verletClusterLists.template traverseClusters<true>(_clusterTraverseFunctor);
  }

 private:
  PairwiseFunctor *_functor;
};
}  // namespace autopas
