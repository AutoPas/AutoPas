/**
 * @file VerletClustersTraversal.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

#include "autopas/containers/TraversalInterface.h"
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
class VerletClustersTraversal : public TraversalInterface,
                                public VerletClustersTraversalInterface<typename ParticleCell::ParticleType> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor of the VerletClustersTraversal.
   * @param pairwiseFunctor The functor to use for the traveral.
   */
  explicit VerletClustersTraversal(PairwiseFunctor *pairwiseFunctor) : _functor(pairwiseFunctor) {}

  TraversalOption getTraversalType() const override { return TraversalOption::verletClusters; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }
  bool getUseNewton3() const override { return useNewton3; }
  bool isApplicable() const override { return (dataLayout == DataLayoutOption::aos) and not useNewton3; }

  void initTraversal() override {}
  void endTraversal() override {}

  void traverseParticlePairs() override {
    VerletClusterLists<Particle> &verletClusterLists = *(this->_verletClusterLists);

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
