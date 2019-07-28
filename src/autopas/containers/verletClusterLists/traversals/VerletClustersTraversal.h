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
 * Traversal for VerletClusterLists. Does not support newton 3.
 * @tparam ParticleCell
 * @tparam PairwiseFunctor The type of the functor.
 * @tparam dataLayout The data layout to use. Currently, only AoS is supported.
 * @tparam useNewton3 If newton 3 should be used. Currently, only false is supported.
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
class VerletClustersTraversal : public TraversalInterface,
                                public VerletClustersTraversalInterface<typename ParticleCell::ParticleType> {
  using Particle = typename ParticleCell::ParticleType;
  static constexpr size_t clusterSize = VerletClusterLists<Particle>::clusterSize;

 public:
  /**
   * Constructor of the VerletClustersTraversal.
   * @param pairwiseFunctor The functor to use for the traveral.
   */
  explicit VerletClustersTraversal(PairwiseFunctor *pairwiseFunctor) : _functor(pairwiseFunctor) {}

  TraversalOption getTraversalType() const override { return TraversalOption::verletClusters; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }
  bool getUseNewton3() const override { return useNewton3; }

  bool isApplicable() const override {
    return (dataLayout == DataLayoutOption::aos || dataLayout == DataLayoutOption::soa) and not useNewton3;
  }

  void initTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;
    // TODO: This and endTraversal(): implementation is the same for all traversals, so pull it in interface
    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
    auto &towers = clusterList.getTowers();
    const auto numTowers = towers.size();

    // TODO: Parallelize
    for (size_t index = 0; index < numTowers; index++) {
      towers[index].loadSoA(_functor);
    }
  }

  void endTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
    auto &towers = clusterList.getTowers();
    const auto numTowers = towers.size();

    // TODO: Parallelize
    for (size_t index = 0; index < numTowers; index++) {
      towers[index].extractSoA(_functor);
    }
  }

  void traverseParticlePairs() override {
    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;

    const auto _clusterTraverseFunctor = [this](internal::Cluster<Particle, clusterSize> &cluster) {
      traverseSingleCluster(cluster);
      for (auto *neighborCluster : cluster.getNeighbors()) {
        traverseNeighborClusters(cluster, *neighborCluster);
      }
    };

    clusterList.template traverseClusters<true>(_clusterTraverseFunctor);
  }

 private:
  void traverseSingleCluster(internal::Cluster<Particle, clusterSize> &cluster) {
    switch (dataLayout) {
      case DataLayoutOption::aos:
        traverseSingleClusterAoS(cluster);
        break;
      case DataLayoutOption::soa:
        traverseSingleClusterSoA(cluster);
        break;
      default:
        autopas::utils::ExceptionHandler::exception(
            "Wrong data layout of VerletClustersTraversal. Only AoS and SoA are supported!");
    }
  }

  void traverseSingleClusterAoS(internal::Cluster<Particle, clusterSize> &cluster) {
    for (size_t i = 0; i < clusterSize; i++) {
      for (size_t j = i + 1; j < clusterSize; j++) {
        auto &iParticle = cluster.getParticle(i);
        auto &jParticle = cluster.getParticle(j);
        _functor->AoSFunctor(iParticle, jParticle, useNewton3);
        if (not useNewton3) _functor->AoSFunctor(jParticle, iParticle, useNewton3);
      }
    }
  }

  void traverseSingleClusterSoA(internal::Cluster<Particle, clusterSize> &cluster) {
    _functor->SoAFunctor(cluster.getSoAView(), useNewton3);
  }

  void traverseNeighborClusters(internal::Cluster<Particle, clusterSize> &firstCluster,
                                internal::Cluster<Particle, clusterSize> &secondCluster) {
    switch (dataLayout) {
      case DataLayoutOption::aos:
        traverseNeighborClustersAoS(firstCluster, secondCluster);
        break;
      case DataLayoutOption::soa:
        traverseNeighborClustersSoA(firstCluster, secondCluster);
        break;
      default:
        autopas::utils::ExceptionHandler::exception(
            "Wrong data layout of VerletClustersTraversal. Only AoS and SoA are supported!");
    }
  }

  void traverseNeighborClustersAoS(internal::Cluster<Particle, clusterSize> &firstCluster,
                                   internal::Cluster<Particle, clusterSize> &secondCluster) {
    for (size_t i = 0; i < clusterSize; i++) {
      for (size_t j = 0; j < clusterSize; j++) {
        auto &iParticle = firstCluster.getParticle(i);
        auto &jParticle = secondCluster.getParticle(j);
        _functor->AoSFunctor(iParticle, jParticle, useNewton3);
      }
    }
  }

  void traverseNeighborClustersSoA(internal::Cluster<Particle, clusterSize> &firstCluster,
                                   internal::Cluster<Particle, clusterSize> &secondCluster) {
    _functor->SoAFunctor(firstCluster.getSoAView(), secondCluster.getSoAView(), useNewton3);
  }

 private:
  PairwiseFunctor *_functor;
};
}  // namespace autopas
