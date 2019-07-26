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

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;

    auto numClusters = clusterList.getNumClusters();
    const auto &aosToSoaMap = clusterList.getClusterIndexMap();

    _clusterSoAs.resize(numClusters);

    const auto _clusterTraverseFunctor = [this, &aosToSoaMap](Particle *clusterStart,
                                                              std::vector<Particle *> &clusterNeighborList) {
      constexpr auto clusterSize = VerletClusterLists<Particle>::clusterSize;
      auto currentClusterIndex = aosToSoaMap.at(clusterStart);
      FullParticleCell<Particle> cell{};
      cell.reserve(clusterSize);
      for (size_t i = 0; i < clusterSize; i++) {
        cell.addParticle(*(clusterStart + i));
      }
      SoA<typename Particle::SoAArraysType> &soa = _clusterSoAs[currentClusterIndex];
      soa.resizeArrays(clusterSize);
      _functor->SoALoader(cell, soa);
    };

    clusterList.template traverseClusters<true>(_clusterTraverseFunctor);
  }

  void endTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;

    const auto &aosToSoaMap = clusterList.getClusterIndexMap();

    const auto _clusterTraverseFunctor = [this, &aosToSoaMap](Particle *clusterStart,
                                                              std::vector<Particle *> &clusterNeighborList) {
      constexpr auto clusterSize = VerletClusterLists<Particle>::clusterSize;
      auto currentClusterIndex = aosToSoaMap.at(clusterStart);
      FullParticleCell<Particle> cell{};
      cell.reserve(clusterSize);
      for (size_t i = 0; i < clusterSize; i++) {
        cell.addParticle(*(clusterStart + i));
      }
      SoA<typename Particle::SoAArraysType> &soa = _clusterSoAs[currentClusterIndex];
      _functor->SoAExtractor(cell, soa);
      for (size_t i = 0; i < clusterSize; i++) {
        *(clusterStart + i) = cell[i];
      }
    };

    clusterList.template traverseClusters<true>(_clusterTraverseFunctor);
  }

  void traverseParticlePairs() override {
    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;

    const auto &aosToSoaMap = clusterList.getClusterIndexMap();

    const auto _clusterTraverseFunctor = [this, &aosToSoaMap](Particle *clusterStart,
                                                              std::vector<Particle *> &clusterNeighborList) {
      for (auto neighborClusterStart : clusterNeighborList) {
        // self pair
        if (clusterStart == neighborClusterStart) {
          traverseSingleCluster(clusterStart, aosToSoaMap);
        } else {
          traverseNeighborClusters(clusterStart, neighborClusterStart, aosToSoaMap);
        }
      }
    };

    clusterList.template traverseClusters<true>(_clusterTraverseFunctor);
  }

 private:
  void traverseSingleCluster(Particle *clusterStart,
                             const std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap) {
    switch (dataLayout) {
      case DataLayoutOption::aos:
        traverseSingleClusterAoS(clusterStart);
        break;
      case DataLayoutOption::soa:
        traverseSingleClusterSoA(clusterStart, aosToSoaMap);
        break;
      default:
        autopas::utils::ExceptionHandler::exception(
            "Wrong data layout of VerletClustersTraversal. Only AoS and SoA are supported!");
    }
  }

  void traverseSingleClusterAoS(Particle *clusterStart) {
    constexpr auto clusterSize = VerletClusterLists<Particle>::clusterSize;
    for (size_t i = 0; i < clusterSize; i++) {
      for (size_t j = i + 1; j < clusterSize; j++) {
        Particle *iParticle = clusterStart + i;
        Particle *jParticle = clusterStart + j;
        _functor->AoSFunctor(*iParticle, *jParticle, useNewton3);
        if (not useNewton3) _functor->AoSFunctor(*jParticle, *iParticle, useNewton3);
      }
    }
  }

  void traverseSingleClusterSoA(Particle *clusterStart,
                                const std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap) {
    _functor->SoAFunctor(_clusterSoAs[aosToSoaMap.at(clusterStart)], useNewton3);
  }

  void traverseNeighborClusters(Particle *firstClusterStart, Particle *secondClusterStart,
                                const std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap) {
    switch (dataLayout) {
      case DataLayoutOption::aos:
        traverseNeighborClustersAoS(firstClusterStart, secondClusterStart);
        break;
      case DataLayoutOption::soa:
        traverseNeighborClustersSoA(firstClusterStart, secondClusterStart, aosToSoaMap);
        break;
      default:
        autopas::utils::ExceptionHandler::exception(
            "Wrong data layout of VerletClustersTraversal. Only AoS and SoA are supported!");
    }
  }

  void traverseNeighborClustersAoS(Particle *firstClusterStart, Particle *secondClusterStart) {
    constexpr auto clusterSize = VerletClusterLists<Particle>::clusterSize;
    for (size_t i = 0; i < clusterSize; i++) {
      for (size_t j = 0; j < clusterSize; j++) {
        Particle *iParticle = firstClusterStart + i;
        Particle *jParticle = secondClusterStart + j;
        _functor->AoSFunctor(*iParticle, *jParticle, useNewton3);
      }
    }
  }

  void traverseNeighborClustersSoA(Particle *firstClusterStart, Particle *secondClusterStart,
                                   const std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap) {
    _functor->SoAFunctor(_clusterSoAs[aosToSoaMap.at(firstClusterStart)],
                         _clusterSoAs[aosToSoaMap.at(secondClusterStart)], useNewton3);
  }

 private:
  PairwiseFunctor *_functor;

  std::vector<SoA<typename Particle::SoAArraysType>> _clusterSoAs;
};
}  // namespace autopas
