/**
 * @file VerletClustersColoringTraversal.h
 * @author humig
 * @date 27.06.19
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CBasedTraversal.h"
#include "autopas/containers/verletClusterLists/traversals/VerletClustersTraversalInterface.h"

namespace autopas {

/**
 * A traversal for VerletClusterLists that uses a coloring over the grids of the container.
 *
 * The traversal uses a 2D coloring with a stride of x=3, y=2, so 3*2=6 colors.
 *
 * When disabling newton 3, interactions inside a cluster are still calculated using newton 3.
 *
 * @tparam ParticleCell
 * @tparam PairwiseFunctor
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
class VerletClustersColoringTraversal : public CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                                        public VerletClustersTraversalInterface<typename ParticleCell::ParticleType> {
 private:
  using Particle = typename ParticleCell::ParticleType;
  static constexpr size_t clusterSize = VerletClusterLists<Particle>::clusterSize;

  /**
   * Each base step looks like this:
   *    X C N  Colors:  1 2 3
   *    N N N           4 5 6
   * Where C is the current cell, N are the neighbor cells that is worked on, and X is not worked on. The neighbor list
   * with newton 3 of the VerletClusterLists container is build in a way that the neighbor lists already contain only
   * the neighbor clusters of these cells.
   * @see VerletClusterLists::updateVerletLists(bool)
   */
  static constexpr std::array<unsigned long, 3> _stride{3ul, 2ul, 1ul};

  /**
   * Helper method to iterate over one color cell.
   * @param xColorCell The x coordinate of the cell.
   * @param yColorCell The y coordinate of the cell.
   * @param zColorCell The z coordinate of the cell.
   * @param towersPerColoringCell The number of grids that every cell has in every dimension.
   */
  void processColorCell(unsigned long xColorCell, unsigned long yColorCell, unsigned long zColorCell,
                        int towersPerColoringCell);

  /**
   * Helper method to traverse two neighbor clusters.
   * @param clusterStart The first cluster.
   * @param neighborClusterStart The second cluster.
   * @param clusterSize The size of the cluster.
   */
  void traverseClusterPairAoS(internal::Cluster<Particle, clusterSize> &cluster,
                              internal::Cluster<Particle, clusterSize> &neighborCluster);

  void traverseClusterPairSoA(internal::Cluster<Particle, clusterSize> &cluster,
                              internal::Cluster<Particle, clusterSize> &neighborCluster);

  void traverseClusterSoA(internal::Cluster<Particle, clusterSize> &cluster);

  void traverseClusterAoS(internal::Cluster<Particle, clusterSize> &cluster);

 public:
  /**
   * Constructor of the VerletClustersTraversal.
   * @param pairwiseFunctor The functor to use for the traveral.
   */
  explicit VerletClustersColoringTraversal(PairwiseFunctor *pairwiseFunctor)
      : CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>({0, 0, 0}, pairwiseFunctor, 0, {}),
        _functor(pairwiseFunctor) {}

  TraversalOption getTraversalType() const override { return TraversalOption::verletClustersColoring; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }
  bool getUseNewton3() const override { return useNewton3; }
  bool isApplicable() const override {
    return (dataLayout == DataLayoutOption::aos || dataLayout == DataLayoutOption::soa);
  }

  void initTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
    clusterList.loadParticlesIntoSoAs(_functor);
  }

  void endTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
    clusterList.extractParticlesFromSoAs(_functor);
  }

  void traverseParticlePairs() override {
    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;

    const auto towersPerColoringCell = clusterList.getInteractionLengthInTowers();
    std::array<unsigned long, 3> coloringCellsPerDim{};
    for (int i = 0; i < 3; i++) {
      coloringCellsPerDim[i] =
          static_cast<unsigned long>(std::ceil(clusterList.getTowersPerDimension()[i] / (double)towersPerColoringCell));
    }

    auto loopBody = [this, towersPerColoringCell](unsigned long x, unsigned long y, unsigned long z) {
      processColorCell(x, y, z, towersPerColoringCell);
    };

    // We are only doing a 2D coloring.
    if (coloringCellsPerDim[2] != 1) {
      autopas::utils::ExceptionHandler::exception(
          "VerletClusterColoringTraversal: Coloring should only be 2D, not in z-direction!");
    }

    // localStride is necessary because stride is constexpr and cTraversal() wants a const &
    auto localStride = _stride;
    this->cTraversal(std::forward<decltype(loopBody)>(loopBody), coloringCellsPerDim, localStride);
  }

 private:
  PairwiseFunctor *_functor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processColorCell(
    unsigned long xColorCell, unsigned long yColorCell, unsigned long zColorCell, int towersPerColoringCell) {
  // We are only doing a 2D coloring.
  if (zColorCell != 0) {
    autopas::utils::ExceptionHandler::exception("Coloring should only be 2D, not in z-direction!");
  }

  auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
  auto &towers = clusterList.getTowers();
  const auto cellsPerDim = clusterList.getTowersPerDimension();

  for (int yInner = 0; yInner < towersPerColoringCell; yInner++) {
    for (int xInner = 0; xInner < towersPerColoringCell; xInner++) {
      const auto y = yColorCell * towersPerColoringCell + yInner;
      const auto x = xColorCell * towersPerColoringCell + xInner;

      // Not every coloring cell has to have gridsPerColoringCell grids in every direction.
      if (x >= cellsPerDim[0] or y >= cellsPerDim[1]) {
        continue;
      }
      auto gridIndex1D = VerletClusterMaths::index1D(x, y, cellsPerDim);

      auto &currentTower = towers[gridIndex1D];
      for (auto &cluster : currentTower.getClusters()) {
        if constexpr (dataLayout == DataLayoutOption::aos) {
          traverseClusterAoS(cluster);
        } else {
          traverseClusterSoA(cluster);
        }

        for (auto *neighborCluster : cluster.getNeighbors()) {
          if constexpr (dataLayout == DataLayoutOption::aos) {
            traverseClusterPairAoS(cluster, *neighborCluster);
          } else {
            traverseClusterPairSoA(cluster, *neighborCluster);
          }
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseClusterAoS(
    internal::Cluster<Particle, clusterSize> &cluster) {
  for (size_t i = 0; i < clusterSize; i++) {
    // Always use newton 3 for interactions within one cluster.
    for (size_t j = i + 1; j < clusterSize; j++) {
      _functor->AoSFunctor(cluster.getParticle(i), cluster.getParticle(j), true);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseClusterSoA(
    internal::Cluster<Particle, clusterSize> &cluster) {
  _functor->SoAFunctor(cluster.getSoAView(), useNewton3);
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseClusterPairAoS(
    internal::Cluster<Particle, clusterSize> &cluster, internal::Cluster<Particle, clusterSize> &neighborCluster) {
  for (size_t i = 0; i < clusterSize; i++) {
    for (size_t j = 0; j < clusterSize; j++) {
      _functor->AoSFunctor(cluster.getParticle(i), neighborCluster.getParticle(j), useNewton3);
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseClusterPairSoA(
    internal::Cluster<Particle, clusterSize> &cluster, internal::Cluster<Particle, clusterSize> &neighborCluster) {
  _functor->SoAFunctor(cluster.getSoAView(), neighborCluster.getSoAView(), useNewton3);
}

}  // namespace autopas
