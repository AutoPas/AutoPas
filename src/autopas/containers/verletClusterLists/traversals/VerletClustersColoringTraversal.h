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
 * It only supports AoS at the moment.
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
  using Particle = typename ParticleCell::ParticleType;
  typedef typename VerletClusterMaths::index_t index_t;
  /**
   * Each base step looks like this:
   *    X C N  Colors:  1 2 3
   *    N N N           4 5 6
   * Where C is the current cell, N are the neighbor cells that is worked on, and X is not worked on. The neighbor list
   * with newton 3 of the VerletClusterLists container is build in a way that the neighbor lists already contain only
   * the neighbor clusters of these cells.
   * @see VerletClusterLists::updateVerletLists(bool)
   */
  static constexpr std::array<unsigned long, 3> stride{3ul, 2ul, 1ul};

 private:
  /**
   * Helper method to iterate over one color cell.
   * @param xColorCell The x coordinate of the cell.
   * @param yColorCell The y coordinate of the cell.
   * @param zColorCell The z coordinate of the cell.
   * @param gridsPerColoringCell The number of grids that every cell has in every dimension.
   */
  void processColorCell(unsigned long xColorCell, unsigned long yColorCell, unsigned long zColorCell,
                        int gridsPerColoringCell);

  /**
   * Helper method to traverse two neighbor clusters.
   * @param clusterStart The first cluster.
   * @param neighborClusterStart The second cluster.
   * @param clusterSize The size of the cluster.
   */
  void traverseClusterPair(Particle *clusterStart, Particle *neighborClusterStart, int clusterSize);

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
  bool isApplicable() const override { return (dataLayout == DataLayoutOption::aos); }

  void initTraversal() override {}

  void endTraversal() override {}

  void traverseParticlePairs() override {
    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;

    double gridSideLength = clusterList.getGridSideLength();
    double cutoff = clusterList.getCutoff();
    double gridsPerColoringCell = std::ceil(cutoff / gridSideLength);
    std::array<unsigned long, 3> coloringCellsPerDim{};
    for (int i = 0; i < 3; i++) {
      coloringCellsPerDim[i] =
          static_cast<unsigned long>(std::ceil(clusterList.getCellsPerDimension()[i] / gridsPerColoringCell));
    }

    auto loopBody = [this, gridsPerColoringCell](unsigned long x, unsigned long y, unsigned long z) {
      processColorCell(x, y, z, gridsPerColoringCell);
    };

    const auto end = ArrayMath::sub(coloringCellsPerDim, {0ul, 0ul, 0ul});
    // We are only doing a 2D coloring.
    if (end[2] != 1) {
      autopas::utils::ExceptionHandler::exception("Coloring should only be 2D, not in z-direction!");
    }

    // localStride is necessary because stride is constexpr and cTraversal() wants a const &
    auto localStride = stride;
    this->cTraversal(std::forward<decltype(loopBody)>(loopBody), end, localStride);
  }

 private:
  PairwiseFunctor *_functor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processColorCell(
    unsigned long xColorCell, unsigned long yColorCell, unsigned long zColorCell, int gridsPerColoringCell) {
  // We are only doing a 2D coloring.
  if (zColorCell != 0) {
    autopas::utils::ExceptionHandler::exception("Coloring should only be 2D, not in z-direction!");
  }

  auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
  const auto cellsPerDim = clusterList.getCellsPerDimension();
  auto &grids = clusterList.getGrids();
  const auto &neighborLists = clusterList.getNeighborLists();
  auto clusterSize = clusterList.getClusterSize();

  for (int yInner = 0; yInner < gridsPerColoringCell; yInner++) {
    for (int xInner = 0; xInner < gridsPerColoringCell; xInner++) {
      unsigned long y = yColorCell * gridsPerColoringCell + yInner;
      unsigned long x = xColorCell * gridsPerColoringCell + xInner;

      // Not every coloring cell has to have gridsPerColoringCell grids in every direction.
      if (x >= cellsPerDim[0] || y >= cellsPerDim[1]) {
        continue;
      }
      unsigned long gridIndex1D = VerletClusterMaths::index1D(x, y, cellsPerDim);

      auto &currentGrid = grids[gridIndex1D];
      auto numClusters = currentGrid.numParticles() / clusterSize;
      for (unsigned long currentCluster = 0; currentCluster < numClusters; currentCluster++) {
        const auto &clusterNeighborList = neighborLists.at(gridIndex1D).at(currentCluster);
        Particle *clusterStart = &currentGrid[currentCluster * clusterSize];
        for (auto neighborClusterStart : clusterNeighborList) {
          traverseClusterPair(clusterStart, neighborClusterStart, clusterSize);
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseClusterPair(
    Particle *clusterStart, Particle *neighborClusterStart, int clusterSize) {
  if (neighborClusterStart == clusterStart) {
    for (int i = 0; i < clusterSize; i++) {
      for (int j = i + 1; j < clusterSize; j++) {
        _functor->AoSFunctor(*(clusterStart + i), *(neighborClusterStart + j), true);
      }
    }
  } else {
    for (int i = 0; i < clusterSize; i++) {
      for (int j = 0; j < clusterSize; j++) {
        _functor->AoSFunctor(*(clusterStart + i), *(neighborClusterStart + j), useNewton3);
      }
    }
  }
}

}  // namespace autopas
