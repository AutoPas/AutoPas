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
  using Particle = typename ParticleCell::ParticleType;
  using index_t = typename VerletClusterMaths::index_t;
  /**
   * This stride is determined by the way the neighbor list with newton 3 of the VerletClusterLists container is build.
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
  void traverseClusterPairAoS(Particle *clusterStart, Particle *neighborClusterStart, int clusterSize);

  void traverseClusterPairSoA(index_t gridIndex, index_t neighborGridIndex, index_t clusterNum,
                              index_t neighborClusterNum, int clusterSize);

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

  void initTraversal(std::vector<ParticleCell> &cells) override {}
  void endTraversal(std::vector<ParticleCell> &cells) override {}

  void initClusterTraversal() override {
    if constexpr (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
    auto &grids = clusterList.getGrids();
    _gridSoAs.resize(grids.size());
    for (index_t gridIndex = 0; gridIndex < grids.size(); gridIndex++) {
      // Load particles into SoA
      auto &grid = grids[gridIndex];
      _functor->SoALoader(grid, _gridSoAs[gridIndex]);

      // Build _clusterToGridIndexMap
      const index_t numClustersInGrid = grid.numParticles() / clusterList.getClusterSize();
      for (index_t clusterIndex = 0; clusterIndex < numClustersInGrid; clusterIndex++) {
        Particle *clusterStart = &grid[clusterIndex * clusterList.getClusterSize()];
        _clusterToGridIndexMap[clusterStart] = {gridIndex, clusterIndex};
      }
    }
  }

  void endClusterTraversal() override {
    if constexpr (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
    auto &grids = clusterList.getGrids();
    for (index_t i = 0; i < grids.size(); i++) {
      _functor->SoAExtractor(grids[i], _gridSoAs[i]);
    }
  }

  /**
   * @copydoc VerletClustersTraversalInterface::traverseParticlePairs
   */
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
    assert(end[2] == 1);

    // localStride is necessary because stride is constexpr and cTraversal() wants a const &
    auto localStride = stride;
    this->cTraversal(std::forward<decltype(loopBody)>(loopBody), end, localStride);
  }

 private:
  PairwiseFunctor *_functor;
  /**
   * The SoAs for each grid.
   */
  std::vector<SoA<typename Particle::SoAArraysType>> _gridSoAs;
  /**
   * A map from the pointer to the start of a cluster to a pair consisting of the grid index and the cluster index of
   * the cluster in the SoA. The n-th cluster in a grid has cluster index n.
   */
  std::unordered_map<Particle *, std::pair<index_t, index_t>> _clusterToGridIndexMap;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processColorCell(
    unsigned long xColorCell, unsigned long yColorCell, unsigned long zColorCell, int gridsPerColoringCell) {
  // We are only doing a 2D coloring.
  assert(zColorCell == 0);

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
      for (unsigned int currentCluster = 0; currentCluster < numClusters; currentCluster++) {
        const auto &clusterNeighborList = neighborLists.at(gridIndex1D).at(currentCluster);
        [[maybe_unused]] Particle *clusterStart = &currentGrid[currentCluster * clusterSize];
        for (auto neighborClusterStart : clusterNeighborList) {
          if constexpr (dataLayout == DataLayoutOption::aos) {
            traverseClusterPairAoS(clusterStart, neighborClusterStart, clusterSize);
          } else {
            auto [neighborGridIndex, neighborClusterIndex] = _clusterToGridIndexMap[neighborClusterStart];
            traverseClusterPairSoA(gridIndex1D, neighborGridIndex, currentCluster, neighborClusterIndex, clusterSize);
          }
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseClusterPairAoS(
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

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseClusterPairSoA(
    index_t gridIndex, index_t neighborGridIndex, index_t clusterNum, index_t neighborClusterNum, int clusterSize) {
  if (gridIndex == neighborGridIndex and clusterNum == neighborClusterNum) {
    _functor->SoAFunctor(_gridSoAs[gridIndex], clusterNum * clusterSize, clusterSize, useNewton3);
  } else {
    _functor->SoAFunctor(_gridSoAs[gridIndex], clusterNum * clusterSize, clusterSize, _gridSoAs[neighborGridIndex],
                         neighborClusterNum * clusterSize, clusterSize, useNewton3);
  }
}

}  // namespace autopas
