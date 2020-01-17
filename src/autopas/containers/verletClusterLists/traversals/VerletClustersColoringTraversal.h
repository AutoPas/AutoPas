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
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VerletClustersColoringTraversal : public CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                                        public VerletClustersTraversalInterface<typename ParticleCell::ParticleType> {
 private:
  using Particle = typename ParticleCell::ParticleType;

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

  void traverseClusterPairSoA(Particle *clusterStart, Particle *neighborClusterStart);

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
    if constexpr (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
    auto &grids = clusterList.getGrids();
    const auto clusterSize = clusterList.getClusterSize();
    _gridSoAs.resize(grids.size());
    for (size_t gridIndex = 0; gridIndex < grids.size(); gridIndex++) {
      // Load particles into SoA
      auto &grid = grids[gridIndex];
      _functor->SoALoader(grid, _gridSoAs[gridIndex]);

      // Build _clusterToGridIndexMap
      const size_t numClustersInGrid = grid.numParticles() / clusterSize;
      for (size_t clusterIndex = 0; clusterIndex < numClustersInGrid; clusterIndex++) {
        Particle *clusterStart = &grid[clusterIndex * clusterSize];
        auto clusterStartIndex = clusterIndex * clusterSize;
        auto clusterEndIndex = clusterStartIndex + clusterSize;
        // Emplace SoAView on cluster at key clusterStart
        _clusterToSoAViewMap.emplace(std::piecewise_construct, std::forward_as_tuple(clusterStart),
                                     std::forward_as_tuple(&_gridSoAs[gridIndex], clusterStartIndex, clusterEndIndex));
      }
    }
  }

  void endTraversal() override {
    if constexpr (dataLayout != DataLayoutOption::soa) return;

    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;
    auto &grids = clusterList.getGrids();
    for (size_t i = 0; i < grids.size(); i++) {
      _functor->SoAExtractor(grids[i], _gridSoAs[i]);
    }
  }

  void traverseParticlePairs() override {
    auto &clusterList = *VerletClustersTraversalInterface<Particle>::_verletClusterLists;

    const auto gridSideLength = clusterList.getGridSideLength();
    const auto cutoff = clusterList.getCutoff();
    const auto skin = clusterList.getSkin();
    const auto gridsPerColoringCell = std::ceil((cutoff + skin) / gridSideLength);
    std::array<unsigned long, 3> coloringCellsPerDim{};
    for (int i = 0; i < 3; i++) {
      coloringCellsPerDim[i] =
          static_cast<unsigned long>(std::ceil(clusterList.getCellsPerDimension()[i] / gridsPerColoringCell));
    }

    auto loopBody = [this, gridsPerColoringCell](unsigned long x, unsigned long y, unsigned long z) {
      processColorCell(x, y, z, gridsPerColoringCell);
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
  /**
   * The SoAs for each grid.
   */
  std::vector<SoA<typename Particle::SoAArraysType>, AlignedAllocator<SoA<typename Particle::SoAArraysType>>> _gridSoAs;
  /**
   * A map from the pointer to the start of a cluster to a SoAView on the cluster.
   */
  std::unordered_map<Particle *, SoAView<typename Particle::SoAArraysType>> _clusterToSoAViewMap;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
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
      const auto y = yColorCell * gridsPerColoringCell + yInner;
      const auto x = xColorCell * gridsPerColoringCell + xInner;

      // Not every coloring cell has to have gridsPerColoringCell grids in every direction.
      if (x >= cellsPerDim[0] or y >= cellsPerDim[1]) {
        continue;
      }
      auto gridIndex1D = VerletClusterMaths::index1D(x, y, cellsPerDim);

      auto &currentGrid = grids[gridIndex1D];
      auto numClusters = currentGrid.numParticles() / clusterSize;
      for (unsigned long currentCluster = 0; currentCluster < numClusters; currentCluster++) {
        const auto &clusterNeighborList = neighborLists.at(gridIndex1D).at(currentCluster);
        Particle *clusterStart = &currentGrid[currentCluster * clusterSize];
        for (auto neighborClusterStart : clusterNeighborList) {
          if constexpr (dataLayout == DataLayoutOption::aos) {
            traverseClusterPairAoS(clusterStart, neighborClusterStart, clusterSize);
          } else {
            traverseClusterPairSoA(clusterStart, neighborClusterStart);
          }
        }
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseClusterPairAoS(
    Particle *clusterStart, Particle *neighborClusterStart, int clusterSize) {
  const bool isClusterInteractionWithItself = neighborClusterStart == clusterStart;
  for (int i = 0; i < clusterSize; i++) {
    if (isClusterInteractionWithItself) {
      // Always use newton 3 for interactions within one cluster.
      for (int j = i + 1; j < clusterSize; j++) {
        _functor->AoSFunctor(*(clusterStart + i), *(neighborClusterStart + j), true);
      }
    } else {
      // Calculate interactions between two different clusters.
      for (int j = 0; j < clusterSize; j++) {
        _functor->AoSFunctor(*(clusterStart + i), *(neighborClusterStart + j), useNewton3);
      }
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void VerletClustersColoringTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseClusterPairSoA(
    Particle *clusterStart, Particle *neighborClusterStart) {
  auto clusterView = _clusterToSoAViewMap[clusterStart];

  const bool isClusterInteractionWithItself = clusterStart == neighborClusterStart;
  if (isClusterInteractionWithItself) {
    _functor->SoAFunctor(clusterView, useNewton3);
  } else {
    auto neighborClusterView = _clusterToSoAViewMap[neighborClusterStart];
    _functor->SoAFunctor(clusterView, neighborClusterView, useNewton3);
  }
}

}  // namespace autopas
