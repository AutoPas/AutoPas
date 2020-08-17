/**
 * @file VCLC06Traversal.h
 * @author humig
 * @date 27.06.19
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CBasedTraversal.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"

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
class VCLC06Traversal : public CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                        public VCLTraversalInterface<typename ParticleCell::ParticleType> {
 private:
  using Particle = typename ParticleCell::ParticleType;

  /**
   * Each base step looks like this:
   *    X C N  Colors:  1 2 3
   *    N N N           4 5 6
   * Where C is the current cell, N are the neighbor cells that is worked on, and X is not worked on. The neighbor list
   * with newton 3 of the VerletClusterLists container is build in a way that the neighbor lists already contain only
   * the neighbor clusters of these cells.s
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

 public:
  /**
   * Constructor of the VCLClusterIterationTraversal.
   * @param pairwiseFunctor The functor to use for the traversal.
   * @param clusterSize Number of particles per cluster.
   */
  explicit VCLC06Traversal(PairwiseFunctor *pairwiseFunctor, size_t clusterSize)
      : CBasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>({0, 0, 0}, pairwiseFunctor, 0, {}),
        _functor(pairwiseFunctor),
        _clusterFunctor(pairwiseFunctor, clusterSize) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_c06; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  [[nodiscard]] bool isApplicable() const override {
    return (dataLayout == DataLayoutOption::aos || dataLayout == DataLayoutOption::soa);
  }

  void initTraversal() override {
    if constexpr (dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->loadParticlesIntoSoAs(_functor);
    }
  }

  void endTraversal() override {
    if constexpr (dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->extractParticlesFromSoAs(_functor);
    }
  }

  void traverseParticlePairs() override {
    auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;

    const auto towersPerColoringCell = clusterList.getNumTowersPerInteractionLength();
    std::array<unsigned long, 2> coloringCellsPerDim{};
    for (int i = 0; i < 2; i++) {
      coloringCellsPerDim[i] =
          static_cast<unsigned long>(std::ceil(clusterList.getTowersPerDimension()[i] / (double)towersPerColoringCell));
    }

    auto loopBody = [this, towersPerColoringCell](unsigned long x, unsigned long y, unsigned long z) {
      processColorCell(x, y, z, towersPerColoringCell);
    };

    // localStride is necessary because stride is constexpr and cTraversal() wants a const &
    auto localStride = _stride;
    this->cTraversal(std::forward<decltype(loopBody)>(loopBody), {coloringCellsPerDim[0], coloringCellsPerDim[1], 1},
                     localStride);
  }

 private:
  PairwiseFunctor *_functor;
  internal::VCLClusterFunctor<Particle, PairwiseFunctor, dataLayout, useNewton3> _clusterFunctor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
void VCLC06Traversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::processColorCell(
    unsigned long xColorCell, unsigned long yColorCell, unsigned long zColorCell, int towersPerColoringCell) {
  // We are only doing a 2D coloring.
  if (zColorCell != 0) {
    autopas::utils::ExceptionHandler::exception("Coloring should only be 2D, not in z-direction!");
  }

  auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;
  const auto towersPerDim = clusterList.getTowersPerDimension();

  for (int yInner = 0; yInner < towersPerColoringCell; yInner++) {
    for (int xInner = 0; xInner < towersPerColoringCell; xInner++) {
      const auto y = yColorCell * towersPerColoringCell + yInner;
      const auto x = xColorCell * towersPerColoringCell + xInner;

      // Not every coloring cell has to have gridsPerColoringCell grids in every direction.
      if (x >= towersPerDim[0] or y >= towersPerDim[1]) {
        continue;
      }

      auto &currentTower = clusterList.getTowerAtCoordinates(x, y);
      for (auto &cluster : currentTower.getClusters()) {
        _clusterFunctor.traverseCluster(cluster);

        for (auto *neighborClusterPtr : cluster.getNeighbors()) {
          _clusterFunctor.traverseClusterPair(cluster, *neighborClusterPtr);
        }
      }
    }
  }
}

}  // namespace autopas
