/**
 * @file VCLC08Traversal.h
 * @author stagrz
 * @date 30.11.2025
 */

#pragma once

#include "autopas/containers/cellTraversals/ColorBasedTraversal.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"

namespace autopas {

/**
 * A traversal for VerletClusterLists that uses a coloring over the grids of the container.
 *
 * The traversal uses a 3D coloring with a stride of x=2, y=2, z=1 so 2*2*1=4 colors.
 *
 * When disabling newton 3, interactions inside a cluster are still calculated using newton 3.
 *
 * @tparam ParticleCell
 * @tparam PairwiseFunctor
 */
template <class ParticleCell, class PairwiseFunctor>
class VCLC08Traversal : public ColorBasedTraversal<ParticleCell, PairwiseFunctor>,
                        public VCLTraversalInterface<typename ParticleCell::ParticleType> {
 private:
  using ParticleType = typename ParticleCell::ParticleType;

  /**
   * Each base step looks like this:
   *     C N  Colors:  1 2
   *     N N           3 4
   * Where C is the current cell, N are the neighbor cells that is worked on, and X is not worked on. The neighbor list
   * with newton 3 of the VerletClusterLists container is build in a way that the neighbor lists already contain only
   * the neighbor clusters of these cells.s
   */
  static constexpr std::array<unsigned long, 3> _stride{2ul, 2ul, 1ul};

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
   * @param dataLayout The data layout with which this traversal should be initialized.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit VCLC08Traversal(PairwiseFunctor *pairwiseFunctor, size_t clusterSize, DataLayoutOption dataLayout,
                           bool useNewton3)
      : ColorBasedTraversal<ParticleCell, PairwiseFunctor>({0, 0, 0}, pairwiseFunctor, 0, {}, dataLayout, useNewton3),
        _functor(pairwiseFunctor),
        _clusterFunctor(pairwiseFunctor, clusterSize, dataLayout, useNewton3) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_c08; }

  [[nodiscard]] bool isApplicable() const override {
    return (this->_dataLayout == DataLayoutOption::aos || this->_dataLayout == DataLayoutOption::soa);
  }

  void initTraversal() override {
    if (this->_dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<ParticleType>::_verletClusterLists->loadParticlesIntoSoAs(_functor);
    }
  }

  void endTraversal() override {
    if (this->_dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<ParticleType>::_verletClusterLists->extractParticlesFromSoAs(_functor);
    }
  }

  void traverseParticles() override {
    auto &clusterList = *VCLTraversalInterface<ParticleType>::_verletClusterLists;

    const auto towersPerColoringCell = clusterList.getNumTowersPerInteractionLength();
    std::array<unsigned long, 2> coloringCellsPerDim{};
    for (int i = 0; i < 2; i++) {
      coloringCellsPerDim[i] =
          static_cast<unsigned long>(std::ceil(clusterList.getTowersPerDimension()[i] / (double)towersPerColoringCell));
    }

    auto loopBody = [this, towersPerColoringCell](unsigned long x, unsigned long y, unsigned long z) {
      processColorCell(x, y, z, towersPerColoringCell);
    };

    // localStride is necessary because stride is constexpr and colorTraversal() wants a const &
    auto localStride = _stride;
    this->colorTraversal(std::forward<decltype(loopBody)>(loopBody),
                         {coloringCellsPerDim[0], coloringCellsPerDim[1], 1}, localStride);
  }

  /**
   * @copydoc autopas::CellTraversal::setSortingThreshold()
   * This traversal does not use the CellFunctor, so the function has no effect here
   */
  void setSortingThreshold(size_t sortingThreshold) override {}

 private:
  PairwiseFunctor *_functor;
  internal::VCLClusterFunctor<ParticleType, PairwiseFunctor> _clusterFunctor;
};

template <class ParticleCell, class PairwiseFunctor>
void VCLC08Traversal<ParticleCell, PairwiseFunctor>::processColorCell(unsigned long xColorCell,
                                                                      unsigned long yColorCell,
                                                                      unsigned long zColorCell,
                                                                      int towersPerColoringCell) {
  auto &clusterList = *VCLTraversalInterface<ParticleType>::_verletClusterLists;
  const auto towersPerDim = clusterList.getTowersPerDimension();
  if (zColorCell != 0) {
    autopas::utils::ExceptionHandler::exception("Coloring should only be 2D, not in z-direction!");
  }
  for (int yInner = 0; yInner < towersPerColoringCell; yInner++) {
    for (int xInner = 0; xInner < towersPerColoringCell; xInner++) {
      const auto y = yColorCell * towersPerColoringCell + yInner;
      const auto x = xColorCell * towersPerColoringCell + xInner;

      // Not every coloring cell has to have gridsPerColoringCell grids in every direction.
      if (x >= towersPerDim[0] || y >= towersPerDim[1]) {
        continue;
      }
      auto &currentTower = clusterList.getTowerByIndex(x, y);
      for (auto clusterIter = this->_useNewton3 ? currentTower.getClusters().begin()
                                                : currentTower.getFirstOwnedCluster();
           clusterIter <
           (this->_useNewton3 ? currentTower.getClusters().end() : currentTower.getFirstTailHaloCluster());
           ++clusterIter) {
        const auto isHaloCluster =
            clusterIter < currentTower.getFirstOwnedCluster() or clusterIter >= currentTower.getFirstTailHaloCluster();
        _clusterFunctor.processCluster(*clusterIter, isHaloCluster);
      }
    }
  }
}

}  // namespace autopas
