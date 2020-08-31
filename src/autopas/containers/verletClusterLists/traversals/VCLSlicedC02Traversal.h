/**
 * @file VerletClustersCSlicedTraversal.h
 * @author fischerv
 * @date 11 Aug 2020
 */

#pragma once

#include "autopas/containers/cellPairTraversals/SlicedC02BasedTraversal.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletClusterLists/traversals/VCLClusterFunctor.h"
#include "autopas/containers/verletClusterLists/traversals/VCLTraversalInterface.h"

namespace autopas {

/**
 * This traversal splits the domain into slices along the longer dimension among x and y.
 * The sliced are divided into two colors which are separately processed to prevent race
 * conditions.
 *
 * @tparam ParticleCell
 * @tparam PairwiseFunctor
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class VCLSlicedC02Traversal : public SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>,
                              public VCLTraversalInterface<typename ParticleCell::ParticleType> {
 private:
  using Particle = typename ParticleCell::ParticleType;

  PairwiseFunctor *_functor;
  internal::VCLClusterFunctor<Particle, PairwiseFunctor, dataLayout, useNewton3> _clusterFunctor;

  void processBaseStep(unsigned long x, unsigned long y) {
    auto &clusterList = *VCLTraversalInterface<Particle>::_verletClusterLists;
    auto &currentTower = clusterList.getTowerAtCoordinates(x, y);
    for (auto &cluster : currentTower.getClusters()) {
      _clusterFunctor.traverseCluster(cluster);
      for (auto *neighborCluster : cluster.getNeighbors()) {
        _clusterFunctor.traverseClusterPair(cluster, *neighborCluster);
      }
    }
  }

 public:
  /**
   * Constructor of the VCLSlicedC02Traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor to use for the traveral.
   * @param interactionLength Interaction length (cutoff + skin).
   * @param cellLength cell length.
   * @param clusterSize the number of particles per cluster.
   */
  explicit VCLSlicedC02Traversal(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                                 const double interactionLength, const std::array<double, 3> &cellLength,
                                 size_t clusterSize)
      : SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>(dims, pairwiseFunctor,
                                                                                       interactionLength, cellLength),
        _functor(pairwiseFunctor),
        _clusterFunctor(pairwiseFunctor, clusterSize) {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::vcl_sliced_c02; }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

  void initTraversal() override {
    if constexpr (dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->loadParticlesIntoSoAs(_functor);
    }
    SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::initTraversal();
  }

  void endTraversal() override {
    if constexpr (dataLayout == DataLayoutOption::soa) {
      VCLTraversalInterface<Particle>::_verletClusterLists->extractParticlesFromSoAs(_functor);
    }
    SlicedC02BasedTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::endTraversal();
  }

  void traverseParticlePairs() override {
    this->template cSlicedTraversal</*allCells*/ true>(
        [&](unsigned long x, unsigned long y, unsigned long z) { processBaseStep(x, y); });
  }
};
}  // namespace autopas
