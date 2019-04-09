/**
 * @file VerletClusterClusterTraversal.h
 * @author jspahl
 * @date 25.3.19
 */

#pragma once

#include <vector>
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#if defined(AUTOPAS_CUDA)
#include "cuda_runtime.h"
#endif

namespace autopas {

/**
 * This Traversal is used to interact all clusters in VerletClusterCluster Container
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class VerletClusterClusterTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor for the VerletClusterClusterTraversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  VerletClusterClusterTraversal(PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>({2, 1, 1}),
        _cellFunctor(CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout,
                                 useNewton3, true>(pairwiseFunctor)) {}

  TraversalOption getTraversalType() override { return TraversalOption::ClusterToClusterVerlet; }

  bool isApplicable() override {
#if defined(AUTOPAS_CUDA)
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    return (not(DataLayout == DataLayoutOption::cuda)) || nDevices > 0;
#else
    return true;
#endif
  }

  void initTraversal(std::vector<ParticleCell> &cells) override {}

  void endTraversal(std::vector<ParticleCell> &cells) override {}

  /**
   * This function interacts all cells with the other cells with their index in neighborCellIds
   * @param cells containing the particles
   * @param neighborCellIds Stores the neighbor ids for each cell in cells
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells, std::vector<std::vector<size_t>> &neighborCellIds);

 private:
  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout, useNewton3, true>
      _cellFunctor;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void VerletClusterClusterTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells, std::vector<std::vector<size_t>> &neighborCellIds) {
  for (size_t i = 0; i < cells.size(); ++i) {
    for (auto &j : neighborCellIds[i]) {
      _cellFunctor.processCellPair(cells[i], cells[j]);
    }
    _cellFunctor.processCell(cells[i]);
  }
}

}  // namespace autopas
