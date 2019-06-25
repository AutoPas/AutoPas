/**
 * @file VerletClustersTraversal.h
 * @author humig
 * @date 20.06.19
 */

#pragma once

#include "VerletClustersTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"

namespace autopas {

/**
 * Traversal for VerletClusterLists.
 * @tparam ParticleCell
 * @tparam PairwiseFunctor The type of the functor.
 * @tparam dataLayout The data layout to use. Currently, only AoS is supported.
 * @tparam useNewton3 If newton 3 should be used. Currently, only false is supported.
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption dataLayout, bool useNewton3>
class VerletClustersTraversal : public CellPairTraversal<ParticleCell, dataLayout, useNewton3>,
                                public VerletClustersTraversalInterface<typename ParticleCell::ParticleType> {
  using Particle = typename ParticleCell::ParticleType;
  using index_t = typename VerletClusterMaths::index_t;

  using Super = VerletClustersTraversalInterface<Particle>;

 public:
  /**
   * Constructor of the VerletClustersTraversal.
   * @param pairwiseFunctor The functor to use for the traveral.
   */
  explicit VerletClustersTraversal(PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell, dataLayout, useNewton3>({0, 0, 0}), _functor(pairwiseFunctor) {}

  TraversalOption getTraversalType() const override { return TraversalOption::verletClusters; }

  DataLayoutOption getDataLayout() const override { return dataLayout; }
  bool getUseNewton3() const override { return useNewton3; }
  bool isApplicable() const override {
    return (dataLayout == DataLayoutOption::aos || dataLayout == DataLayoutOption::soa) && not useNewton3;
  }

  void initTraversal(std::vector<ParticleCell> &cells) override {}
  void endTraversal(std::vector<ParticleCell> &cells) override {}

  void initClusterTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;

    std::array<index_t, 3> cellsPerDim = Super::_cellsPerDim;
    index_t numClusters = Super::_numClusters;
    int clusterSize = Super::_clusterSize;
    std::vector<FullParticleCell<Particle>> &grids = *Super::_grids;
    std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap = *Super::_aosToSoaMap;

    _clusterSoAs.resize(numClusters);
    // iterate over all clusters
#if defined(AUTOPAS_OPENMP)
    // @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    for (index_t x = 0; x < cellsPerDim[0]; x++) {
      for (index_t y = 0; y < cellsPerDim[1]; y++) {
        FullParticleCell<Particle> cell{};
        cell.reserve(clusterSize);
        index_t index = VerletClusterMaths::index1D(x, y, cellsPerDim);
        auto &grid = grids[index];

        const index_t numClustersInGrid = grid.numParticles() / clusterSize;
        for (index_t clusterIndex = 0; clusterIndex < numClustersInGrid; clusterIndex++) {
          Particle *clusterStart = &grid[clusterIndex * clusterSize];
          index_t currentClusterIndex = aosToSoaMap[clusterStart];

          // actual loop body
          cell.clear();
          for (int i = 0; i < clusterSize; i++) {
            cell.addParticle(*(clusterStart + i));
          }
          SoA<typename Particle::SoAArraysType> &soa = _clusterSoAs[currentClusterIndex];
          soa.resizeArrays(clusterSize);
          _functor->SoALoader(cell, soa);
        }
      }
    }
  }

  void endClusterTraversal() override {
    if (dataLayout != DataLayoutOption::soa) return;

    std::array<index_t, 3> cellsPerDim = Super::_cellsPerDim;
    int clusterSize = Super::_clusterSize;
    std::vector<FullParticleCell<Particle>> &grids = *Super::_grids;
    std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap = *Super::_aosToSoaMap;

    // iterate over all clusters
#if defined(AUTOPAS_OPENMP)
    // @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    for (index_t x = 0; x < cellsPerDim[0]; x++) {
      for (index_t y = 0; y < cellsPerDim[1]; y++) {
        FullParticleCell<Particle> cell{};
        cell.reserve(clusterSize);
        index_t index = VerletClusterMaths::index1D(x, y, cellsPerDim);
        auto &grid = grids[index];

        const index_t numClustersInGrid = grid.numParticles() / clusterSize;
        for (index_t clusterIndex = 0; clusterIndex < numClustersInGrid; clusterIndex++) {
          Particle *clusterStart = &grid[clusterIndex * clusterSize];
          index_t currentClusterIndex = aosToSoaMap[clusterStart];

          // actual loop body
          cell.clear();
          for (int i = 0; i < clusterSize; i++) {
            cell.addParticle(*(clusterStart + i));
          }
          SoA<typename Particle::SoAArraysType> &soa = _clusterSoAs[currentClusterIndex];
          _functor->SoAExtractor(cell, soa);
          for (int i = 0; i < clusterSize; i++) {
            *(clusterStart + i) = cell[i];
          }
        }
      }
    }
  }

  /**
   * @copydoc VerletClustersTraversalInterface::traverseParticlePairs
   */
  void traverseParticlePairs() override {
    std::array<index_t, 3> cellsPerDim = Super::_cellsPerDim;
    int clusterSize = Super::_clusterSize;
    std::vector<FullParticleCell<Particle>> &grids = *Super::_grids;
    std::vector<std::vector<std::vector<Particle *>>> &neighborLists = *Super::_neighborLists;
    std::unordered_map<Particle *, VerletClusterMaths::index_t> aosToSoaMap = *Super::_aosToSoaMap;

    const index_t end_x = cellsPerDim[0];
    const index_t end_y = cellsPerDim[1];

#if defined(AUTOPAS_OPENMP)
    // @todo: find sensible chunksize
#pragma omp parallel for schedule(dynamic) collapse(2)
#endif
    for (index_t x = 0; x < end_x; x++) {
      for (index_t y = 0; y < end_y; y++) {
        index_t index = VerletClusterMaths::index1D(x, y, cellsPerDim);
        auto &grid = grids[index];
        auto &gridNeighborList = neighborLists[index];

        const index_t numClustersInGrid = grid.numParticles() / clusterSize;
        for (index_t clusterIndex = 0; clusterIndex < numClustersInGrid; clusterIndex++) {
          Particle *clusterStart = &grid[clusterIndex * clusterSize];
          for (auto neighborClusterStart : gridNeighborList[clusterIndex]) {
            // self pair
            if (clusterStart == neighborClusterStart) {
              traverseSingleCluster(clusterStart, clusterSize, aosToSoaMap);
            } else {
              traverseNeighborClusters(clusterStart, neighborClusterStart, clusterSize, aosToSoaMap);
            }
          }
        }
      }
    }
  }

 private:
  void traverseSingleCluster(Particle *clusterStart, int clusterSize,
                             std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap) {
    switch (dataLayout) {
      case DataLayoutOption::aos:
        traverseSingleClusterAoS(clusterStart, clusterSize);
        break;
      case DataLayoutOption::soa:
        traverseSingleClusterSoA(clusterStart, aosToSoaMap);
        break;
      default:
        autopas::utils::ExceptionHandler::exception(
            "Wrong data layout of VerletClustersTraversal. Only AoS and SoA are supported!");
    }
  }

  void traverseSingleClusterAoS(Particle *clusterStart, int clusterSize) {
    for (int i = 0; i < clusterSize; i++) {
      for (int j = i + 1; j < clusterSize; j++) {
        Particle *iParticle = clusterStart + i;
        Particle *jParticle = clusterStart + j;
        _functor->AoSFunctor(*iParticle, *jParticle, useNewton3);
        if (not useNewton3) _functor->AoSFunctor(*jParticle, *iParticle, useNewton3);
      }
    }
  }

  void traverseSingleClusterSoA(Particle *clusterStart,
                                std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap) {
    _functor->SoAFunctor(_clusterSoAs[aosToSoaMap[clusterStart]], useNewton3);
  }

  void traverseNeighborClusters(Particle *firstClusterStart, Particle *secondClusterStart, int clusterSize,
                                std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap) {
    switch (dataLayout) {
      case DataLayoutOption::aos:
        traverseNeighborClustersAoS(firstClusterStart, secondClusterStart, clusterSize);
        break;
      case DataLayoutOption::soa:
        traverseNeighborClustersSoA(firstClusterStart, secondClusterStart, aosToSoaMap);
        break;
      default:
        autopas::utils::ExceptionHandler::exception(
            "Wrong data layout of VerletClustersTraversal. Only AoS and SoA are supported!");
    }
  }

  void traverseNeighborClustersAoS(Particle *firstClusterStart, Particle *secondClusterStart, int clusterSize) {
    for (int i = 0; i < clusterSize; i++) {
      for (int j = 0; j < clusterSize; j++) {
        Particle *iParticle = firstClusterStart + i;
        Particle *jParticle = secondClusterStart + j;
        _functor->AoSFunctor(*iParticle, *jParticle, useNewton3);
      }
    }
  }

  void traverseNeighborClustersSoA(Particle *firstClusterStart, Particle *secondClusterStart,
                                   std::unordered_map<Particle *, VerletClusterMaths::index_t> &aosToSoaMap) {
    _functor->SoAFunctor(_clusterSoAs[aosToSoaMap[firstClusterStart]], _clusterSoAs[aosToSoaMap[secondClusterStart]],
                         useNewton3);
  }

 private:
  PairwiseFunctor *_functor;

  std::vector<SoA<typename Particle::SoAArraysType>> _clusterSoAs;
};
}  // namespace autopas
