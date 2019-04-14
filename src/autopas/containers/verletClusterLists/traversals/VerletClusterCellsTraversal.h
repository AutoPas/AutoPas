/**
 * @file VerletClusterCellsTraversal.h
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
class VerletClusterCellsTraversal : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor for the VerletClusterClusterTraversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  VerletClusterCellsTraversal(PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>({2, 1, 1}),
        _functor(pairwiseFunctor),
        _cellFunctor(CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout,
                                 useNewton3, true>(pairwiseFunctor)) {}

  TraversalOption getTraversalType() override { return TraversalOption::ClusterToClusterVerlet; }

  bool isApplicable() override {
    int nDevices = 0;
#if defined(AUTOPAS_CUDA)
    cudaGetDeviceCount(&nDevices);
#endif
    return nDevices > 0;
  }

  void rebuild(const std::array<unsigned long, 3> &dims, std::vector<ParticleCell> &cells) {
    this->_cellsPerDimension = dims;
    /*
      for (index_t i = 0; i < cells.size(); ++i) {
        for (index_t j = i + 1; j < cells.size(); ++j) {
          if (boxesOverlap(_boundingBoxes[i], _boundingBoxes[j])) {
            _neighborCellIds[i].push_back(j);
          }
        }
      }
      */
  };

  void initTraversal(std::vector<ParticleCell> &cells) override {}

  void endTraversal(std::vector<ParticleCell> &cells) override {}

  /**
   * This function interacts all cells with the other cells with their index in neighborCellIds
   * @param cells containing the particles
   * @param neighborCellIds Stores the neighbor ids for each cell in cells
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells, std::vector<std::vector<size_t>> &neighborCellIds);

 private:
  void traverseCellPairsCPU(std::vector<ParticleCell> &cells, std::vector<std::vector<size_t>> &neighborCellIds);

  void traverseCellPairsGPU(std::vector<ParticleCell> &cells, std::vector<std::vector<size_t>> &neighborCellIds);

  PairwiseFunctor *_functor;

  /**
   * CellFunctor to be used for the traversal defining the interaction between two cells.
   */
  CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout, useNewton3, true>
      _cellFunctor;

  // id of neighbor clusters of a clusters
  std::vector<std::vector<size_t>> _neighborCellIds;

  /**
   * SoA Storage cell containing SoAs and device Memory
   */
  ParticleCell _storageCell;

  utils::CudaDeviceVector<unsigned int> _neighborMatrix;

  size_t _clusterSize;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void VerletClusterCellsTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairsCPU(
    std::vector<ParticleCell> &cells, std::vector<std::vector<size_t>> &neighborCellIds) {
  switch (DataLayout) {
    case DataLayoutOption::aos: {
      traverseCellPairsCPU(cells, neighborCellIds);
      return;
    }
    case DataLayoutOption::soa: {
      traverseCellPairsCPU(cells, neighborCellIds);
      return;
    }
    case DataLayoutOption::cuda: {
      traverseCellPairsGPU(cells, neighborCellIds);
      return;
    }
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void VerletClusterCellsTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells, std::vector<std::vector<size_t>> &neighborCellIds) {
  for (size_t i = 0; i < cells.size(); ++i) {
    for (auto &j : neighborCellIds[i]) {
      _cellFunctor.processCellPair(cells[i], cells[j]);
    }
    _cellFunctor.processCell(cells[i]);
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void VerletClusterCellsTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairsGPU(
    std::vector<ParticleCell> &cells, std::vector<std::vector<size_t>> &neighborCellIds) {
  size_t size = 0;
  for (auto &cell : cells) {
    _functor->SoALoader(cell, _storageCell._particleSoABuffer, size);
    size += cell.numParticles();
  }
  _clusterSize = cells[0].numParticles();

  _functor->deviceSoALoader(_storageCell._particleSoABuffer, _storageCell._particleSoABufferDevice);

  if (!_functor->getCudaWrapper()) {
    _functor->CudaFunctor(_storageCell._particleSoABufferDevice, useNewton3);
    return;
  }

  size_t maxNumNeighbors = 0;
  for (auto &cellNeighbor : neighborCellIds) {
    maxNumNeighbors = std::max(cellNeighbor.size(), maxNumNeighbors);
  }
  std::vector<unsigned int> neighborMatrix;
  ++maxNumNeighbors;

  for (auto &cellNeighbors : neighborCellIds) {
    const size_t rest = maxNumNeighbors - cellNeighbors.size();
    for (auto &nid : cellNeighbors) {
      neighborMatrix.push_back(nid);
    }
    for (size_t i = 0; i < rest; ++i) {
      neighborMatrix.push_back(UINT_MAX);
    }
  }
  _neighborMatrix.copyHostToDevice(neighborMatrix.size(), neighborMatrix.data());
  utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());

  _functor->getCudaWrapper()->CellVerletTraversalNoN3Wrapper(
      LJFunctorCudaSoA<typename ParticleCell::ParticleType::ParticleFloatingPointType>(
          0,
          _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::posX>().get(),
          _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::posY>().get(),
          _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::posZ>().get(),
          _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::forceX>()
              .get(),
          _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::forceY>()
              .get(),
          _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::forceZ>()
              .get()),
      cells.size(), _clusterSize, _neighborMatrix.get(), 0);
  utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());

  // Extract
  _functor->deviceSoAExtractor(_storageCell._particleSoABuffer, _storageCell._particleSoABufferDevice);
  utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());

  size = 0;
  for (auto &cell : cells) {
    _functor->SoAExtractor(cell, _storageCell._particleSoABuffer, size);
    size += cell.numParticles();
  }
}

}  // namespace autopas
