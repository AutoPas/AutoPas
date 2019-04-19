/**
 * @file VerletClusterCellsTraversal.h
 * @author jspahl
 * @date 25.3.19
 */

#pragma once

#include <vector>
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/containers/cellPairTraversals/VerletClusterTraversalInterface.h"
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
class VerletClusterCellsTraversal : public CellPairTraversal<ParticleCell>,
                                    public VerletClusterTraversalInterface<ParticleCell> {
  using Particle = typename ParticleCell::ParticleType;

 public:
  /**
   * Constructor for the VerletClusterClusterTraversal.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  VerletClusterCellsTraversal(PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>({2, 1, 1}),
        _functor(pairwiseFunctor),
        _cellFunctor(CellFunctor<typename ParticleCell::ParticleType, ParticleCell, PairwiseFunctor, DataLayout,
                                 useNewton3, true>(pairwiseFunctor)),
        _neighborMatrixDim(0),
        _clusterSize(0) {}

  TraversalOption getTraversalType() override { return TraversalOption::verletClusterCellsTraversal; }

  bool isApplicable() override {
    int nDevices = 0;
#if defined(AUTOPAS_CUDA)
    cudaGetDeviceCount(&nDevices);
#endif
    return nDevices > 0;
  }

  void rebuild(const std::array<unsigned long, 3> &dims, unsigned int clusterSize, std::vector<ParticleCell> &cells,
               std::vector<std::array<typename Particle::ParticleFloatingPointType, 6>> boundingBoxes,
               typename Particle::ParticleFloatingPointType distance) override {
    this->_cellsPerDimension = dims;
    _clusterSize = clusterSize;
    const size_t cellsSize = cells.size();
    _neighborCellIds.clear();
    _neighborCellIds.resize(cellsSize, {});

    if (DataLayout == DataLayoutOption::aos or DataLayout == DataLayoutOption::soa) {
      for (size_t i = 0; i < cellsSize; ++i) {
        for (size_t j = i + 1; j < cellsSize; ++j) {
          if (boxesOverlap(boundingBoxes[i], boundingBoxes[j], distance)) {
            _neighborCellIds[i].push_back(j);
          }
        }
      }
      return;
    } else if (DataLayout == DataLayoutOption::cuda) {
      for (size_t i = 0; i < cellsSize; ++i) {
        for (size_t j = i + 1; j < cellsSize; ++j) {
          if (boxesOverlap(boundingBoxes[i], boundingBoxes[j], distance)) {
            if (useNewton3) {
              // load balance greedy
              if (_neighborCellIds[i].size() < _neighborCellIds[j].size())
                _neighborCellIds[i].push_back(j);
              else
                _neighborCellIds[j].push_back(i);

            } else {
              _neighborCellIds[i].push_back(j);
              _neighborCellIds[j].push_back(i);
            }
          }
        }
        // self
        if (not useNewton3) _neighborCellIds[i].push_back(i);
        _neighborMatrixDim = std::max(_neighborCellIds[i].size(), _neighborMatrixDim);
      }
      ++_neighborMatrixDim;

      std::vector<unsigned int> neighborMatrix;
      for (auto &neighborCellId : _neighborCellIds) {
        size_t j = 0;
        for (; j < neighborCellId.size(); ++j) {
          neighborMatrix.push_back(neighborCellId[j]);
        }

        for (; j < _neighborMatrixDim; ++j) {
          neighborMatrix.push_back(UINT_MAX);
        }
      }
#ifdef AUTOPAS_CUDA
      _neighborMatrix.copyHostToDevice(neighborMatrix.size(), neighborMatrix.data());
      utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());
#endif
      return;
    }
  }

  void initTraversal(std::vector<ParticleCell> &cells) override {
    switch (DataLayout) {
      case DataLayoutOption::aos: {
        return;
      }
      case DataLayoutOption::soa: {
        for (size_t i = 0; i < cells.size(); ++i) {
          _functor->SoALoader(cells[i], cells[i]._particleSoABuffer);
        }

        return;
      }
      case DataLayoutOption::cuda: {
        size_t size = 0;
        for (auto &cell : cells) {
          _functor->SoALoader(cell, _storageCell._particleSoABuffer, size);
          size += cell.numParticles();
        }

        _functor->deviceSoALoader(_storageCell._particleSoABuffer, _storageCell._particleSoABufferDevice);
#ifdef AUTOPAS_CUDA
        utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());
#endif
        return;
      }
    }
  }

  void endTraversal(std::vector<ParticleCell> &cells) override {
    switch (DataLayout) {
      case DataLayoutOption::aos: {
        return;
      }
      case DataLayoutOption::soa: {
        for (size_t i = 0; i < cells.size(); ++i) {
          _functor->SoAExtractor(cells[i], cells[i]._particleSoABuffer);
        }

        return;
      }
      case DataLayoutOption::cuda: {
        // Extract
        _functor->deviceSoAExtractor(_storageCell._particleSoABuffer, _storageCell._particleSoABufferDevice);
#ifdef AUTOPAS_CUDA
        utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());
#endif
        size_t size = 0;
        for (auto &cell : cells) {
          _functor->SoAExtractor(cell, _storageCell._particleSoABuffer, size);
          size += cell.numParticles();
        }
        return;
      }
    }
  }

  /**
   * This function interacts all cells with the other cells with their index in neighborCellIds
   * @param cells containing the particles
   * @param neighborCellIds Stores the neighbor ids for each cell in cells
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells) override {
    switch (DataLayout) {
      case DataLayoutOption::aos: {
        traverseCellPairsCPU(cells);
        return;
      }
      case DataLayoutOption::soa: {
        traverseCellPairsCPU(cells);
        return;
      }
      case DataLayoutOption::cuda: {
        traverseCellPairsGPU(cells);
        return;
      }
    }
  }

 private:
  void traverseCellPairsCPU(std::vector<ParticleCell> &cells);

  void traverseCellPairsGPU(std::vector<ParticleCell> &cells);

  /**
   * Returns true if the two boxes are within distance
   * @param box1
   * @param box2
   * @param distance betwwen the boxes to return true
   * @return true if the boxes are overlapping
   */
  inline bool boxesOverlap(std::array<typename Particle::ParticleFloatingPointType, 6> &box1,
                           std::array<typename Particle::ParticleFloatingPointType, 6> &box2,
                           typename Particle::ParticleFloatingPointType distance) {
    for (int i = 0; i < 3; ++i) {
      if (box1[0 + i] - distance > box2[3 + i] || box1[3 + i] + distance < box2[0 + i]) return false;
    }
    return true;
  }

  /**
   * Pairwise functor used in this traversal
   */
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

  size_t _neighborMatrixDim;
  utils::CudaDeviceVector<unsigned int> _neighborMatrix;

  size_t _clusterSize;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void VerletClusterCellsTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairsCPU(
    std::vector<ParticleCell> &cells) {
  for (size_t i = 0; i < cells.size(); ++i) {
    for (auto &j : _neighborCellIds[i]) {
      _cellFunctor.processCellPair(cells[i], cells[j]);
    }
    _cellFunctor.processCell(cells[i]);
  }
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void VerletClusterCellsTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairsGPU(
    std::vector<ParticleCell> &cells) {
#ifdef AUTOPAS_CUDA
  if (!_functor->getCudaWrapper()) {
    _functor->CudaFunctor(_storageCell._particleSoABufferDevice, useNewton3);
    return;
  }

  auto cudaSoA = _functor->createFunctorCudaSoA(_storageCell._particleSoABufferDevice);

  if (useNewton3) {
    _functor->getCudaWrapper()->CellVerletTraversalN3Wrapper(cudaSoA.get(), cells.size(), _clusterSize,
                                                             _neighborMatrixDim, _neighborMatrix.get(), 0);
  } else {
    _functor->getCudaWrapper()->CellVerletTraversalNoN3Wrapper(cudaSoA.get(), cells.size(), _clusterSize,
                                                               _neighborMatrixDim, _neighborMatrix.get(), 0);
  }
  utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());
#endif
}

}  // namespace autopas
