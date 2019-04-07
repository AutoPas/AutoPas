/**
 * @file VerletClusterClusterCudaTraversal.h
 * @author jspahl
 * @date 25.03.2019
 */

#pragma once

#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/pairwiseFunctors/CellFunctor.h"
#include "autopas/utils/CudaDeviceVector.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"
#if defined(AUTOPAS_CUDA)
#include "autopas/pairwiseFunctors/LJFunctor.h"
#include "autopas/pairwiseFunctors/LJFunctorCuda.cuh"
#include "autopas/utils/CudaExceptionHandler.h"
#include "autopas/utils/CudaStreamHandler.h"
#endif
namespace autopas {

/**
 * This class provides the c01 traversal on the GPU.
 *
 * The traversal calculates all cells in parallel
 *
 * @tparam ParticleCell the type of cells
 * @tparam PairwiseFunctor The functor that defines the interaction of two particles.
 * @tparam DataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
class VerletClusterClusterCuda : public CellPairTraversal<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   */
  explicit VerletClusterClusterCuda(const std::array<unsigned long, 3> &dims, PairwiseFunctor *pairwiseFunctor)
      : CellPairTraversal<ParticleCell>(dims), _functor(pairwiseFunctor) {}

  TraversalOption getTraversalType() override { return TraversalOption::dummyTraversal; }

  /**
   * Cuda traversal is only usable if using a GPU.
   * @return true if applicable
   */
  bool isApplicable() override {
#if defined(AUTOPAS_CUDA)
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    return (DataLayout == DataLayoutOption::cuda) && (nDevices > 0);
#else
    return false;
#endif
  }

  DataLayoutOption requiredDataLayout() override { return DataLayoutOption::aos; }

  /**
   * This function interacts all cells with the other cells with their index in neighborCellIds
   * @param cells containing the particles
   * @param neighborCellIds Stores the neighbor ids for each cell in cells
   */
  void traverseCellPairs(std::vector<ParticleCell> &cells, std::vector<std::vector<size_t>> &neighborCellIds);

 private:
  /**
   * Pairwise Functor to be used
   */
  PairwiseFunctor *_functor;

  /**
   * SoA Storage cell for all cells and device Memory
   */
  ParticleCell _storageCell;

  /**
   * neighbor cells
   */
  utils::CudaDeviceVector<unsigned int> _neighborCellIds;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3>
void VerletClusterClusterTraversal<ParticleCell, PairwiseFunctor, DataLayout, useNewton3>::traverseCellPairs(
    std::vector<ParticleCell> &cells, std::vector<std::vector<size_t>> &neighborCellIds) {
  if (not this->isApplicable()) {
    utils::ExceptionHandler::exception(
        "The Cuda traversal cannot work with Data Layouts other than DataLayoutOption::cuda!");
  }
  unsigned int maxSize = 0;
  for (auto cell : neighborCellIds) {
    std::max(maxSize, cell.size());
  }
  std::vector<unsigned int> neighborMatrix(maxSize * neighborCellIds.size() + 1);

  size_t i = 0;
  for (auto cellNs : neighborCellIds) {
    size_t nid;
    for (nid = 0; nid < cellNs.size(); ++nid) {
      neighborMatrix[i++] = cellNs[nid];
    }
    for (; pid < maxSize; ++nid) {
      neighborMatrix[i++] = UINT_MAX;
    }
  }

#if defined(AUTOPAS_CUDA)
  // Load
  for (size_t i = 0; i < cells.size(); ++i) {
    _functor->SoALoader(cells[i], _storageCell._particleSoABuffer, _storageCell.numParticles());
  }

  _functor->getCudaWrapper()->setNumThreads(cells.front().size());
  _functor->deviceSoALoader(_storageCell._particleSoABuffer, _storageCell._particleSoABufferDevice);
  _neighborCellIds.copyHostToDevice(neighborMatrix.size(), neighborMatrix.data());

  utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());

  // calculate
  if (useNewton3) {
  } else {
    _functor->getCudaWrapper()->CellVerletTraversalNoN3Wrapper(
        _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::posX>().get(),
        _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::posY>().get(),
        _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::posZ>().get(),
        _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::forceX>().get(),
        _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::forceY>().get(),
        _storageCell._particleSoABufferDevice.template get<ParticleCell::ParticleType::AttributeNames::forceZ>().get(),
        cells.size(), _nonHaloCells.get(), _deviceCellSizes.size(), _deviceCellSizes.get(), _deviceCellOffsets.size(),
        _deviceCellOffsets.get(), 0);
  }
  utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());
  // Extract
  _functor->deviceSoAExtractor(_storageCell._particleSoABuffer, _storageCell._particleSoABufferDevice);
  for (size_t i = 0; i < cells.size(); ++i) {
    _functor->SoAExtractor(cells[i], _storageCell._particleSoABuffer, cells[i].size());
  }
  utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());
#else
  utils::ExceptionHandler::exception("Autopas was compiled without cuda support!");
#endif
}

}  // namespace autopas
