/**
 * @file LCC01CudaTraversal.h
 * @author jspahl
 * @date 11.03.2019
 */

#pragma once

#include "LCTraversalInterface.h"
#include "autopas/containers/cellPairTraversals/CellPairTraversal.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/utils/CudaDeviceVector.h"
#include "autopas/utils/ThreeDimensionalMapping.h"
#include "autopas/utils/WrapOpenMP.h"
#if defined(AUTOPAS_CUDA)
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
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
class LCC01CudaTraversal : public CellPairTraversal<ParticleCell>, public LCTraversalInterface<ParticleCell> {
 public:
  /**
   * Constructor of the c01 traversal.
   * @param dims The dimensions of the cellblock, i.e. the number of cells in x,
   * y and z direction.
   * @param pairwiseFunctor The functor that defines the interaction of two particles.
   * @param interactionLength The interaction length.
   * @param cellLengths The lengths of the cell.
   */
  explicit LCC01CudaTraversal(const std::array<uint64_t, 3> &dims, PairwiseFunctor *pairwiseFunctor,
                              double interactionLength, std::array<double, 3> cellLengths)
      : CellPairTraversal<ParticleCell>(dims),
        _functor{pairwiseFunctor},
        _interactionLength{interactionLength},
        _cellLengths{cellLengths} {
    computeOffsets();
  }

  /**
   * Computes pairs used in processBaseCell()
   */
  void computeOffsets();

  void traverseParticlePairs() override;

  void initTraversal() override {}

  void endTraversal() override {}

  [[nodiscard]] TraversalOption getTraversalType() const override { return TraversalOption::lc_c01_cuda; }

  /**
   * Cuda traversal is only usable if using a GPU.
   * @return
   */
  [[nodiscard]] bool isApplicable() const override {
#if defined(AUTOPAS_CUDA)
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    // cuda only:
    bool applicable = dataLayout == DataLayoutOption::cuda;
    // only if we actually have some devices:
    applicable &= nDevices > 0;
    for (auto cellLength : _cellLengths) {
      // we only interact neighboring cells, therefore the cell length has to be bigger than the interaction length.
      /// @todo reenable when https://github.com/AutoPas/AutoPas/issues/417 is done.
      applicable &= cellLength >= _interactionLength;
    }
    // currently newton3 support for this traversal is buggy
    /// @todo reenable when fixed, see: https://github.com/AutoPas/AutoPas/issues/420
    applicable &= not useNewton3;
    return applicable;
#else
    return false;
#endif
  }

  [[nodiscard]] DataLayoutOption getDataLayout() const override { return dataLayout; }

  [[nodiscard]] bool getUseNewton3() const override { return useNewton3; }

 private:
  /**
   * Pairs for processBaseCell().
   */
  std::vector<int> _cellOffsets;

  /**
   * Pairwise Functor to be used
   */
  PairwiseFunctor *_functor;

  /**
   * SoA Storage cell for all cells and device Memory
   */
  ParticleCell _storageCell;

  /**
   * Non Halo Cell Ids
   */
  utils::CudaDeviceVector<unsigned int> _nonHaloCells;

  /**
   * device cell sizes storage
   */
  utils::CudaDeviceVector<size_t> _deviceCellSizes;

  /**
   * Interaction length.
   */
  double _interactionLength;

  /**
   * Cell lengths.
   */
  std::array<double, 3> _cellLengths;
};

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC01CudaTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::computeOffsets() {
  for (int z = -1; z <= 1; ++z) {
    for (int y = -1; y <= 1; ++y) {
      for (int x = -1; x <= 1; ++x) {
        int offset = (z * this->_cellsPerDimension[1] + y) * this->_cellsPerDimension[0] + x;
        if (not useNewton3) {
          // all offsets for useNewton3 == false
          _cellOffsets.push_back(offset);
        } else {
          // only positive offsets for useNewton3 == true
          /// @todo this is wrong! especially for halo cells and has to be adapted for correct calculations!
          /// see https://github.com/AutoPas/AutoPas/issues/420
          if (offset > 0) {
            _cellOffsets.push_back(offset);
          }
        }
      }
    }
  }

  std::vector<unsigned int> nonHaloCells((this->_cellsPerDimension[0] - 2) * (this->_cellsPerDimension[1] - 2) *
                                         (this->_cellsPerDimension[2] - 2));
  const uint64_t end_y = this->_cellsPerDimension[1] - 1;
  const uint64_t end_z = this->_cellsPerDimension[2] - 1;
  const uint64_t length_x = this->_cellsPerDimension[0] - 2;

  auto it = nonHaloCells.begin();
  for (uint64_t z = 1; z < end_z; ++z) {
    for (uint64_t y = 1; y < end_y; ++y) {
      std::iota(it, it + length_x, utils::ThreeDimensionalMapping::threeToOneD(1ul, y, z, this->_cellsPerDimension));
      it += length_x;
    }
  }
#if defined(AUTOPAS_CUDA)
  if (dataLayout == DataLayoutOption::cuda) {
    if (_functor->getCudaWrapper())
      _functor->getCudaWrapper()->loadLinkedCellsOffsets(_cellOffsets.size(), _cellOffsets.data());
    _nonHaloCells.copyHostToDevice(nonHaloCells.size(), nonHaloCells.data());
  }
#endif
}

template <class ParticleCell, class PairwiseFunctor, DataLayoutOption::Value dataLayout, bool useNewton3>
inline void LCC01CudaTraversal<ParticleCell, PairwiseFunctor, dataLayout, useNewton3>::traverseParticlePairs() {
  if (not(dataLayout == DataLayoutOption::cuda)) {
    utils::ExceptionHandler::exception(
        "The Cuda traversal cannot work with Data Layouts other than DataLayoutOption::cuda!");
  }
#if defined(AUTOPAS_CUDA)
  auto &cells = *(this->_cells);
  // load CUDA SOA
  std::vector<size_t> cellSizePartialSum = {0};
  size_t maxParticlesInCell = 0;

  for (size_t i = 0; i < cells.size(); ++i) {
    _functor->SoALoader(cells[i], _storageCell._particleSoABuffer, cellSizePartialSum.back());
    const size_t size = cells[i].numParticles();

    maxParticlesInCell = std::max(maxParticlesInCell, size);

    cellSizePartialSum.push_back(cellSizePartialSum.back() + size);
  }
  if (maxParticlesInCell == 0) {
    return;
  }

  unsigned int requiredThreads = ((maxParticlesInCell - 1) / 32 + 1) * 32;

  _deviceCellSizes.copyHostToDevice(cellSizePartialSum.size(), cellSizePartialSum.data());
  _functor->deviceSoALoader(_storageCell._particleSoABuffer, _storageCell._particleSoABufferDevice);

  // wait for copies to be done
  utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());

  auto cudaSoA = _functor->createFunctorCudaSoA(_storageCell._particleSoABufferDevice);

  if (useNewton3) {
    _functor->getCudaWrapper()->LinkedCellsTraversalN3Wrapper(cudaSoA.get(), requiredThreads, _nonHaloCells.size(),
                                                              _nonHaloCells.get(), _deviceCellSizes.size(),
                                                              _deviceCellSizes.get(), 0);
  } else {
    _functor->getCudaWrapper()->LinkedCellsTraversalNoN3Wrapper(cudaSoA.get(), requiredThreads, _nonHaloCells.size(),
                                                                _nonHaloCells.get(), _deviceCellSizes.size(),
                                                                _deviceCellSizes.get(), 0);
  }
  utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());

  // Extract
  _functor->deviceSoAExtractor(_storageCell._particleSoABuffer, _storageCell._particleSoABufferDevice);

  for (size_t i = 0; i < cells.size(); ++i) {
    _functor->SoAExtractor(cells[i], _storageCell._particleSoABuffer, cellSizePartialSum[i]);
  }

  utils::CudaExceptionHandler::checkErrorCode(cudaDeviceSynchronize());
#endif
}

}  // namespace autopas
