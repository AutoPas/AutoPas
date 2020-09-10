/**
 * @file EmptyCudaWrapper.cuh
 * @author seckler
 * @date 17.08.20
 */

#pragma once

#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/utils/CudaSoA.h"

/**
 * Empty Wrapper for cuda support.
 */
template <typename floatingPointType>
class EmptyCudaWrapper : public autopas::CudaWrapperInterface<floatingPointType> {
  using floatType = floatingPointType;

 public:
  void setNumThreads(int num_threads) override {}
  void loadConstants(autopas::FunctorCudaConstants<floatType> *constants) override {}
  void SoAFunctorNoN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1, cudaStream_t stream) override {}
  void SoAFunctorNoN3PairWrapper(autopas::FunctorCudaSoA<floatType> *cell1, autopas::FunctorCudaSoA<floatType> *cell2,
                                 cudaStream_t stream) override {}
  void SoAFunctorN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1, cudaStream_t stream) override {}
  void SoAFunctorN3PairWrapper(autopas::FunctorCudaSoA<floatType> *cell1, autopas::FunctorCudaSoA<floatType> *cell2,
                               cudaStream_t stream) override {}
  void LinkedCellsTraversalNoN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1, unsigned int reqThreads,
                                       unsigned int cids_size, unsigned int *cids, unsigned int cellSizes_size,
                                       size_t *cellSizes, cudaStream_t stream) override {}
  void LinkedCellsTraversalN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1, unsigned int reqThreads,
                                     unsigned int cids_size, unsigned int *cids, unsigned int cellSizes_size,
                                     size_t *cellSizes, cudaStream_t stream) override {}
  void CellVerletTraversalNoN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1Base, unsigned int ncells,
                                      unsigned int clusterSize, unsigned int others_size, unsigned int *other_ids,
                                      cudaStream_t stream) override {}
  void CellVerletTraversalN3Wrapper(autopas::FunctorCudaSoA<floatType> *cell1Base, unsigned int ncells,
                                    unsigned int clusterSize, unsigned int others_size, unsigned int *other_ids,
                                    cudaStream_t stream) override {}
  void loadLinkedCellsOffsets(unsigned int offsets_size, int *offsets) override {}
  bool isAppropriateClusterSize(unsigned int clusterSize) const override { return true; }
};