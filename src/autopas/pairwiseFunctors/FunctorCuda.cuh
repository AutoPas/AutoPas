/**
 * @file FunctorCuda.cuh
 *
 * @date 18.04.2019
 * @author jspahl
 */

#pragma once

#include "cuda_runtime.h"

namespace autopas {

/**
 * Wraps vectors of size 3 with the required precision
 * @tparam T floating point Type
 */
template <typename T>
struct vec3 {
  using Type = T;
};

/**
 * Specialized version for float.
 */
template <>
struct vec3<float> {
  using Type = float3;
};

/**
 * Specialized version for double.
 */
template <>
struct vec3<double> {
  using Type = double3;
};

/**
 * Wraps vectors of size 4 with the required precision
 * @tparam T floating point Type
 */
template <typename T>
struct vec4 {
  using Type = T;
};

/**
 * Specialized version for float.
 */
template <>
struct vec4<float> {
  using Type = float4;
};

/**
 * Specialized version for float.
 */
template <>
struct vec4<double> {
  using Type = double4;
};

/**
 * Class to specify SoAs for cuda.
 * @tparam floatingPointType
 */
template <typename floatingPointType>
class FunctorCudaSoA {};

/**
 * Class to specify constants.
 * @tparam floatingPointType
 */
template <typename floatingPointType>
class FunctorCudaConstants {};

/**
 * Interface for Cuda wrapper classes.
 * @tparam floatingPointType
 */
template <typename floatingPointType>
class CudaWrapperInterface {
 public:
  using floatType = floatingPointType;

  virtual void setNumThreads(int num_threads) = 0;
  virtual void loadConstants(FunctorCudaConstants<floatType> *constants) = 0;

  virtual void SoAFunctorNoN3Wrapper(FunctorCudaSoA<floatType> *cell1, cudaStream_t stream = 0) = 0;
  virtual void SoAFunctorNoN3PairWrapper(FunctorCudaSoA<floatType> *cell1, FunctorCudaSoA<floatType> *cell2,
                                         cudaStream_t stream) = 0;

  virtual void SoAFunctorN3Wrapper(FunctorCudaSoA<floatType> *cell1, cudaStream_t stream = 0) = 0;
  virtual void SoAFunctorN3PairWrapper(FunctorCudaSoA<floatType> *cell1, FunctorCudaSoA<floatType> *cell2,
                                       cudaStream_t stream) = 0;

  virtual void LinkedCellsTraversalNoN3Wrapper(FunctorCudaSoA<floatType> *cell1, unsigned int reqThreads,
                                               unsigned int cids_size, unsigned int *cids, unsigned int cellSizes_size,
                                               size_t *cellSizes, cudaStream_t stream) = 0;

  virtual void LinkedCellsTraversalN3Wrapper(FunctorCudaSoA<floatType> *cell1, unsigned int reqThreads,
                                             unsigned int cids_size, unsigned int *cids, unsigned int cellSizes_size,
                                             size_t *cellSizes, cudaStream_t stream) = 0;

  virtual void CellVerletTraversalNoN3Wrapper(FunctorCudaSoA<floatType> *cell1Base, unsigned int ncells,
                                              unsigned int clusterSize, unsigned int others_size,
                                              unsigned int *other_ids, cudaStream_t stream) = 0;

  virtual void CellVerletTraversalN3Wrapper(FunctorCudaSoA<floatType> *cell1Base, unsigned int ncells,
                                            unsigned int clusterSize, unsigned int others_size, unsigned int *other_ids,
                                            cudaStream_t stream) = 0;

  virtual void loadLinkedCellsOffsets(unsigned int offsets_size, int *offsets) = 0;

  virtual bool isAppropriateClusterSize(unsigned int clusterSize) const = 0;
};

}  // namespace autopas
