/**
 * @file LJFunctorCuda.cu
 *
 * @date 9.1.2019
 * @author jspahl
 */
#include <iostream>
#include "autopas/pairwiseFunctors/LJFunctorCuda.cuh"
#include "autopas/utils/CudaExceptionHandler.h"
#include "autopas/utils/ExceptionHandler.h"
#include "math_constants.h"

namespace autopas {

#define CREATESWITCHCASE(blockSize, gridSize, function, params)                \
  case blockSize:                                                              \
    function<floatType, blockSize><<<gridSize, blockSize, 0, stream>>> params; \
    break;

#define CREATESWITCHCASES(gridSize, function, params) \
  CREATESWITCHCASE(32, gridSize, function, params)    \
  CREATESWITCHCASE(64, gridSize, function, params)    \
  CREATESWITCHCASE(96, gridSize, function, params)    \
  CREATESWITCHCASE(128, gridSize, function, params)   \
  CREATESWITCHCASE(160, gridSize, function, params)   \
  CREATESWITCHCASE(192, gridSize, function, params)   \
  CREATESWITCHCASE(224, gridSize, function, params)   \
  CREATESWITCHCASE(256, gridSize, function, params)   \
  CREATESWITCHCASE(288, gridSize, function, params)   \
  CREATESWITCHCASE(320, gridSize, function, params)   \
  CREATESWITCHCASE(352, gridSize, function, params)   \
  CREATESWITCHCASE(384, gridSize, function, params)   \
  CREATESWITCHCASE(416, gridSize, function, params)   \
  CREATESWITCHCASE(448, gridSize, function, params)   \
  CREATESWITCHCASE(480, gridSize, function, params)   \
  CREATESWITCHCASE(512, gridSize, function, params)   \
  CREATESWITCHCASE(544, gridSize, function, params)   \
  CREATESWITCHCASE(576, gridSize, function, params)   \
  CREATESWITCHCASE(608, gridSize, function, params)   \
  CREATESWITCHCASE(640, gridSize, function, params)   \
  CREATESWITCHCASE(672, gridSize, function, params)   \
  CREATESWITCHCASE(704, gridSize, function, params)   \
  CREATESWITCHCASE(736, gridSize, function, params)   \
  CREATESWITCHCASE(768, gridSize, function, params)   \
  CREATESWITCHCASE(800, gridSize, function, params)   \
  CREATESWITCHCASE(832, gridSize, function, params)   \
  CREATESWITCHCASE(864, gridSize, function, params)   \
  CREATESWITCHCASE(896, gridSize, function, params)   \
  CREATESWITCHCASE(928, gridSize, function, params)   \
  CREATESWITCHCASE(960, gridSize, function, params)   \
  CREATESWITCHCASE(992, gridSize, function, params)   \
  CREATESWITCHCASE(1024, gridSize, function, params)

__constant__ LJFunctorConstants<float> global_constants_float;
__constant__ LJFunctorConstants<double> global_constants_double;

__constant__ unsigned int linkedCellsOffsetsSize;
__constant__ int linkedCellsOffsets[27];

template <typename T>
__device__ inline LJFunctorConstants<T>& getConstants() {
  return global_constants_float;
}
template <>
__device__ inline LJFunctorConstants<double>& getConstants<double>() {
  return global_constants_double;
}

template <typename T>
__device__ inline T getInfinity() {
  return CUDART_INF_F;
}
template <>
__device__ inline double getInfinity<double>() {
  return CUDART_INF;
}

template <typename floatType>
__device__ inline typename vec3<floatType>::Type bodyBodyF(typename vec3<floatType>::Type i,
                                                           typename vec3<floatType>::Type j,
                                                           typename vec3<floatType>::Type fi) {
  floatType drx = i.x - j.x;
  floatType drz = i.z - j.z;
  floatType dry = i.y - j.y;

  floatType dr2 = drx * drx + dry * dry + drz * drz;

  if (dr2 > getConstants<floatType>().cutoffsquare | dr2 == 0.0) {
    return fi;
  }

  floatType invdr2 = 1. / dr2;
  floatType lj6 = getConstants<floatType>().sigmasquare * invdr2;
  lj6 = lj6 * lj6 * lj6;
  floatType lj12 = lj6 * lj6;
  floatType lj12m6 = lj12 - lj6;
  floatType fac = getConstants<floatType>().epsilon24 * (lj12 + lj12m6) * invdr2;

  fi.x += drx * fac;
  fi.y += dry * fac;
  fi.z += drz * fac;

  return fi;
}

template <typename floatType, bool n3AdditionSafe = false>
__device__ inline typename vec3<floatType>::Type bodyBodyFN3(typename vec3<floatType>::Type i,
                                                             typename vec3<floatType>::Type j,
                                                             typename vec3<floatType>::Type fi,
                                                             typename vec3<floatType>::Type* fj) {
  floatType drx = i.x - j.x;
  floatType drz = i.z - j.z;
  floatType dry = i.y - j.y;

  floatType dr2 = drx * drx + dry * dry + drz * drz;

  if (dr2 > getConstants<floatType>().cutoffsquare) {
    return fi;
  }

  floatType invdr2 = 1. / dr2;
  floatType lj6 = getConstants<floatType>().sigmasquare * invdr2;
  lj6 = lj6 * lj6 * lj6;
  floatType lj12 = lj6 * lj6;
  floatType lj12m6 = lj12 - lj6;
  floatType fac = getConstants<floatType>().epsilon24 * (lj12 + lj12m6) * invdr2;

  floatType dfx = drx * fac;
  floatType dfy = dry * fac;
  floatType dfz = drz * fac;

  fi.x += dfx;
  fi.y += dfy;
  fi.z += dfz;

  if (n3AdditionSafe) {
    fj->x -= dfx;
    fj->y -= dfy;
    fj->z -= dfz;
  } else {
    atomicAdd(&(fj->x), -dfx);
    atomicAdd(&(fj->y), -dfy);
    atomicAdd(&(fj->z), -dfz);
  }
  return fi;
}

template <typename floatType, int block_size>
__global__ void SoAFunctorNoN3(LJFunctorCudaSoA<floatType> cell1) {
  __shared__ typename vec3<floatType>::Type block_pos[block_size];
  int i, tile;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  typename vec3<floatType>::Type myposition = {0, 0, 0};
  typename vec3<floatType>::Type myf = {0, 0, 0};
  if (tid < cell1._size) {
    myposition.x = cell1._posX[tid];
    myposition.y = cell1._posY[tid];
    myposition.z = cell1._posZ[tid];
  }

  for (i = block_size, tile = 0; i < cell1._size; i += block_size, ++tile) {
    int idx = tile * block_size + threadIdx.x;

    block_pos[threadIdx.x] = {cell1._posX[idx], cell1._posY[idx], cell1._posZ[idx]};
    __syncthreads();
    if (tid < cell1._size) {
      for (int j = 0; j < blockDim.x; ++j) {
        myf = bodyBodyF<floatType>(myposition, block_pos[j], myf);
      }
    }
    __syncthreads();
  }
  {
    int idx = tile * block_size + threadIdx.x;
    block_pos[threadIdx.x] = {cell1._posX[idx], cell1._posY[idx], cell1._posZ[idx]};
    __syncthreads();

    const int size = cell1._size - tile * blockDim.x;
    for (int j = 0; j < size; ++j) {
      myf = bodyBodyF<floatType>(myposition, block_pos[j], myf);
    }

    __syncthreads();
  }

  atomicAdd(cell1._forceX + tid, myf.x);
  atomicAdd(cell1._forceY + tid, myf.y);
  atomicAdd(cell1._forceZ + tid, myf.z);
}

template <typename floatType, int block_size>
__global__ void SoAFunctorNoN3Pair(LJFunctorCudaSoA<floatType> cell1, LJFunctorCudaSoA<floatType> cell2) {
  __shared__ typename vec3<floatType>::Type block_pos[block_size];
  int i, tile;
  int tid = blockIdx.x * block_size + threadIdx.x;
  typename vec3<floatType>::Type myposition;
  typename vec3<floatType>::Type myf = {0, 0, 0};

  if (tid < cell1._size) {
    myposition.x = cell1._posX[tid];
    myposition.y = cell1._posY[tid];
    myposition.z = cell1._posZ[tid];
  }

  for (i = 0, tile = 0; i < cell2._size; i += block_size, ++tile) {
    int idx = tile * block_size + threadIdx.x;

    if (idx < cell2._size) block_pos[threadIdx.x] = {cell2._posX[idx], cell2._posY[idx], cell2._posZ[idx]};
    __syncthreads();

    const int size = min(block_size, cell2._size - i);
    for (int j = 0; j < size; ++j) {
      myf = bodyBodyF<floatType>(myposition, block_pos[j], myf);
    }
    __syncthreads();
  }
  atomicAdd(cell1._forceX + tid, myf.x);
  atomicAdd(cell1._forceY + tid, myf.y);
  atomicAdd(cell1._forceZ + tid, myf.z);
}

template <typename floatType, int block_size, bool NMisMultipleBlockSize = false>
__global__ void SoAFunctorN3(LJFunctorCudaSoA<floatType> cell1) {
  __shared__ typename vec3<floatType>::Type cell1_pos_shared[block_size];
  __shared__ typename vec3<floatType>::Type cell1_forces_shared[block_size];
  int tid = blockIdx.x * block_size + threadIdx.x;
  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};
  int i, tile;
  const int mask = block_size - 1;

  if (not NMisMultipleBlockSize && tid < cell1._size) {
    myposition.x = cell1._posX[tid];
    myposition.y = cell1._posY[tid];
    myposition.z = cell1._posZ[tid];
  }

  for (i = 0, tile = 0; tile < blockIdx.x; i += block_size, ++tile) {
    int idx = tile * block_size + threadIdx.x;
    cell1_pos_shared[threadIdx.x] = {cell1._posX[idx], cell1._posY[idx], cell1._posZ[idx]};
    cell1_forces_shared[threadIdx.x] = {0, 0, 0};
    __syncthreads();

    for (int j = 0; j < block_size; ++j) {
      unsigned int offset;
      if ((block_size & (block_size - 1)) == 0) {
        offset = (j + threadIdx.x) & mask;
      } else {
        offset = (j + threadIdx.x) % mask;
      }
      myf = bodyBodyFN3<floatType>(myposition, cell1_pos_shared[offset], myf, cell1_forces_shared + offset);
    }
    __syncthreads();

    atomicAdd(cell1._forceX + idx, cell1_forces_shared[threadIdx.x].x);
    atomicAdd(cell1._forceY + idx, cell1_forces_shared[threadIdx.x].y);
    atomicAdd(cell1._forceZ + idx, cell1_forces_shared[threadIdx.x].z);
    __syncthreads();
  }

  {
    int idx = blockIdx.x * block_size + threadIdx.x;
    cell1_pos_shared[threadIdx.x] = {cell1._posX[idx], cell1._posY[idx], cell1._posZ[idx]};
    cell1_forces_shared[threadIdx.x] = {0, 0, 0};
    __syncthreads();

    for (int j = threadIdx.x - 1; j >= 0; --j) {
      myf = bodyBodyFN3<floatType>(myposition, cell1_pos_shared[j], myf, cell1_forces_shared + j);
    }
    __syncthreads();

    atomicAdd(cell1._forceX + idx, cell1_forces_shared[threadIdx.x].x);
    atomicAdd(cell1._forceY + idx, cell1_forces_shared[threadIdx.x].y);
    atomicAdd(cell1._forceZ + idx, cell1_forces_shared[threadIdx.x].z);
    __syncthreads();
  }

  atomicAdd(cell1._forceX + tid, myf.x);
  atomicAdd(cell1._forceY + tid, myf.y);
  atomicAdd(cell1._forceZ + tid, myf.z);
}

template <typename floatType, int block_size, bool NMisMultipleBlockSize = false>
__global__ void SoAFunctorN3Pair(LJFunctorCudaSoA<floatType> cell1, LJFunctorCudaSoA<floatType> cell2) {
  __shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
  __shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];
  int tid = blockIdx.x * block_size + threadIdx.x;
  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};
  int i, tile;
  const int mask = block_size - 1;

  if (not NMisMultipleBlockSize && tid < cell1._size) {
    myposition.x = cell1._posX[tid];
    myposition.y = cell1._posY[tid];
    myposition.z = cell1._posZ[tid];
  }
  for (i = block_size, tile = 0; i <= cell2._size; i += block_size, ++tile) {
    int idx = tile * block_size + threadIdx.x;
    cell2_pos_shared[threadIdx.x] = {cell2._posX[idx], cell2._posY[idx], cell2._posZ[idx]};
    cell2_forces_shared[threadIdx.x] = {0, 0, 0};
    __syncthreads();

    for (int j = 0; j < block_size; ++j) {
      unsigned int offset;
      if ((block_size & (block_size - 1)) == 0) {
        offset = (j + threadIdx.x) & mask;
      } else {
        offset = (j + threadIdx.x) % mask;
      }
      myf = bodyBodyFN3<floatType, true>(myposition, cell2_pos_shared[offset], myf, cell2_forces_shared + offset);
    }
    __syncthreads();

    atomicAdd(cell2._forceX + idx, cell2_forces_shared[threadIdx.x].x);
    atomicAdd(cell2._forceY + idx, cell2_forces_shared[threadIdx.x].y);
    atomicAdd(cell2._forceZ + idx, cell2_forces_shared[threadIdx.x].z);
    __syncthreads();
  }
  if ((not NMisMultipleBlockSize) && (i > cell2._size)) {
    int idx = tile * block_size + threadIdx.x;
    if (idx < cell2._size) {
      cell2_pos_shared[threadIdx.x] = {cell2._posX[idx], cell2._posY[idx], cell2._posZ[idx]};
      cell2_forces_shared[threadIdx.x] = {0, 0, 0};
    }
    __syncthreads();

    const int size = block_size + cell2._size - i;
    for (int j = 0; j < size; ++j) {
      const int offset = (j + threadIdx.x) % size;
      myf = bodyBodyFN3<floatType>(myposition, cell2_pos_shared[offset], myf, cell2_forces_shared + offset);
    }
    __syncthreads();
    if (idx < cell2._size) {
      atomicAdd(cell2._forceX + idx, cell2_forces_shared[threadIdx.x].x);
      atomicAdd(cell2._forceY + idx, cell2_forces_shared[threadIdx.x].y);
      atomicAdd(cell2._forceZ + idx, cell2_forces_shared[threadIdx.x].z);
      __syncthreads();
    }
  }

  atomicAdd(cell1._forceX + tid, myf.x);
  atomicAdd(cell1._forceY + tid, myf.y);
  atomicAdd(cell1._forceZ + tid, myf.z);
}

template <typename floatType>
void LJFunctorCudaWrapper<floatType>::SoAFunctorNoN3Wrapper(FunctorCudaSoA<floatType>* cell1Base, cudaStream_t stream) {
  LJFunctorCudaSoA<floatType> cell1 = *static_cast<LJFunctorCudaSoA<floatType>*>(cell1Base);

  switch (_num_threads) {
    CREATESWITCHCASES(numRequiredBlocks(cell1._size), SoAFunctorNoN3, (cell1));
    default:
      autopas::utils::ExceptionHandler::exception(std::string("cuda Kernel size not available"));
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template <typename floatType>
void LJFunctorCudaWrapper<floatType>::SoAFunctorNoN3PairWrapper(FunctorCudaSoA<floatType>* cell1Base,
                                                                FunctorCudaSoA<floatType>* cell2Base,
                                                                cudaStream_t stream) {
  LJFunctorCudaSoA<floatType> cell1 = *static_cast<LJFunctorCudaSoA<floatType>*>(cell1Base);
  LJFunctorCudaSoA<floatType> cell2 = *static_cast<LJFunctorCudaSoA<floatType>*>(cell2Base);

  switch (_num_threads) {
    CREATESWITCHCASES(numRequiredBlocks(cell1._size), SoAFunctorNoN3Pair, (cell1, cell2));
    default:
      autopas::utils::ExceptionHandler::exception(std::string("cuda Kernel size not available"));
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template <typename floatType>
void LJFunctorCudaWrapper<floatType>::SoAFunctorN3Wrapper(FunctorCudaSoA<floatType>* cell1Base, cudaStream_t stream) {
  LJFunctorCudaSoA<floatType> cell1 = *static_cast<LJFunctorCudaSoA<floatType>*>(cell1Base);

  switch (_num_threads) {
    CREATESWITCHCASES(numRequiredBlocks(cell1._size), SoAFunctorN3, (cell1));
    default:
      autopas::utils::ExceptionHandler::exception(std::string("cuda Kernel size not available"));
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template <typename floatType>
void LJFunctorCudaWrapper<floatType>::SoAFunctorN3PairWrapper(FunctorCudaSoA<floatType>* cell1Base,
                                                              FunctorCudaSoA<floatType>* cell2Base,
                                                              cudaStream_t stream) {
  LJFunctorCudaSoA<floatType> cell1 = *static_cast<LJFunctorCudaSoA<floatType>*>(cell1Base);
  LJFunctorCudaSoA<floatType> cell2 = *static_cast<LJFunctorCudaSoA<floatType>*>(cell2Base);

  switch (_num_threads) {
    CREATESWITCHCASES(numRequiredBlocks(cell1._size), SoAFunctorN3Pair, (cell1, cell2));
    default:
      autopas::utils::ExceptionHandler::exception(std::string("cuda Kernel size not available"));
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template <typename floatType, int block_size>
__global__ void LinkedCellsTraversalNoN3(LJFunctorCudaSoA<floatType> cell, unsigned int* cids, size_t* cellSizes) {
  int own_cid = cids[blockIdx.x];
  __shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};

  int index = cellSizes[own_cid] + threadIdx.x;
  if (threadIdx.x < (cellSizes[own_cid + 1] - cellSizes[own_cid])) {
    myposition.x = cell._posX[index];
    myposition.y = cell._posY[index];
    myposition.z = cell._posZ[index];
  }

  // other cells
  for (auto other_index = 0; other_index < linkedCellsOffsetsSize; ++other_index) {
    const int other_id = own_cid + linkedCellsOffsets[other_index];
    const size_t cell2Start = cellSizes[other_id];
    const unsigned int sizeCell2 = cellSizes[other_id + 1] - cell2Start;

    cell2_pos_shared[threadIdx.x] = {cell._posX[cell2Start + threadIdx.x], cell._posY[cell2Start + threadIdx.x],
                                     cell._posZ[cell2Start + threadIdx.x]};
    __syncthreads();
    for (int j = 0; j < sizeCell2; ++j) {
      myf = bodyBodyF<floatType>(myposition, cell2_pos_shared[j], myf);
    }
    __syncthreads();
  }
  if (threadIdx.x < (cellSizes[own_cid + 1] - cellSizes[own_cid])) {
    atomicAdd(cell._forceX + index, myf.x);
    atomicAdd(cell._forceY + index, myf.y);
    atomicAdd(cell._forceZ + index, myf.z);
  }
}

template <typename floatType, int block_size>
__global__ void LinkedCellsTraversalN3(LJFunctorCudaSoA<floatType> cell, unsigned int* cids, size_t* cellSizes) {
  int own_cid = cids[blockIdx.x];
  __shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
  __shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];

  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};

  int index = cellSizes[own_cid] + threadIdx.x;
  if (threadIdx.x < (cellSizes[own_cid + 1] - cellSizes[own_cid])) {
    myposition.x = cell._posX[index];
    myposition.y = cell._posY[index];
    myposition.z = cell._posZ[index];
  }
  // other cells
  for (auto other_index = 0; other_index < linkedCellsOffsetsSize; ++other_index) {
    const int other_id = own_cid + linkedCellsOffsets[other_index];
    const size_t cell2Start = cellSizes[other_id];
    const int sizeCell2 = cellSizes[other_id + 1] - cell2Start;

    cell2_pos_shared[threadIdx.x] = {cell._posX[cell2Start + threadIdx.x], cell._posY[cell2Start + threadIdx.x],
                                     cell._posZ[cell2Start + threadIdx.x]};
    cell2_forces_shared[threadIdx.x] = {0, 0, 0};
    __syncthreads();
    for (int j = 0; j < sizeCell2; ++j) {
      const int offset = (j + threadIdx.x) % sizeCell2;
      myf = bodyBodyFN3<floatType, false>(myposition, cell2_pos_shared[offset], myf, cell2_forces_shared + offset);
    }
    __syncthreads();

    atomicAdd(cell._forceX + cell2Start + threadIdx.x, cell2_forces_shared[threadIdx.x].x);
    atomicAdd(cell._forceY + cell2Start + threadIdx.x, cell2_forces_shared[threadIdx.x].y);
    atomicAdd(cell._forceZ + cell2Start + threadIdx.x, cell2_forces_shared[threadIdx.x].z);
    __syncthreads();
  }
  // same cells
  {
    const int cell1Start = cellSizes[own_cid];
    const int sizeCell1 = cellSizes[own_cid + 1] - cell1Start;

    cell2_pos_shared[threadIdx.x] = {cell._posX[cell1Start + threadIdx.x], cell._posY[cell1Start + threadIdx.x],
                                     cell._posZ[cell1Start + threadIdx.x]};
    __syncthreads();
    for (int j = 0; j < sizeCell1; ++j) {
      myf = bodyBodyF<floatType>(myposition, cell2_pos_shared[j], myf);
    }
    __syncthreads();
  }
  if (threadIdx.x < (cellSizes[own_cid + 1] - cellSizes[own_cid])) {
    atomicAdd(cell._forceX + index, myf.x);
    atomicAdd(cell._forceY + index, myf.y);
    atomicAdd(cell._forceZ + index, myf.z);
  }
}

template <typename floatType>
void LJFunctorCudaWrapper<floatType>::LinkedCellsTraversalNoN3Wrapper(FunctorCudaSoA<floatType>* cell1Base,
                                                                      unsigned int reqThreads, unsigned int cids_size,
                                                                      unsigned int* cids, unsigned int cellSizes_size,
                                                                      size_t* cellSizes, cudaStream_t stream) {
  LJFunctorCudaSoA<floatType> cell1 = *static_cast<LJFunctorCudaSoA<floatType>*>(cell1Base);

  switch (reqThreads) {
    CREATESWITCHCASES(cids_size, LinkedCellsTraversalNoN3, (cell1, cids, cellSizes));
    default:
      autopas::utils::ExceptionHandler::exception(
          "Linked Cells NoN3: cuda Kernel size not available for Linked cells available. Too many particles "
          "in a cell. Requested: {}",
          reqThreads);
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template <typename floatType>
void LJFunctorCudaWrapper<floatType>::LinkedCellsTraversalN3Wrapper(FunctorCudaSoA<floatType>* cell1Base,
                                                                    unsigned int reqThreads, unsigned int cids_size,
                                                                    unsigned int* cids, unsigned int cellSizes_size,
                                                                    size_t* cellSizes, cudaStream_t stream) {
  LJFunctorCudaSoA<floatType> cell1 = *static_cast<LJFunctorCudaSoA<floatType>*>(cell1Base);

  switch (reqThreads) {
    CREATESWITCHCASES(cids_size, LinkedCellsTraversalN3, (cell1, cids, cellSizes));
    default:
      autopas::utils::ExceptionHandler::exception(
          "Linked Cells N3:cuda Kernel size not available for Linked cells available. Too many particles in "
          "a cell. Requested: {}",
          reqThreads);
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template <typename floatType, int block_size>
__global__ void CellVerletTraversalNoN3(LJFunctorCudaSoA<floatType> cell, const unsigned int others_size,
                                        unsigned int* other_ids) {
  __shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};

  unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
  myposition.x = cell._posX[index];
  myposition.y = cell._posY[index];
  myposition.z = cell._posZ[index];

  // other cells
  for (auto other_index = others_size * blockIdx.x; other_ids[other_index] < UINT_MAX; ++other_index) {
    const size_t cell2Start = blockDim.x * other_ids[other_index];

    const size_t own_particle = cell2Start + threadIdx.x;
    cell2_pos_shared[threadIdx.x] = {cell._posX[own_particle], cell._posY[own_particle], cell._posZ[own_particle]};
    __syncthreads();
    for (int j = 0; j < block_size; ++j) {
      myf = bodyBodyF<floatType>(myposition, cell2_pos_shared[j], myf);
    }
    __syncthreads();
  }
  atomicAdd(cell._forceX + index, myf.x);
  atomicAdd(cell._forceY + index, myf.y);
  atomicAdd(cell._forceZ + index, myf.z);
}

template <typename floatType>
void LJFunctorCudaWrapper<floatType>::CellVerletTraversalNoN3Wrapper(FunctorCudaSoA<floatType>* cell1Base,
                                                                     unsigned int ncells, unsigned int clusterSize,
                                                                     unsigned int others_size, unsigned int* other_ids,
                                                                     cudaStream_t stream) {
  LJFunctorCudaSoA<floatType> cell1 = *static_cast<LJFunctorCudaSoA<floatType>*>(cell1Base);
  switch (clusterSize) {
    CREATESWITCHCASES(ncells, CellVerletTraversalNoN3, (cell1, others_size, other_ids));
    default:
      autopas::utils::ExceptionHandler::exception(
          "cuda Kernel size not available for Verlet cells available. Too many particles in a cell. "
          "Requested: {}",
          clusterSize);
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template <typename floatType, int block_size>
__global__ void CellVerletTraversalN3(LJFunctorCudaSoA<floatType> cell, unsigned int others_size,
                                      unsigned int* other_ids) {
  const unsigned int mask = block_size - 1;

  __shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
  __shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];

  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};

  int index = blockIdx.x * blockDim.x + threadIdx.x;
  myposition.x = cell._posX[index];
  myposition.y = cell._posY[index];
  myposition.z = cell._posZ[index];

  // other cells
  for (auto other_index = others_size * blockIdx.x; other_ids[other_index] != UINT_MAX; ++other_index) {
    const unsigned int cell2Start = blockDim.x * other_ids[other_index];

    cell2_pos_shared[threadIdx.x] = {cell._posX[cell2Start + threadIdx.x], cell._posY[cell2Start + threadIdx.x],
                                     cell._posZ[cell2Start + threadIdx.x]};
    cell2_forces_shared[threadIdx.x] = {0, 0, 0};
    __syncthreads();
    for (int j = 0; j < block_size; ++j) {
      unsigned int offset = 0;
      if ((block_size & (block_size - 1)) == 0) {
        offset = (j + threadIdx.x) & mask;
      } else {
        offset = (j + threadIdx.x) % block_size;
      }
      myf = bodyBodyFN3<floatType, false>(myposition, cell2_pos_shared[offset], myf, cell2_forces_shared + offset);
    }
    __syncthreads();

    atomicAdd(cell._forceX + cell2Start + threadIdx.x, cell2_forces_shared[threadIdx.x].x);
    atomicAdd(cell._forceY + cell2Start + threadIdx.x, cell2_forces_shared[threadIdx.x].y);
    atomicAdd(cell._forceZ + cell2Start + threadIdx.x, cell2_forces_shared[threadIdx.x].z);
    __syncthreads();
  }

  // same cluster without N3
  {
    const unsigned int cellStart = blockIdx.x * blockDim.x;

    cell2_pos_shared[threadIdx.x] = {cell._posX[cellStart + threadIdx.x], cell._posY[cellStart + threadIdx.x],
                                     cell._posZ[cellStart + threadIdx.x]};
    __syncthreads();
    for (int j = 0; j < block_size; ++j) {
      myf = bodyBodyF<floatType>(myposition, cell2_pos_shared[j], myf);
    }
    __syncthreads();
  }
  atomicAdd(cell._forceX + index, myf.x);
  atomicAdd(cell._forceY + index, myf.y);
  atomicAdd(cell._forceZ + index, myf.z);
}

template <typename floatType>
void LJFunctorCudaWrapper<floatType>::CellVerletTraversalN3Wrapper(FunctorCudaSoA<floatType>* cell1Base,
                                                                   unsigned int ncells, unsigned int clusterSize,
                                                                   unsigned int others_size, unsigned int* other_ids,
                                                                   cudaStream_t stream) {
  LJFunctorCudaSoA<floatType> cell1 = *static_cast<LJFunctorCudaSoA<floatType>*>(cell1Base);
  switch (clusterSize) {
    CREATESWITCHCASES(ncells, CellVerletTraversalN3, (cell1, others_size, other_ids));
    default:
      autopas::utils::ExceptionHandler::exception(
          "cuda Kernel size not available for Verlet cells available. Too many particles in a cell. "
          "Requested: {}",
          clusterSize);
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template <>
void LJFunctorCudaWrapper<float>::loadConstants(FunctorCudaConstants<float>* constants) {
  LJFunctorConstants<float>* c = static_cast<LJFunctorConstants<float>*>(constants);

  cudaMemcpyToSymbol(global_constants_float, c, sizeof(LJFunctorConstants<float>));
}

template <>
void LJFunctorCudaWrapper<double>::loadConstants(FunctorCudaConstants<double>* constants) {
  LJFunctorConstants<double>* c = static_cast<LJFunctorConstants<double>*>(constants);
  autopas::utils::CudaExceptionHandler::checkErrorCode(
      cudaMemcpyToSymbol(global_constants_double, c, sizeof(LJFunctorConstants<double>)));
}

template <typename T>
void LJFunctorCudaWrapper<T>::loadConstants(FunctorCudaConstants<T>* constants) {
  autopas::utils::ExceptionHandler::exception("Cuda constants with unknown Type loaded");
}

template <typename T>
void LJFunctorCudaWrapper<T>::loadLinkedCellsOffsets(unsigned int offsets_size, int* offsets) {
  if (offsets_size > 27) {
    autopas::utils::ExceptionHandler::exception(
        "LJFunctorCudaWrapper does not support linked cells with >27 neighbors");
  }
  autopas::utils::CudaExceptionHandler::checkErrorCode(
      cudaMemcpyToSymbol(linkedCellsOffsetsSize, &offsets_size, sizeof(unsigned int)));
  autopas::utils::CudaExceptionHandler::checkErrorCode(
      cudaMemcpyToSymbol(linkedCellsOffsets, offsets, offsets_size * sizeof(int)));
}

template class LJFunctorCudaWrapper<float>;
template class LJFunctorCudaWrapper<double>;

}  // namespace autopas
