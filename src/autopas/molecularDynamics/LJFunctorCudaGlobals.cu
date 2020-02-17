/**
 * @file LJFunctorCuda.cu
 *
 * @date 26.4.2019
 * @author jspahl
 */
#include <iostream>

#include "LJFunctorCudaGlobals.cuh"
#include "autopas/molecularDynamics/LJFunctorCuda.cuh"
#include "autopas/utils/CudaExceptionHandler.h"
#include "autopas/utils/ExceptionHandler.h"
#include "math_constants.h"

namespace autopas {

/**
 * Macro to create a switch case for starting the kernel function with the given block and grid size as well as
 * parameters
 */
#define CREATESWITCHCASE(blockSize, gridSize, sharedMemSize, function, params)             \
  case blockSize:                                                                          \
    function<floatType, blockSize><<<gridSize, blockSize, sharedMemSize, stream>>> params; \
    break;

/**
 * Macro to create cases for all possible block dimensions
 */
#define CREATESWITCHCASES(gridSize, sharedMemSize, function, params) \
  CREATESWITCHCASE(32, gridSize, sharedMemSize, function, params)    \
  CREATESWITCHCASE(64, gridSize, sharedMemSize, function, params)    \
  CREATESWITCHCASE(96, gridSize, sharedMemSize, function, params)    \
  CREATESWITCHCASE(128, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(160, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(192, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(224, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(256, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(288, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(320, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(352, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(384, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(416, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(448, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(480, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(512, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(544, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(576, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(608, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(640, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(672, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(704, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(736, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(768, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(800, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(832, gridSize, sharedMemSize, function, params)   \
  CREATESWITCHCASE(864, gridSize, sharedMemSize, function, params)

/**
 * global constant with float precision
 */
__constant__ LJFunctorConstants<float> global_constants_float;
/**
 * global constant with double precision
 */
__constant__ LJFunctorConstants<double> global_constants_double;

/**
 * number of offsets in linkedCellsOffsets array
 */
__constant__ unsigned int linkedCellsOffsetsSize;
/**
 * offsets to neighbor cells when using linked cells
 */
__constant__ int linkedCellsOffsets[27];

/**
 * Get constant storage depending on template parameter
 * @return constants
 */
template <typename T>
__device__ inline LJFunctorConstants<T> &getConstants() {
  return global_constants_float;
}
template <>
__device__ inline LJFunctorConstants<double> &getConstants<double>() {
  return global_constants_double;
}

/**
 * returns Infinity in the requested precision
 * @return infinity
 */
template <typename T>
__device__ inline T getInfinity() {
  return CUDART_INF_F;
}
template <>
__device__ inline double getInfinity<double>() {
  return CUDART_INF;
}

/**
 * Calculates the next power of two higher than or equal to the given number.
 * This function is taken from: http://locklessinc.com/articles/next_pow2/
 * @param x
 * @return
 */
__device__ constexpr uint32_t nextPowerOfTwo(uint32_t x) {
  x -= 1;
  x |= (x >> 1);
  x |= (x >> 2);
  x |= (x >> 4);
  x |= (x >> 8);
  x |= (x >> 16);

  return x + 1;
}

template <typename floatType, int blockSize>
__device__ inline void reduceGlobalsShared(typename vec4<floatType>::Type *globals_shared) {
  constexpr uint32_t nextPowerOfTwoVal = nextPowerOfTwo(blockSize);
  for (uint32_t s = nextPowerOfTwoVal / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s and threadIdx.x + s < blockSize) {
      globals_shared[threadIdx.x].x += globals_shared[threadIdx.x + s].x;
      globals_shared[threadIdx.x].y += globals_shared[threadIdx.x + s].y;
      globals_shared[threadIdx.x].z += globals_shared[threadIdx.x + s].z;
      globals_shared[threadIdx.x].w += globals_shared[threadIdx.x + s].w;
    }
    __syncthreads();
  }
}

/**
 * Calculates the lennard-jones interactions between two particles
 * @tparam floatType flaot number type used in this function
 * @tparam halfGlobals true iff globals should be halfed on addition instead of post processing
 * @param i position of particle i
 * @param j position of particle j
 * @param fi force of particle i
 * @return updated force of particle i
 */
template <typename floatType, bool halfGlobals = false>
__device__ inline typename vec3<floatType>::Type bodyBodyF(typename vec3<floatType>::Type i,
                                                           typename vec4<floatType>::Type j,
                                                           typename vec3<floatType>::Type fi,
                                                           typename vec4<floatType>::Type &globalsi, bool iIsOwn) {
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

  floatType fix = drx * fac;
  floatType fiy = dry * fac;
  floatType fiz = drz * fac;

  if (iIsOwn) {
    if (halfGlobals) {
      globalsi.x += drx * fix * 0.5;
      globalsi.y += dry * fiy * 0.5;
      globalsi.z += drz * fiz * 0.5;
      globalsi.w += 0.5 * (getConstants<floatType>().epsilon24 * lj12m6 + getConstants<floatType>().shift6);
    } else {
      globalsi.x += drx * fix;
      globalsi.y += dry * fiy;
      globalsi.z += drz * fiz;
      globalsi.w += getConstants<floatType>().epsilon24 * lj12m6 + getConstants<floatType>().shift6;
    }
  }

  fi.x += fix;
  fi.y += fiy;
  fi.z += fiz;

  return fi;
}

/**
 * Calculates the lennard-jones interactions between two particles with newton3
 * @tparam floatType flaot number type used in this function
 * @tparam n3AdditionSafe true iff addition to the other particle is threadsafe
 * @param i position of particle i
 * @param j position of particle j
 * @param fi force of particle i
 * @param fj [in/out] force of particle j to update
 * @return updated force of particle i
 */
template <typename floatType, bool n3AdditionSafe = false>
__device__ inline typename vec3<floatType>::Type bodyBodyFN3(
    typename vec3<floatType>::Type i, typename vec4<floatType>::Type j, typename vec3<floatType>::Type fi,
    typename vec3<floatType>::Type *fj, typename vec4<floatType>::Type &globalsi, floatType iOwnMask) {
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

  iOwnMask += j.w;

  globalsi.x += drx * dfx * iOwnMask;
  globalsi.y += dry * dfy * iOwnMask;
  globalsi.z += drz * dfz * iOwnMask;
  globalsi.w += (getConstants<floatType>().epsilon24 * lj12m6 + getConstants<floatType>().shift6) * iOwnMask;

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

/**
 * Calculates all interactions within a single cell
 * @param cell1 particle storage
 */
template <typename floatType, int block_size>
__global__ void SoAFunctorNoN3(LJFunctorCudaGlobalsSoA<floatType> cell1) {
  __shared__ typename vec4<floatType>::Type block_pos[block_size];

  int i, tile;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};
  bool isOwned = false;
  if (tid < cell1._size) {
    myposition.x = cell1._posX[tid];
    myposition.y = cell1._posY[tid];
    myposition.z = cell1._posZ[tid];
    isOwned = cell1._owned[tid] == (floatType)1.0;
  }

  typename vec4<floatType>::Type myglobals = {0, 0, 0, 0};

  for (i = block_size, tile = 0; i < cell1._size; i += block_size, ++tile) {
    int idx = tile * block_size + threadIdx.x;

    block_pos[threadIdx.x] = {cell1._posX[idx], cell1._posY[idx], cell1._posZ[idx]};
    __syncthreads();

    if (tid < cell1._size) {
      for (int j = 0; j < block_size; ++j) {
        myf = bodyBodyF<floatType>(myposition, block_pos[j], myf, myglobals, isOwned);
      }
    }
    __syncthreads();
  }
  {
    int idx = tile * block_size + threadIdx.x;
    block_pos[threadIdx.x] = {cell1._posX[idx], cell1._posY[idx], cell1._posZ[idx]};

    __syncthreads();

    const int size = cell1._size - tile * block_size;
    for (int j = 0; j < size; ++j) {
      myf = bodyBodyF<floatType>(myposition, block_pos[j], myf, myglobals, isOwned);
    }

    __syncthreads();
  }

  // reduce globals
  block_pos[threadIdx.x] = myglobals;
  __syncthreads();

  reduceGlobalsShared<floatType, block_size>(block_pos);
  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(cell1._globals, block_pos[0].x);
    atomicAdd(cell1._globals + 1, block_pos[0].y);
    atomicAdd(cell1._globals + 2, block_pos[0].z);
    atomicAdd(cell1._globals + 3, block_pos[0].w);
  }
  __syncthreads();

  atomicAdd(cell1._forceX + tid, myf.x);
  atomicAdd(cell1._forceY + tid, myf.y);
  atomicAdd(cell1._forceZ + tid, myf.z);
}

/**
 * Calculates all interactions between two cells
 * @tparam floatType
 * @tparam block_size cuda block size
 * @param cell1 first particle storage
 * @paran cell2 second particle storage
 */
template <typename floatType, int block_size>
__global__ void SoAFunctorNoN3Pair(LJFunctorCudaGlobalsSoA<floatType> cell1, LJFunctorCudaGlobalsSoA<floatType> cell2) {
  __shared__ typename vec4<floatType>::Type block_pos[block_size];
  int i, tile;
  int tid = blockIdx.x * block_size + threadIdx.x;
  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};
  typename vec4<floatType>::Type myglobals = {0, 0, 0, 0};
  bool isOwned = false;

  if (tid < cell1._size) {
    myposition.x = cell1._posX[tid];
    myposition.y = cell1._posY[tid];
    myposition.z = cell1._posZ[tid];
    isOwned = cell1._owned[tid] == (floatType)1.0;
  }

  for (i = 0, tile = 0; i < cell2._size; i += block_size, ++tile) {
    int idx = tile * block_size + threadIdx.x;

    if (idx < cell2._size) block_pos[threadIdx.x] = {cell2._posX[idx], cell2._posY[idx], cell2._posZ[idx]};
    __syncthreads();

    const int size = min(block_size, cell2._size - i);
    for (int j = 0; j < size; ++j) {
      myf = bodyBodyF<floatType>(myposition, block_pos[j], myf, myglobals, isOwned);
    }
    __syncthreads();
  }

  // reduce globals
  block_pos[threadIdx.x] = myglobals;
  __syncthreads();

  reduceGlobalsShared<floatType, block_size>(block_pos);
  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(cell1._globals, block_pos[0].x);
    atomicAdd(cell1._globals + 1, block_pos[0].y);
    atomicAdd(cell1._globals + 2, block_pos[0].z);
    atomicAdd(cell1._globals + 3, block_pos[0].w);
  }
  __syncthreads();

  atomicAdd(cell1._forceX + tid, myf.x);
  atomicAdd(cell1._forceY + tid, myf.y);
  atomicAdd(cell1._forceZ + tid, myf.z);
}

/**
 * Calculates all interactions within a single cell using newton 3
 * @tparam floatType
 * @tparam block_size cuda block size
 * @param cell1 particle data
 */
template <typename floatType, int block_size, bool NMisMultipleBlockSize = false>
__global__ void SoAFunctorN3(LJFunctorCudaGlobalsSoA<floatType> cell1) {
  __shared__ typename vec4<floatType>::Type cell1_pos_shared[block_size];
  __shared__ typename vec3<floatType>::Type cell1_forces_shared[block_size];

  int tid = blockIdx.x * block_size + threadIdx.x;
  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};
  typename vec4<floatType>::Type myglobals = {0, 0, 0, 0};
  floatType isOwned = 0.0;

  int i, tile;
  const int mask = block_size - 1;

  if (not NMisMultipleBlockSize && tid < cell1._size) {
    myposition.x = cell1._posX[tid];
    myposition.y = cell1._posY[tid];
    myposition.z = cell1._posZ[tid];
    isOwned = cell1._owned[tid] * (floatType)0.5;
  }

  for (i = 0, tile = 0; tile < blockIdx.x; i += block_size, ++tile) {
    int idx = tile * block_size + threadIdx.x;
    cell1_pos_shared[threadIdx.x] = {cell1._posX[idx], cell1._posY[idx], cell1._posZ[idx],
                                     cell1._owned[idx] * (floatType)0.5};
    cell1_forces_shared[threadIdx.x] = {0, 0, 0};
    __syncthreads();

    for (int j = 0; j < block_size; ++j) {
      unsigned int offset;
      // use bitwise and if equivalent to modulo (for block_size = 2^n)
      if ((block_size & (block_size - 1)) == 0) {
        offset = (j + threadIdx.x) & mask;
      } else {
        offset = (j + threadIdx.x) % block_size;
      }
      myf = bodyBodyFN3<floatType>(myposition, cell1_pos_shared[offset], myf, cell1_forces_shared + offset, myglobals,
                                   isOwned);
    }
    __syncthreads();

    atomicAdd(cell1._forceX + idx, cell1_forces_shared[threadIdx.x].x);
    atomicAdd(cell1._forceY + idx, cell1_forces_shared[threadIdx.x].y);
    atomicAdd(cell1._forceZ + idx, cell1_forces_shared[threadIdx.x].z);
    __syncthreads();
  }

  {
    int idx = blockIdx.x * block_size + threadIdx.x;
    cell1_pos_shared[threadIdx.x] = {cell1._posX[idx], cell1._posY[idx], cell1._posZ[idx],
                                     cell1._owned[idx] * (floatType)0.5};
    cell1_forces_shared[threadIdx.x] = {0, 0, 0};
    __syncthreads();

    for (int j = threadIdx.x - 1; j >= 0; --j) {
      myf = bodyBodyFN3<floatType>(myposition, cell1_pos_shared[j], myf, cell1_forces_shared + j, myglobals, isOwned);
    }
    __syncthreads();

    atomicAdd(cell1._forceX + idx, cell1_forces_shared[threadIdx.x].x);
    atomicAdd(cell1._forceY + idx, cell1_forces_shared[threadIdx.x].y);
    atomicAdd(cell1._forceZ + idx, cell1_forces_shared[threadIdx.x].z);
    __syncthreads();
  }

  // reduce globals
  cell1_pos_shared[threadIdx.x] = myglobals;
  __syncthreads();

  reduceGlobalsShared<floatType, block_size>(cell1_pos_shared);
  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(cell1._globals, cell1_pos_shared[0].x);
    atomicAdd(cell1._globals + 1, cell1_pos_shared[0].y);
    atomicAdd(cell1._globals + 2, cell1_pos_shared[0].z);
    atomicAdd(cell1._globals + 3, cell1_pos_shared[0].w);
  }

  atomicAdd(cell1._forceX + tid, myf.x);
  atomicAdd(cell1._forceY + tid, myf.y);
  atomicAdd(cell1._forceZ + tid, myf.z);
}

/**
 * Calculates all interactions between two cells using newton3
 * @tparam floatType
 * @tparam block_size cuda block size
 * @param cell1 first particle storage
 * @paran cell2 second particle storage
 */
template <typename floatType, int block_size, bool NMisMultipleBlockSize = false>
__global__ void SoAFunctorN3Pair(LJFunctorCudaGlobalsSoA<floatType> cell1, LJFunctorCudaGlobalsSoA<floatType> cell2) {
  __shared__ typename vec4<floatType>::Type cell2_pos_shared[block_size];
  __shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];

  int tid = blockIdx.x * block_size + threadIdx.x;
  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};
  typename vec4<floatType>::Type myglobals = {0, 0, 0, 0};

  int i, tile;
  const int mask = block_size - 1;

  floatType isOwned;
  if (not NMisMultipleBlockSize && tid < cell1._size) {
    myposition.x = cell1._posX[tid];
    myposition.y = cell1._posY[tid];
    myposition.z = cell1._posZ[tid];
    isOwned = cell1._owned[tid] * (floatType)0.5;
  }
  for (i = block_size, tile = 0; i <= cell2._size; i += block_size, ++tile) {
    int idx = tile * block_size + threadIdx.x;
    cell2_pos_shared[threadIdx.x] = {cell2._posX[idx], cell2._posY[idx], cell2._posZ[idx],
                                     cell2._owned[idx] * (floatType)0.5};
    cell2_forces_shared[threadIdx.x] = {0, 0, 0};
    __syncthreads();

    for (int j = 0; j < block_size; ++j) {
      unsigned int offset;
      // use bitwise and if equivalent to modulo (for block_size = 2^n)
      if ((block_size & (block_size - 1)) == 0) {
        offset = (j + threadIdx.x) & mask;
      } else {
        offset = (j + threadIdx.x) % block_size;
      }
      myf = bodyBodyFN3<floatType, true>(myposition, cell2_pos_shared[offset], myf, cell2_forces_shared + offset,
                                         myglobals, isOwned);
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
      cell2_pos_shared[threadIdx.x] = {cell2._posX[idx], cell2._posY[idx], cell2._posZ[idx],
                                       cell2._owned[idx] * (floatType)0.5};
      cell2_forces_shared[threadIdx.x] = {0, 0, 0};
    }
    __syncthreads();

    const int size = block_size + cell2._size - i;
    for (int j = 0; j < size; ++j) {
      const int offset = (j + threadIdx.x) % size;
      myf = bodyBodyFN3<floatType>(myposition, cell2_pos_shared[offset], myf, cell2_forces_shared + offset, myglobals,
                                   isOwned);
    }
    __syncthreads();
    if (idx < cell2._size) {
      atomicAdd(cell2._forceX + idx, cell2_forces_shared[threadIdx.x].x);
      atomicAdd(cell2._forceY + idx, cell2_forces_shared[threadIdx.x].y);
      atomicAdd(cell2._forceZ + idx, cell2_forces_shared[threadIdx.x].z);
    }
    __syncthreads();
  }

  // reduce globals
  cell2_pos_shared[threadIdx.x] = myglobals;
  __syncthreads();

  reduceGlobalsShared<floatType, block_size>(cell2_pos_shared);
  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(cell1._globals, cell2_pos_shared[0].x);
    atomicAdd(cell1._globals + 1, cell2_pos_shared[0].y);
    atomicAdd(cell1._globals + 2, cell2_pos_shared[0].z);
    atomicAdd(cell1._globals + 3, cell2_pos_shared[0].w);
  }

  atomicAdd(cell1._forceX + tid, myf.x);
  atomicAdd(cell1._forceY + tid, myf.y);
  atomicAdd(cell1._forceZ + tid, myf.z);
}

/**
 * Calls SoAFunctorNoN3 with the correct block size
 * @tparam floatType
 * @param cell1Base particle storage
 * @param stream cuda stream to start kernel on
 */
template <typename floatType>
void LJFunctorCudaGlobalsWrapper<floatType>::SoAFunctorNoN3Wrapper(FunctorCudaSoA<floatType> *cell1Base,
                                                                   cudaStream_t stream) {
  LJFunctorCudaGlobalsSoA<floatType> cell1 = *static_cast<LJFunctorCudaGlobalsSoA<floatType> *>(cell1Base);

  switch (_num_threads) {
    CREATESWITCHCASES(numRequiredBlocks(cell1._size), 0, SoAFunctorNoN3, (cell1));
    default:
      autopas::utils::ExceptionHandler::exception(std::string("cuda Kernel size not available"));
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

/**
 * Calls SoAFunctorNoN3Pair with the correct block size
 * @tparam floatType
 * @param cell1Base first particle storage
 * @param cell2Base  second particle storage
 * @param stream cuda stream to start kernel on
 */
template <typename floatType>
void LJFunctorCudaGlobalsWrapper<floatType>::SoAFunctorNoN3PairWrapper(FunctorCudaSoA<floatType> *cell1Base,
                                                                       FunctorCudaSoA<floatType> *cell2Base,
                                                                       cudaStream_t stream) {
  LJFunctorCudaGlobalsSoA<floatType> cell1 = *static_cast<LJFunctorCudaGlobalsSoA<floatType> *>(cell1Base);
  LJFunctorCudaGlobalsSoA<floatType> cell2 = *static_cast<LJFunctorCudaGlobalsSoA<floatType> *>(cell2Base);

  switch (_num_threads) {
    CREATESWITCHCASES(numRequiredBlocks(cell1._size), 0, SoAFunctorNoN3Pair, (cell1, cell2));
    default:
      autopas::utils::ExceptionHandler::exception(std::string("cuda Kernel size not available"));
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

/**
 * Calls SoAFunctorN3 with the correct block size
 * @tparam floatType
 * @param cell1Base particle storage
 * @param stream cuda stream to start kernel on
 */
template <typename floatType>
void LJFunctorCudaGlobalsWrapper<floatType>::SoAFunctorN3Wrapper(FunctorCudaSoA<floatType> *cell1Base,
                                                                 cudaStream_t stream) {
  LJFunctorCudaGlobalsSoA<floatType> cell1 = *static_cast<LJFunctorCudaGlobalsSoA<floatType> *>(cell1Base);

  switch (_num_threads) {
    CREATESWITCHCASES(numRequiredBlocks(cell1._size), 0, SoAFunctorN3, (cell1));
    default:
      autopas::utils::ExceptionHandler::exception(std::string("cuda Kernel size not available"));
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

/**
 * Calls SoAFunctorN3Pair with the correct block size
 * @tparam floatType
 * @param cell1Base first particle storage
 * @param cell2Base  second particle storage
 * @param stream cuda stream to start kernel on
 */
template <typename floatType>
void LJFunctorCudaGlobalsWrapper<floatType>::SoAFunctorN3PairWrapper(FunctorCudaSoA<floatType> *cell1Base,
                                                                     FunctorCudaSoA<floatType> *cell2Base,
                                                                     cudaStream_t stream) {
  LJFunctorCudaGlobalsSoA<floatType> cell1 = *static_cast<LJFunctorCudaGlobalsSoA<floatType> *>(cell1Base);
  LJFunctorCudaGlobalsSoA<floatType> cell2 = *static_cast<LJFunctorCudaGlobalsSoA<floatType> *>(cell2Base);

  switch (_num_threads) {
    CREATESWITCHCASES(numRequiredBlocks(cell1._size), 0, SoAFunctorN3Pair, (cell1, cell2));
    default:
      autopas::utils::ExceptionHandler::exception(std::string("cuda Kernel size not available"));
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

/**
 * Calculates interactions using linked cells
 * @tparam floatType
 * @tparam block_size cuda block size
 * @param cell particle storage
 * @param cids cell ids to calculate interactions for
 * @param cellSizes sizes of the cells by id
 */
template <typename floatType, int block_size>
__global__ void LinkedCellsTraversalNoN3(LJFunctorCudaGlobalsSoA<floatType> cell, unsigned int *cids,
                                         size_t *cellSizes) {
  unsigned int own_cid = cids[blockIdx.x];
  __shared__ typename vec4<floatType>::Type cell2_pos_shared[block_size];
  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};
  typename vec4<floatType>::Type myglobals = {0, 0, 0, 0};
  bool isOwned = false;

  int index = cellSizes[own_cid] + threadIdx.x;
  if (threadIdx.x < (cellSizes[own_cid + 1] - cellSizes[own_cid])) {
    myposition.x = cell._posX[index];
    myposition.y = cell._posY[index];
    myposition.z = cell._posZ[index];
    isOwned = cell._owned[index] == (floatType)1.0;
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
      myf = bodyBodyF<floatType>(myposition, cell2_pos_shared[j], myf, myglobals, isOwned);
    }
    __syncthreads();
  }
  if (threadIdx.x < (cellSizes[own_cid + 1] - cellSizes[own_cid])) {
    atomicAdd(cell._forceX + index, myf.x);
    atomicAdd(cell._forceY + index, myf.y);
    atomicAdd(cell._forceZ + index, myf.z);
  }

  // reduce globals
  cell2_pos_shared[threadIdx.x] = myglobals;
  __syncthreads();

  reduceGlobalsShared<floatType, block_size>(cell2_pos_shared);
  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(cell._globals, cell2_pos_shared[0].x);
    atomicAdd(cell._globals + 1, cell2_pos_shared[0].y);
    atomicAdd(cell._globals + 2, cell2_pos_shared[0].z);
    atomicAdd(cell._globals + 3, cell2_pos_shared[0].w);
  }
}

/**
 * Calculates interactions using linked cells with newton3
 * @tparam floatType
 * @tparam block_size cuda block size
 * @param cell particle storage
 * @param cids cell ids to calculate interactions for
 * @param cellSizes sizes of the cells by id
 */
template <typename floatType, int block_size>
__global__ void LinkedCellsTraversalN3(LJFunctorCudaGlobalsSoA<floatType> cell, unsigned int *cids, size_t *cellSizes) {
  unsigned int own_cid = cids[blockIdx.x];
  __shared__ typename vec4<floatType>::Type cell2_pos_shared[block_size];
  __shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];

  typename vec3<floatType>::Type myposition = {getInfinity<floatType>(), getInfinity<floatType>(),
                                               getInfinity<floatType>()};
  typename vec3<floatType>::Type myf = {0, 0, 0};
  typename vec4<floatType>::Type myglobals = {0, 0, 0, 0};
  floatType isOwned = 0.0;

  int index = cellSizes[own_cid] + threadIdx.x;
  if (threadIdx.x < (cellSizes[own_cid + 1] - cellSizes[own_cid])) {
    myposition.x = cell._posX[index];
    myposition.y = cell._posY[index];
    myposition.z = cell._posZ[index];
    isOwned = (cell._owned[index]) * (floatType)0.5;
  }
  // other cells
  for (auto other_index = 0; other_index < linkedCellsOffsetsSize; ++other_index) {
    const int other_id = own_cid + linkedCellsOffsets[other_index];
    const size_t cell2Start = cellSizes[other_id];
    const int sizeCell2 = cellSizes[other_id + 1] - cell2Start;

    cell2_pos_shared[threadIdx.x] = {cell._posX[cell2Start + threadIdx.x], cell._posY[cell2Start + threadIdx.x],
                                     cell._posZ[cell2Start + threadIdx.x],
                                     (cell._owned[cell2Start + threadIdx.x]) * (floatType)0.5};
    cell2_forces_shared[threadIdx.x] = {0, 0, 0};
    __syncthreads();
    for (int j = 0; j < sizeCell2; ++j) {
      const int offset = (j + threadIdx.x) % sizeCell2;
      myf = bodyBodyFN3<floatType, false>(myposition, cell2_pos_shared[offset], myf, cell2_forces_shared + offset,
                                          myglobals, isOwned);
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
      myf = bodyBodyF<floatType, true>(myposition, cell2_pos_shared[j], myf, myglobals, isOwned != 0.0);
    }
    __syncthreads();
  }

  // reduce globals
  cell2_pos_shared[threadIdx.x] = myglobals;
  __syncthreads();

  reduceGlobalsShared<floatType, block_size>(cell2_pos_shared);
  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(cell._globals, cell2_pos_shared[0].x);
    atomicAdd(cell._globals + 1, cell2_pos_shared[0].y);
    atomicAdd(cell._globals + 2, cell2_pos_shared[0].z);
    atomicAdd(cell._globals + 3, cell2_pos_shared[0].w);
  }

  if (threadIdx.x < (cellSizes[own_cid + 1] - cellSizes[own_cid])) {
    atomicAdd(cell._forceX + index, myf.x);
    atomicAdd(cell._forceY + index, myf.y);
    atomicAdd(cell._forceZ + index, myf.z);
  }
}

/**
 * Calls LinkedCellsTraversalNoN3
 * @tparam floatType
 * @param cell1Base particle storage
 * @param reqThreads max number of particles over all cells
 * @param cids_size number of cells to traverse
 * @param cids ids of the cells to traverse
 * @param cellSizes_size number of cells
 * @param cellSizes sizes for all cells
 * @param stream cuda stream to start the kernel on
 */
template <typename floatType>
void LJFunctorCudaGlobalsWrapper<floatType>::LinkedCellsTraversalNoN3Wrapper(FunctorCudaSoA<floatType> *cell1Base,
                                                                             unsigned int reqThreads,
                                                                             unsigned int cids_size, unsigned int *cids,
                                                                             unsigned int cellSizes_size,
                                                                             size_t *cellSizes, cudaStream_t stream) {
  LJFunctorCudaGlobalsSoA<floatType> cell1 = *static_cast<LJFunctorCudaGlobalsSoA<floatType> *>(cell1Base);

  switch (reqThreads) {
    CREATESWITCHCASES(cids_size, 0, LinkedCellsTraversalNoN3, (cell1, cids, cellSizes));
    default:
      autopas::utils::ExceptionHandler::exception(
          "Linked Cells NoN3: cuda Kernel size not available for Linked cells. Too many particles "
          "in a cell. Requested: {}",
          reqThreads);
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

/**
 * Calls LinkedCellsTraversalN3
 * @tparam floatType
 * @param cell1Base particle storage
 * @param reqThreads max number of particles over all cells
 * @param cids_size number of cells to traverse
 * @param cids ids of the cells to traverse
 * @param cellSizes_size number of cells
 * @param cellSizes sizes for all cells
 * @param stream cuda stream to start the kernel on
 */
template <typename floatType>
void LJFunctorCudaGlobalsWrapper<floatType>::LinkedCellsTraversalN3Wrapper(FunctorCudaSoA<floatType> *cell1Base,
                                                                           unsigned int reqThreads,
                                                                           unsigned int cids_size, unsigned int *cids,
                                                                           unsigned int cellSizes_size,
                                                                           size_t *cellSizes, cudaStream_t stream) {
  LJFunctorCudaGlobalsSoA<floatType> cell1 = *static_cast<LJFunctorCudaGlobalsSoA<floatType> *>(cell1Base);

  switch (reqThreads) {
    CREATESWITCHCASES(cids_size, 0, LinkedCellsTraversalN3, (cell1, cids, cellSizes));
    default:
      autopas::utils::ExceptionHandler::exception(
          "Linked Cells N3:cuda Kernel size not available for Linked cells. Too many particles in "
          "a cell. Requested: {}",
          reqThreads);
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

/**
 * Calculates interactions between equal sized clusters using neighbor lists
 * @tparam floatType
 * @tparam block_size cuda block size
 * @param cell particle storage
 * @param others_size max number of neighbors to address other ids
 * @param other_ids array of neighbor cells ids for a cell i (others_size * i) terminated by UINT_MAX
 */
template <typename floatType, int block_size>
__global__ void CellVerletTraversalNoN3(LJFunctorCudaGlobalsSoA<floatType> cell, const unsigned int others_size,
                                        unsigned int *other_ids) {
  __shared__ typename vec4<floatType>::Type cell2_pos_shared[block_size];
  typename vec3<floatType>::Type myf = {0, 0, 0};
  typename vec4<floatType>::Type myglobals = {0, 0, 0, 0};

  unsigned int index = blockIdx.x * block_size + threadIdx.x;
  typename vec3<floatType>::Type myposition = {cell._posX[index], cell._posY[index], cell._posZ[index]};
  bool isOwned = cell._owned[index] == (floatType)1.0;

  // other cells
  unsigned int cid;
  for (auto other_index = others_size * blockIdx.x; (cid = other_ids[other_index]) < UINT_MAX; ++other_index) {
    const size_t own_particle = block_size * cid + threadIdx.x;
    cell2_pos_shared[threadIdx.x] = {cell._posX[own_particle], cell._posY[own_particle], cell._posZ[own_particle]};
    __syncthreads();
    for (int j = 0; j < block_size; ++j) {
      myf = bodyBodyF<floatType>(myposition, cell2_pos_shared[j], myf, myglobals, isOwned);
    }
    __syncthreads();
  }

  // reduce globals
  cell2_pos_shared[threadIdx.x] = myglobals;
  __syncthreads();

  reduceGlobalsShared<floatType, block_size>(cell2_pos_shared);
  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(cell._globals, cell2_pos_shared[0].x);
    atomicAdd(cell._globals + 1, cell2_pos_shared[0].y);
    atomicAdd(cell._globals + 2, cell2_pos_shared[0].z);
    atomicAdd(cell._globals + 3, cell2_pos_shared[0].w);
  }

  atomicAdd(cell._forceX + index, myf.x);
  atomicAdd(cell._forceY + index, myf.y);
  atomicAdd(cell._forceZ + index, myf.z);
}

/**
 * Calls CellVerletTraversalNoN3
 * @tparam floatType
 * @param cell1Base particle storage
 * @param ncells number of clusters
 * @param clusterSize size of the clusters
 * @param others_size maximum number of neighbors for a cluster
 * @param other_ids neighbor ids
 * @param stream cuda stream to start the kernel on
 */
template <typename floatType>
void LJFunctorCudaGlobalsWrapper<floatType>::CellVerletTraversalNoN3Wrapper(
    FunctorCudaSoA<floatType> *cell1Base, unsigned int ncells, unsigned int clusterSize, unsigned int others_size,
    unsigned int *other_ids, cudaStream_t stream) {
  LJFunctorCudaGlobalsSoA<floatType> cell1 = *static_cast<LJFunctorCudaGlobalsSoA<floatType> *>(cell1Base);
  switch (clusterSize) {
    CREATESWITCHCASES(ncells, 0, CellVerletTraversalNoN3, (cell1, others_size, other_ids));
    default:
      autopas::utils::ExceptionHandler::exception(
          "cuda Kernel size not available for Verlet cells. Too many particles in a cell. "
          "Requested: {}",
          clusterSize);
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

/**
 * Calculates interactions between equal sized clusters using neighbor lists using newton 3
 * @tparam floatType
 * @tparam block_size cuda block size
 * @param cell particle storage
 * @param others_size max number of neighbors to address other ids
 * @param other_ids array of neighbor cells ids for a cell i (others_size * i) terminated by UINT_MAX
 */
template <typename floatType, int block_size>
__global__ void CellVerletTraversalN3(LJFunctorCudaGlobalsSoA<floatType> cell, unsigned int others_size,
                                      unsigned int *other_ids) {
  const unsigned int mask = block_size - 1;

  __shared__ typename vec4<floatType>::Type cell2_pos_shared[block_size];
  __shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];

  typename vec3<floatType>::Type myf = {0, 0, 0};
  typename vec4<floatType>::Type myglobals = {0, 0, 0, 0};
  floatType isOwned = 0.0;

  int index = blockIdx.x * block_size + threadIdx.x;
  typename vec3<floatType>::Type myposition = {cell._posX[index], cell._posY[index], cell._posZ[index]};
  isOwned = cell._owned[index] * (floatType)0.5;

  // other cells
  unsigned int cid;
  for (auto other_index = others_size * blockIdx.x; (cid = other_ids[other_index]) != UINT_MAX; ++other_index) {
    const unsigned int cell2Start = block_size * cid;

    cell2_pos_shared[threadIdx.x] = {cell._posX[cell2Start + threadIdx.x], cell._posY[cell2Start + threadIdx.x],
                                     cell._posZ[cell2Start + threadIdx.x],
                                     cell._owned[cell2Start + threadIdx.x] * (floatType)0.5};
    cell2_forces_shared[threadIdx.x] = {0, 0, 0};
    __syncthreads();
    for (int j = 0; j < block_size; ++j) {
      unsigned int offset = 0;
      if ((block_size & (block_size - 1)) == 0) {
        offset = (j + threadIdx.x) & mask;
      } else {
        offset = (j + threadIdx.x) % block_size;
      }
      myf = bodyBodyFN3<floatType, false>(myposition, cell2_pos_shared[offset], myf, cell2_forces_shared + offset,
                                          myglobals, isOwned);
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
      myf = bodyBodyF<floatType, true>(myposition, cell2_pos_shared[j], myf, myglobals, isOwned != 0.0);
    }
    __syncthreads();
  }

  // reduce globals
  cell2_pos_shared[threadIdx.x] = myglobals;
  __syncthreads();

  reduceGlobalsShared<floatType, block_size>(cell2_pos_shared);
  __syncthreads();

  if (threadIdx.x == 0) {
    atomicAdd(cell._globals, cell2_pos_shared[0].x);
    atomicAdd(cell._globals + 1, cell2_pos_shared[0].y);
    atomicAdd(cell._globals + 2, cell2_pos_shared[0].z);
    atomicAdd(cell._globals + 3, cell2_pos_shared[0].w);
  }

  atomicAdd(cell._forceX + index, myf.x);
  atomicAdd(cell._forceY + index, myf.y);
  atomicAdd(cell._forceZ + index, myf.z);
}

/**
 * Calls CellVerletTraversalN3
 * @tparam floatType
 * @param cell1Base particle storage
 * @param ncells number of clusters
 * @param clusterSize size of the clusters
 * @param others_size maximum number of neighbors for a cluster
 * @param other_ids neighbor ids
 * @param stream cuda stream to start the kernel on
 */
template <typename floatType>
void LJFunctorCudaGlobalsWrapper<floatType>::CellVerletTraversalN3Wrapper(FunctorCudaSoA<floatType> *cell1Base,
                                                                          unsigned int ncells, unsigned int clusterSize,
                                                                          unsigned int others_size,
                                                                          unsigned int *other_ids,
                                                                          cudaStream_t stream) {
  LJFunctorCudaGlobalsSoA<floatType> cell1 = *static_cast<LJFunctorCudaGlobalsSoA<floatType> *>(cell1Base);
  switch (clusterSize) {
    CREATESWITCHCASES(ncells, 0, CellVerletTraversalN3, (cell1, others_size, other_ids));
    default:
      autopas::utils::ExceptionHandler::exception(
          "cuda Kernel size not available for Verlet cells. Too many particles in a cell. "
          "Requested: {}",
          clusterSize);
      break;
  }
  autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

/**
 * float version to load constants to constant memory
 * @param constants LJ FP32 constants
 */
template <>
void LJFunctorCudaGlobalsWrapper<float>::loadConstants(FunctorCudaConstants<float> *constants) {
  LJFunctorConstants<float> *c = static_cast<LJFunctorConstants<float> *>(constants);

  cudaMemcpyToSymbol(global_constants_float, c, sizeof(LJFunctorConstants<float>));
}

/**
 * double version to load constants to constant memory
 * @param constants LJ FP64 constants
 */
template <>
void LJFunctorCudaGlobalsWrapper<double>::loadConstants(FunctorCudaConstants<double> *constants) {
  LJFunctorConstants<double> *c = static_cast<LJFunctorConstants<double> *>(constants);
  autopas::utils::CudaExceptionHandler::checkErrorCode(
      cudaMemcpyToSymbol(global_constants_double, c, sizeof(LJFunctorConstants<double>)));
}

/**
 * constant loader for other data types throws exception
 * @param constants
 */
template <typename T>
void LJFunctorCudaGlobalsWrapper<T>::loadConstants(FunctorCudaConstants<T> *constants) {
  autopas::utils::ExceptionHandler::exception("Cuda constants with unknown Type loaded");
}

/**
 * Load offsets needed by linked cells traversal to constant memory
 * @param offsets_size [1, .. 27]
 * @param offsets offsets to neighbor cells
 */
template <typename T>
void LJFunctorCudaGlobalsWrapper<T>::loadLinkedCellsOffsets(unsigned int offsets_size, int *offsets) {
  if (offsets_size > 27) {
    autopas::utils::ExceptionHandler::exception(
        "LJFunctorCudaGlobalsWrapper does not support linked cells with >27 neighbors");
  }
  autopas::utils::CudaExceptionHandler::checkErrorCode(
      cudaMemcpyToSymbol(linkedCellsOffsetsSize, &offsets_size, sizeof(unsigned int)));
  autopas::utils::CudaExceptionHandler::checkErrorCode(
      cudaMemcpyToSymbol(linkedCellsOffsets, offsets, offsets_size * sizeof(int)));
}

// Instantiation of Wrapper with float precision
template class LJFunctorCudaGlobalsWrapper<float>;
// Instantiation of Wrapper with double precision
template class LJFunctorCudaGlobalsWrapper<double>;

}  // namespace autopas
