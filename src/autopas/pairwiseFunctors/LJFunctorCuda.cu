/**
 * @file LJFunctorCuda.h
 *
 * @date 9.1.2019
 * @author jspahl
 */
#include "autopas/pairwiseFunctors/LJFunctorCuda.cuh"
#include <iostream>
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/CudaExceptionHandler.h"
#include "math_constants.h"

namespace autopas {

__constant__ constants<float> global_constants_float;
__constant__ constants<double> global_constants_double;

template<typename T>
__device__ inline constants<T>& getConstants() {
	return global_constants_float;
}
template<>
__device__ inline constants<double>& getConstants<double>() {
	return global_constants_double;
}

template<typename T>
__device__ inline T getInfinity() {
	return CUDART_INF_F;
}
template<>
__device__ inline double getInfinity<double>() {
	return CUDART_INF;
}

template<typename floatType>
__device__
inline typename vec3<floatType>::Type bodyBodyF(
		typename vec3<floatType>::Type i, typename vec3<floatType>::Type j,
		typename vec3<floatType>::Type fi) {
	floatType drx = i.x - j.x;
	floatType dry = i.y - j.y;
	floatType drz = i.z - j.z;

	floatType dr2 = drx * drx + dry * dry + drz * drz;

	if (dr2 > getConstants<floatType>().cutoffsquare | dr2 == 0.0) {
		return fi;
	}

	floatType invdr2 = 1. / dr2;
	floatType lj6 = getConstants<floatType>().sigmasquare * invdr2;
	lj6 = lj6 * lj6 * lj6;
	floatType lj12 = lj6 * lj6;
	floatType lj12m6 = lj12 - lj6;
	floatType fac = getConstants<floatType>().epsilon24 * (lj12 + lj12m6)
			* invdr2;

	fi.x += drx * fac;
	fi.y += dry * fac;
	fi.z += drz * fac;

	return fi;
}

template<typename floatType, bool n3AdditionSafe = false>
__device__
inline typename vec3<floatType>::Type bodyBodyFN3(
		typename vec3<floatType>::Type i, typename vec3<floatType>::Type j,
		typename vec3<floatType>::Type fi, typename vec3<floatType>::Type* fj) {
	floatType drx = i.x - j.x;
	floatType dry = i.y - j.y;
	floatType drz = i.z - j.z;

	floatType dr2 = drx * drx + dry * dry + drz * drz;

	if (dr2 > getConstants<floatType>().cutoffsquare) {
		return fi;
	}

	floatType invdr2 = 1. / dr2;
	floatType lj6 = getConstants<floatType>().sigmasquare * invdr2;
	lj6 = lj6 * lj6 * lj6;
	floatType lj12 = lj6 * lj6;
	floatType lj12m6 = lj12 - lj6;
	floatType fac = getConstants<floatType>().epsilon24 * (lj12 + lj12m6)
			* invdr2;

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

template<typename floatType, int block_size>
__global__
void SoAFunctorNoN3(LJFunctorCudaSoA<floatType> cell1) {
	__shared__ typename vec3<floatType>::Type block_pos[block_size];
	int i, tile;
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	typename vec3<floatType>::Type myposition = { 0, 0, 0 };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };
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

template<typename floatType, int block_size>
__global__
void SoAFunctorNoN3Pair(LJFunctorCudaSoA<floatType> cell1,
		LJFunctorCudaSoA<floatType> cell2) {
	__shared__ typename vec3<floatType>::Type block_pos[block_size];
	int i, tile;
	int tid = blockIdx.x * block_size + threadIdx.x;
	typename vec3<floatType>::Type myposition;
	typename vec3<floatType>::Type myf = { 0, 0, 0 };

	if (tid < cell1._size) {
		myposition.x = cell1._posX[tid];
		myposition.y = cell1._posY[tid];
		myposition.z = cell1._posZ[tid];
	}

	for (i = 0, tile = 0; i < cell2._size; i += block_size, ++tile) {
		int idx = tile * block_size + threadIdx.x;

		if (idx < cell2._size)
			block_pos[threadIdx.x] = {cell2._posX[idx], cell2._posY[idx], cell2._posZ[idx]};
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

template<typename floatType, int block_size, bool NMisMultipleBlockSize = false>
__global__
void SoAFunctorN3(LJFunctorCudaSoA<floatType> cell1) {
	static_assert((block_size & (block_size - 1)) == 0, "block size must be power of 2");
	__shared__ typename vec3<floatType>::Type cell1_pos_shared[block_size];
	__shared__ typename vec3<floatType>::Type cell1_forces_shared[block_size];
	int tid = blockIdx.x * block_size + threadIdx.x;
	typename vec3<floatType>::Type myposition = { getInfinity<floatType>(),
			getInfinity<floatType>(), getInfinity<floatType>() };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };
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
		cell1_forces_shared[threadIdx.x] = {0,0,0};
		__syncthreads();

		for (int j = 0; j < block_size; ++j) {
			const int offset = (j + threadIdx.x) & mask;
			myf = bodyBodyFN3<floatType>(myposition, cell1_pos_shared[offset],
					myf, cell1_forces_shared + offset);
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
		cell1_forces_shared[threadIdx.x] = {0,0,0};
		__syncthreads();

		for (int j = threadIdx.x - 1; j >= 0; --j) {
			myf = bodyBodyFN3<floatType>(myposition, cell1_pos_shared[j], myf,
					cell1_forces_shared + j);
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

template<typename floatType, int block_size, bool NMisMultipleBlockSize = false>
__global__
void SoAFunctorN3Pair(LJFunctorCudaSoA<floatType> cell1,
		LJFunctorCudaSoA<floatType> cell2) {
	static_assert((block_size & (block_size - 1)) == 0, "block size must be power of 2");
	__shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
	__shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];
	int tid = blockIdx.x * block_size + threadIdx.x;
	typename vec3<floatType>::Type myposition = { getInfinity<floatType>(),
			getInfinity<floatType>(), getInfinity<floatType>() };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };
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
		cell2_forces_shared[threadIdx.x] = {0,0,0};
		__syncthreads();

		for (int j = 0; j < block_size; ++j) {
			const int offset = (j + threadIdx.x) & mask;
			myf = bodyBodyFN3<floatType, true>(myposition,
					cell2_pos_shared[offset], myf,
					cell2_forces_shared + offset);
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
			cell2_forces_shared[threadIdx.x] = {0,0,0};
		}
		__syncthreads();

		const int size = block_size + cell2._size - i;
		for (int j = 0; j < size; ++j) {
			const int offset = (j + threadIdx.x) % size;
			myf = bodyBodyFN3<floatType>(myposition, cell2_pos_shared[offset],
					myf, cell2_forces_shared + offset);
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

template<typename floatType>
void CudaWrapper::SoAFunctorNoN3Wrapper(LJFunctorCudaSoA<floatType> cell1,
		cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		SoAFunctorNoN3<floatType, 32> <<<numRequiredBlocks(cell1._size), 32, 0,
				stream>>>(cell1);
		break;
	case 64:
		SoAFunctorNoN3<floatType, 64> <<<numRequiredBlocks(cell1._size), 64, 0,
				stream>>>(cell1);
		break;
	case 96:
		SoAFunctorNoN3<floatType, 96> <<<numRequiredBlocks(cell1._size), 96, 0,
				stream>>>(cell1);
		break;
	case 128:
		SoAFunctorNoN3<floatType, 128> <<<numRequiredBlocks(cell1._size), 128,
				0, stream>>>(cell1);
		break;
	case 256:
		SoAFunctorNoN3<floatType, 256> <<<numRequiredBlocks(cell1._size), 256,
				0, stream>>>(cell1);
		break;
	case 512:
		SoAFunctorNoN3<floatType, 512> <<<numRequiredBlocks(cell1._size), 512,
				0, stream>>>(cell1);
		break;
	case 1024:
		SoAFunctorNoN3<floatType, 1024> <<<numRequiredBlocks(cell1._size), 1024,
				0, stream>>>(cell1);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::SoAFunctorNoN3Wrapper<float>(
		LJFunctorCudaSoA<float> cell1, cudaStream_t stream);
template void CudaWrapper::SoAFunctorNoN3Wrapper<double>(
		LJFunctorCudaSoA<double> cell1, cudaStream_t stream);

template<typename floatType>
void CudaWrapper::SoAFunctorNoN3PairWrapper(LJFunctorCudaSoA<floatType> cell1,
		LJFunctorCudaSoA<floatType> cell2, cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		SoAFunctorNoN3Pair<floatType, 32> <<<numRequiredBlocks(cell1._size), 32,
				0, stream>>>(cell1, cell2);
		break;
	case 64:
		SoAFunctorNoN3Pair<floatType, 64> <<<numRequiredBlocks(cell1._size), 64,
				0, stream>>>(cell1, cell2);
		break;
	case 96:
		SoAFunctorNoN3Pair<floatType, 96> <<<numRequiredBlocks(cell1._size), 96,
				0, stream>>>(cell1, cell2);
		break;
	case 128:
		SoAFunctorNoN3Pair<floatType, 128> <<<numRequiredBlocks(cell1._size),
				128, 0, stream>>>(cell1, cell2);
		break;
	case 256:
		SoAFunctorNoN3Pair<floatType, 256> <<<numRequiredBlocks(cell1._size),
				256, 0, stream>>>(cell1, cell2);
		break;
	case 512:
		SoAFunctorNoN3Pair<floatType, 512> <<<numRequiredBlocks(cell1._size),
				512, 0, stream>>>(cell1, cell2);
		break;
	case 1024:
		SoAFunctorNoN3Pair<floatType, 1024> <<<numRequiredBlocks(cell1._size),
				1024, 0, stream>>>(cell1, cell2);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::SoAFunctorNoN3PairWrapper<float>(
		LJFunctorCudaSoA<float> cell1, LJFunctorCudaSoA<float> cell2,
		cudaStream_t stream);
template void CudaWrapper::SoAFunctorNoN3PairWrapper<double>(
		LJFunctorCudaSoA<double> cell1, LJFunctorCudaSoA<double> cell2,
		cudaStream_t stream);

template<typename floatType>
void CudaWrapper::SoAFunctorN3Wrapper(LJFunctorCudaSoA<floatType> cell1,
		cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		SoAFunctorN3<floatType, 32> <<<numRequiredBlocks(cell1._size), 32, 0,
				stream>>>(cell1);
		break;
	case 64:
		SoAFunctorN3<floatType, 64> <<<numRequiredBlocks(cell1._size), 64, 0,
				stream>>>(cell1);
		break;
	case 128:
		SoAFunctorN3<floatType, 128> <<<numRequiredBlocks(cell1._size), 128, 0,
				stream>>>(cell1);
		break;
	case 256:
		SoAFunctorN3<floatType, 256> <<<numRequiredBlocks(cell1._size), 256, 0,
				stream>>>(cell1);
		break;
	case 512:
		SoAFunctorN3<floatType, 512> <<<numRequiredBlocks(cell1._size), 512, 0,
				stream>>>(cell1);
		break;
	case 1024:
		SoAFunctorN3<floatType, 1024> <<<numRequiredBlocks(cell1._size), 1024,
				0, stream>>>(cell1);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::SoAFunctorN3Wrapper<float>(
		LJFunctorCudaSoA<float> cell1, cudaStream_t stream);
template void CudaWrapper::SoAFunctorN3Wrapper<double>(
		LJFunctorCudaSoA<double> cell1, cudaStream_t stream);

template<typename floatType>
void CudaWrapper::SoAFunctorN3PairWrapper(LJFunctorCudaSoA<floatType> cell1,
		LJFunctorCudaSoA<floatType> cell2, cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		SoAFunctorN3Pair<floatType, 32> <<<numRequiredBlocks(cell1._size), 32,
				0, stream>>>(cell1, cell2);
		break;
	case 64:
		SoAFunctorN3Pair<floatType, 64> <<<numRequiredBlocks(cell1._size), 64,
				0, stream>>>(cell1, cell2);
		break;
	case 128:
		SoAFunctorN3Pair<floatType, 128> <<<numRequiredBlocks(cell1._size), 128,
				0, stream>>>(cell1, cell2);
		break;
	case 256:
		SoAFunctorN3Pair<floatType, 256> <<<numRequiredBlocks(cell1._size), 256,
				0, stream>>>(cell1, cell2);
		break;
	case 512:
		SoAFunctorN3Pair<floatType, 512> <<<numRequiredBlocks(cell1._size), 512,
				0, stream>>>(cell1, cell2);
		break;
	case 1024:
		SoAFunctorN3Pair<floatType, 1024> <<<numRequiredBlocks(cell1._size),
				1024, 0, stream>>>(cell1, cell2);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();

}

template void CudaWrapper::SoAFunctorN3PairWrapper<float>(
		LJFunctorCudaSoA<float> cell1, LJFunctorCudaSoA<float> cell2,
		cudaStream_t stream);
template void CudaWrapper::SoAFunctorN3PairWrapper<double>(
		LJFunctorCudaSoA<double> cell1, LJFunctorCudaSoA<double> cell2,
		cudaStream_t stream);

template<typename floatType, int block_size>
__global__
void LinkedCellsTraversalNoN3(LJFunctorCudaSoA<floatType> cell,
		unsigned int* cids, size_t* cellSizes, unsigned int offsets_size,
		int* offsets) {

	int own_cid = cids[blockIdx.x];
	__shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
	typename vec3<floatType>::Type myposition = { getInfinity<floatType>(),
			getInfinity<floatType>(), getInfinity<floatType>() };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };

	int index = cellSizes[own_cid] + threadIdx.x;
	if (threadIdx.x < (cellSizes[own_cid + 1] - cellSizes[own_cid])) {
		myposition.x = cell._posX[index];
		myposition.y = cell._posY[index];
		myposition.z = cell._posZ[index];
	}
	//other cells
	for (auto other_index = 0; other_index < offsets_size; ++other_index) {
		const int other_id = own_cid + offsets[other_index];
		const int cell2Start = cellSizes[other_id];
		const int sizeCell2 = cellSizes[other_id + 1] - cell2Start;

		cell2_pos_shared[threadIdx.x] = {cell._posX[cell2Start+threadIdx.x], cell._posY[cell2Start+threadIdx.x], cell._posZ[cell2Start+threadIdx.x]};
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

template<typename floatType, int block_size>
__global__
void LinkedCellsTraversalN3(LJFunctorCudaSoA<floatType> cell,
		unsigned int* cids, size_t* cellSizes, unsigned int offsets_size,
		int* offsets) {

	int own_cid = cids[blockIdx.x];
	__shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
	__shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];

	typename vec3<floatType>::Type myposition = { getInfinity<floatType>(),
			getInfinity<floatType>(), getInfinity<floatType>() };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };

	int index = cellSizes[own_cid] + threadIdx.x;
	if (threadIdx.x < (cellSizes[own_cid + 1] - cellSizes[own_cid])) {
		myposition.x = cell._posX[index];
		myposition.y = cell._posY[index];
		myposition.z = cell._posZ[index];
	}
	//other cells
	for (auto other_index = 0; other_index < offsets_size; ++other_index) {
		const int other_id = own_cid + offsets[other_index];
		const int cell2Start = cellSizes[other_id];
		const int sizeCell2 = cellSizes[other_id + 1] - cell2Start;

		cell2_pos_shared[threadIdx.x] = {cell._posX[cell2Start+threadIdx.x], cell._posY[cell2Start+threadIdx.x], cell._posZ[cell2Start+threadIdx.x]};
		cell2_forces_shared[threadIdx.x] = {0,0,0};
		__syncthreads();
		for (int j = 0; j < sizeCell2; ++j) {
			const int offset = (j + threadIdx.x) % sizeCell2;
			myf = bodyBodyFN3<floatType, false>(myposition,
					cell2_pos_shared[offset], myf,
					cell2_forces_shared + offset);
		}
		__syncthreads();

		atomicAdd(cell._forceX + cell2Start + threadIdx.x,
				cell2_forces_shared[threadIdx.x].x);
		atomicAdd(cell._forceY + cell2Start + threadIdx.x,
				cell2_forces_shared[threadIdx.x].y);
		atomicAdd(cell._forceZ + cell2Start + threadIdx.x,
				cell2_forces_shared[threadIdx.x].z);
		__syncthreads();
	}
	//same cells
	{
		const int cell1Start = cellSizes[own_cid];
		const int sizeCell1 = cellSizes[own_cid + 1] - cell1Start;

		cell2_pos_shared[threadIdx.x] = {cell._posX[cell1Start+threadIdx.x], cell._posY[cell1Start+threadIdx.x], cell._posZ[cell1Start+threadIdx.x]};
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

template<typename floatType>
void CudaWrapper::LinkedCellsTraversalNoN3Wrapper(
		LJFunctorCudaSoA<floatType> cell1, unsigned int cids_size,
		unsigned int* cids, unsigned int cellSizes_size, size_t* cellSizes,
		unsigned int offsets_size, int* offsets, cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		LinkedCellsTraversalNoN3<floatType, 32> <<<cids_size, 32, 0, stream>>>(
				cell1, cids, cellSizes, offsets_size, offsets);
		break;
	case 64:
		LinkedCellsTraversalNoN3<floatType, 64> <<<cids_size, 64, 0, stream>>>(
				cell1, cids, cellSizes, offsets_size, offsets);
		break;
	case 96:
		LinkedCellsTraversalNoN3<floatType, 96> <<<cids_size, 96, 0, stream>>>(
				cell1, cids, cellSizes, offsets_size, offsets);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				"cuda Kernel size not available for Linked cells available 32, 64, 96. Too many particles in a cell. Requested: {}",
				_num_threads);
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::LinkedCellsTraversalNoN3Wrapper<float>(
		LJFunctorCudaSoA<float> cell1, unsigned int cids_size,
		unsigned int* cids, unsigned int cellSizes_size, size_t* cellSizes,
		unsigned int offsets_size, int* offsets, cudaStream_t stream);

template void CudaWrapper::LinkedCellsTraversalNoN3Wrapper<double>(
		LJFunctorCudaSoA<double> cell1, unsigned int cids_size,
		unsigned int* cids, unsigned int cellSizes_size, size_t* cellSizes,
		unsigned int offsets_size, int* offsets, cudaStream_t stream);

template<typename floatType>
void CudaWrapper::LinkedCellsTraversalN3Wrapper(
		LJFunctorCudaSoA<floatType> cell1, unsigned int cids_size,
		unsigned int* cids, unsigned int cellSizes_size, size_t* cellSizes,
		unsigned int offsets_size, int* offsets, cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		LinkedCellsTraversalN3<floatType, 32> <<<cids_size, 32, 0, stream>>>(
				cell1, cids, cellSizes, offsets_size, offsets);
		break;
	case 64:
		LinkedCellsTraversalN3<floatType, 64> <<<cids_size, 64, 0, stream>>>(
				cell1, cids, cellSizes, offsets_size, offsets);
		break;
	case 96:
		LinkedCellsTraversalN3<floatType, 96> <<<cids_size, 96, 0, stream>>>(
				cell1, cids, cellSizes, offsets_size, offsets);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				"cuda Kernel size not available for Linked cells available 32, 64, 96. Too many particles in a cell. Requested: {}",
				_num_threads);
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::LinkedCellsTraversalN3Wrapper<float>(
		LJFunctorCudaSoA<float> cell1, unsigned int cids_size,
		unsigned int* cids, unsigned int cellSizes_size, size_t* cellSizes,
		unsigned int offsets_size, int* offsets, cudaStream_t stream);

template void CudaWrapper::LinkedCellsTraversalN3Wrapper<double>(
		LJFunctorCudaSoA<double> cell1, unsigned int cids_size,
		unsigned int* cids, unsigned int cellSizes_size, size_t* cellSizes,
		unsigned int offsets_size, int* offsets, cudaStream_t stream);

template<typename floatType, int block_size>
__global__
void CellVerletTraversalNoN3(LJFunctorCudaSoA<floatType> cell,
		unsigned int others_size, unsigned int* other_ids) {

	__shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
	typename vec3<floatType>::Type myposition = { getInfinity<floatType>(),
			getInfinity<floatType>(), getInfinity<floatType>() };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };

	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
	myposition.x = cell._posX[index];
	myposition.y = cell._posY[index];
	myposition.z = cell._posZ[index];

	//other cells
	for (auto other_index = others_size * blockIdx.x;
			other_ids[other_index] < UINT_MAX; ++other_index) {
		const unsigned int cell2Start = blockDim.x * other_ids[other_index];

		cell2_pos_shared[threadIdx.x] = {cell._posX[cell2Start+threadIdx.x], cell._posY[cell2Start+threadIdx.x], cell._posZ[cell2Start+threadIdx.x]};
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

template<typename floatType>
void CudaWrapper::CellVerletTraversalNoN3Wrapper(
		LJFunctorCudaSoA<floatType> cell1, unsigned int ncells,
		unsigned int clusterSize, unsigned int others_size,
		unsigned int* other_ids, cudaStream_t stream) {
	switch (clusterSize) {
	case 32:
		CellVerletTraversalNoN3<floatType, 32> <<<ncells, 32, 0, stream>>>(
				cell1, others_size, other_ids);
		break;
	case 64:
		CellVerletTraversalNoN3<floatType, 64> <<<ncells, 64, 0, stream>>>(
				cell1, others_size, other_ids);
		break;
	case 96:
		CellVerletTraversalNoN3<floatType, 96> <<<ncells, 96, 0, stream>>>(
				cell1, others_size, other_ids);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				"cuda Kernel size not available for Verlet cells available 32, 64, 96. Too many particles in a cell. Requested: {}",
				_num_threads);
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::CellVerletTraversalNoN3Wrapper<float>(
		LJFunctorCudaSoA<float> cell1, unsigned int ncells,
		unsigned int clusterSize, unsigned int others_size,
		unsigned int* other_ids, cudaStream_t stream);

template void CudaWrapper::CellVerletTraversalNoN3Wrapper<double>(
		LJFunctorCudaSoA<double> cell1, unsigned int ncells,
		unsigned int clusterSize, unsigned int others_size,
		unsigned int* other_ids, cudaStream_t stream);

template<typename floatType, int block_size>
__global__
void CellVerletTraversalN3(LJFunctorCudaSoA<floatType> cell,
		unsigned int others_size, unsigned int* other_ids) {
	const unsigned int mask = block_size - 1;

	__shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
	__shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];

	typename vec3<floatType>::Type myposition = { getInfinity<floatType>(),
			getInfinity<floatType>(), getInfinity<floatType>() };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };

	int index = blockIdx.x * blockDim.x + threadIdx.x;
	myposition.x = cell._posX[index];
	myposition.y = cell._posY[index];
	myposition.z = cell._posZ[index];

	//other cells
	for (auto other_index = others_size * blockIdx.x;
			other_ids[other_index] != UINT_MAX; ++other_index) {
		const unsigned int cell2Start = blockDim.x * other_ids[other_index];

		cell2_pos_shared[threadIdx.x] = {cell._posX[cell2Start+threadIdx.x], cell._posY[cell2Start+threadIdx.x], cell._posZ[cell2Start+threadIdx.x]};
		cell2_forces_shared[threadIdx.x] = {0,0,0};
		__syncthreads();
		for (int j = 0; j < block_size; ++j) {
			unsigned int offset = 0;
			if ((block_size & (block_size - 1)) == 0) {
				offset = (j + threadIdx.x) & mask;
			} else {
				offset = (j + threadIdx.x) % block_size;
			}
			myf = bodyBodyFN3<floatType, false>(myposition,
					cell2_pos_shared[offset], myf,
					cell2_forces_shared + offset);
		}
		__syncthreads();

		atomicAdd(cell._forceX + cell2Start + threadIdx.x,
				cell2_forces_shared[threadIdx.x].x);
		atomicAdd(cell._forceY + cell2Start + threadIdx.x,
				cell2_forces_shared[threadIdx.x].y);
		atomicAdd(cell._forceZ + cell2Start + threadIdx.x,
				cell2_forces_shared[threadIdx.x].z);
		__syncthreads();
	}

	//same cluster without N3
	{
		const unsigned int cellStart = blockIdx.x * blockDim.x;

		cell2_pos_shared[threadIdx.x] = {cell._posX[cellStart+threadIdx.x], cell._posY[cellStart+threadIdx.x], cell._posZ[cellStart+threadIdx.x]};
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

template<typename floatType>
void CudaWrapper::CellVerletTraversalN3Wrapper(
		LJFunctorCudaSoA<floatType> cell1, unsigned int ncells,
		unsigned int clusterSize, unsigned int others_size,
		unsigned int* other_ids, cudaStream_t stream) {
	switch (clusterSize) {
	case 32:
		CellVerletTraversalN3<floatType, 32> <<<ncells, 32, 0, stream>>>(cell1,
				others_size, other_ids);
		break;
	case 64:
		CellVerletTraversalN3<floatType, 64> <<<ncells, 64, 0, stream>>>(cell1,
				others_size, other_ids);
		break;
	case 96:
		CellVerletTraversalN3<floatType, 96> <<<ncells, 96, 0, stream>>>(cell1,
				others_size, other_ids);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				"cuda Kernel size not available for Verlet cells available 32, 64, 96. Too many particles in a cell. Requested: {}",
				_num_threads);
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::CellVerletTraversalN3Wrapper<float>(
		LJFunctorCudaSoA<float> cell1, unsigned int ncells,
		unsigned int clusterSize, unsigned int others_size,
		unsigned int* other_ids, cudaStream_t stream);

template void CudaWrapper::CellVerletTraversalN3Wrapper<double>(
		LJFunctorCudaSoA<double> cell1, unsigned int ncells,
		unsigned int clusterSize, unsigned int others_size,
		unsigned int* other_ids, cudaStream_t stream);

template<typename floatType>
void CudaWrapper::loadConstants(floatType cutoffsquare, floatType epsilon24,
		floatType sigmasquare) {

	constants<floatType> c;
	c.cutoffsquare = cutoffsquare;
	c.epsilon24 = epsilon24;
	c.sigmasquare = sigmasquare;
	cudaMemcpyToSymbol(global_constants_double, &c,
			sizeof(constants<floatType> ));
}
template void CudaWrapper::loadConstants<float>(float cutoffsquare,
		float epsilon24, float sigmasquare);
template void CudaWrapper::loadConstants<double>(double cutoffsquare,
		double epsilon24, double sigmasquare);

} // namespace autopas
