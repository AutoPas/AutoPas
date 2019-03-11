/**
 * @file LJFunctorCuda.h
 *
 * @date 9.1.2019
 * @author jspahl
 */
#include "autopas/pairwiseFunctors/LJFunctorCuda.cuh"
#include <iostream>
#include "autopas/utils/CudaExceptionHandler.h"
#include "math_constants.h"
#include "autopas/utils/CudaDeviceVector.h"
#include "autopas/utils/CudaSoA.h"
#include "autopas/particles/Particle.h"

namespace autopas {

__constant__ constants<float> global_constants_float;
__constant__ constants<double> global_constants_double;

template<typename floatType>
class soa {
public:
	floatType* posX;
	floatType* posY;
	floatType* posZ;
	floatType* forceX;
	floatType* forceY;
	floatType* forceZ;
};

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
void SoAFunctorNoN3(int N, floatType* posX, floatType* posY, floatType* posZ,
		floatType* forceX, floatType* forceY, floatType* forceZ) {
	__shared__ typename vec3<floatType>::Type block_pos[block_size];
	int i, tile;
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	typename vec3<floatType>::Type myposition = { 0, 0, 0 };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };
	if (tid < N) {
		myposition.x = posX[tid];
		myposition.y = posY[tid];
		myposition.z = posZ[tid];
	}

	for (i = block_size, tile = 0; i < N; i += block_size, ++tile) {
		int idx = tile * block_size + threadIdx.x;

		block_pos[threadIdx.x] = {posX[idx], posY[idx], posZ[idx]};
		__syncthreads();
		if (tid < N) {
			for (int j = 0; j < blockDim.x; ++j) {
				myf = bodyBodyF<floatType>(myposition, block_pos[j], myf);
			}
		}
		__syncthreads();
	}
	{
		int idx = tile * block_size + threadIdx.x;
		block_pos[threadIdx.x] = {posX[idx], posY[idx], posZ[idx]};
		__syncthreads();

		const int size = N - tile * blockDim.x;
		for (int j = 0; j < size; ++j) {
			myf = bodyBodyF<floatType>(myposition, block_pos[j], myf);
		}

		__syncthreads();
	}

	atomicAdd(forceX + tid, myf.x);
	atomicAdd(forceY + tid, myf.y);
	atomicAdd(forceZ + tid, myf.z);
}

template<typename floatType, int block_size>
__global__
void SoAFunctorNoN3Pair(int N, floatType* posX, floatType* posY,
		floatType* posZ, floatType* forceX, floatType* forceY,
		floatType* forceZ, int M, floatType* posX2, floatType* posY2,
		floatType* posZ2) {
	__shared__ typename vec3<floatType>::Type block_pos[block_size];
	int i, tile;
	int tid = blockIdx.x * block_size + threadIdx.x;
	typename vec3<floatType>::Type myposition;
	typename vec3<floatType>::Type myf = { 0, 0, 0 };

	if (tid < N) {
		myposition.x = posX[tid];
		myposition.y = posY[tid];
		myposition.z = posZ[tid];
	}

	for (i = 0, tile = 0; i < M; i += block_size, ++tile) {
		int idx = tile * block_size + threadIdx.x;

		if (idx < M)
			block_pos[threadIdx.x] = {posX2[idx], posY2[idx], posZ2[idx]};
		__syncthreads();

		const int size = min(block_size, M - i);
		for (int j = 0; j < size; ++j) {
			myf = bodyBodyF<floatType>(myposition, block_pos[j], myf);
		}
		__syncthreads();
	}
	atomicAdd(forceX + tid, myf.x);
	atomicAdd(forceY + tid, myf.y);
	atomicAdd(forceZ + tid, myf.z);
}

template<typename floatType, int block_size, bool NMisMultipleBlockSize = false>
__global__
void SoAFunctorN3(int N, floatType* posX, floatType* posY, floatType* posZ,
		floatType* forceX, floatType* forceY, floatType* forceZ) {
	static_assert((block_size & (block_size - 1)) == 0, "block size must be power of 2");
	__shared__ typename vec3<floatType>::Type cell1_pos_shared[block_size];
	__shared__ typename vec3<floatType>::Type cell1_forces_shared[block_size];
	int tid = blockIdx.x * block_size + threadIdx.x;
	typename vec3<floatType>::Type myposition = { getInfinity<floatType>(),
			getInfinity<floatType>(), getInfinity<floatType>() };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };
	int i, tile;
	const int mask = block_size - 1;

	if (not NMisMultipleBlockSize && tid < N) {
		myposition.x = posX[tid];
		myposition.y = posY[tid];
		myposition.z = posZ[tid];
	}

	for (i = 0, tile = 0; tile < blockIdx.x; i += block_size, ++tile) {
		int idx = tile * block_size + threadIdx.x;
		cell1_pos_shared[threadIdx.x] = {posX[idx], posY[idx], posZ[idx]};
		cell1_forces_shared[threadIdx.x] = {0,0,0};
		__syncthreads();

		for (int j = 0; j < block_size; ++j) {
			const int offset = (j + threadIdx.x) & mask;
			myf = bodyBodyFN3<floatType>(myposition, cell1_pos_shared[offset],
					myf, cell1_forces_shared + offset);
		}
		__syncthreads();

		atomicAdd(forceX + idx, cell1_forces_shared[threadIdx.x].x);
		atomicAdd(forceY + idx, cell1_forces_shared[threadIdx.x].y);
		atomicAdd(forceZ + idx, cell1_forces_shared[threadIdx.x].z);
		__syncthreads();
	}

	{
		int idx = blockIdx.x * block_size + threadIdx.x;
		cell1_pos_shared[threadIdx.x] = {posX[idx], posY[idx], posZ[idx]};
		cell1_forces_shared[threadIdx.x] = {0,0,0};
		__syncthreads();

		for (int j = threadIdx.x - 1; j >= 0; --j) {
			myf = bodyBodyFN3<floatType>(myposition, cell1_pos_shared[j], myf,
					cell1_forces_shared + j);
		}
		__syncthreads();

		atomicAdd(forceX + idx, cell1_forces_shared[threadIdx.x].x);
		atomicAdd(forceY + idx, cell1_forces_shared[threadIdx.x].y);
		atomicAdd(forceZ + idx, cell1_forces_shared[threadIdx.x].z);
		__syncthreads();
	}

	atomicAdd(forceX + tid, myf.x);
	atomicAdd(forceY + tid, myf.y);
	atomicAdd(forceZ + tid, myf.z);
}

template<typename floatType, int block_size, bool NMisMultipleBlockSize = false>
__global__
void SoAFunctorN3Pair(int N, soa<floatType> cell1, int M,
		soa<floatType> cell2) {
	static_assert((block_size & (block_size - 1)) == 0, "block size must be power of 2");
	__shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
	__shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];
	int tid = blockIdx.x * block_size + threadIdx.x;
	typename vec3<floatType>::Type myposition = { getInfinity<floatType>(),
			getInfinity<floatType>(), getInfinity<floatType>() };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };
	int i, tile;
	const int mask = block_size - 1;

	if (not NMisMultipleBlockSize && tid < N) {
		myposition.x = cell1.posX[tid];
		myposition.y = cell1.posY[tid];
		myposition.z = cell1.posZ[tid];
	}
	for (i = block_size, tile = 0; i <= M; i += block_size, ++tile) {
		int idx = tile * block_size + threadIdx.x;
		cell2_pos_shared[threadIdx.x] = {cell2.posX[idx], cell2.posY[idx], cell2.posZ[idx]};
		cell2_forces_shared[threadIdx.x] = {0,0,0};
		__syncthreads();

		for (int j = 0; j < block_size; ++j) {
			const int offset = (j + threadIdx.x) & mask;
			myf = bodyBodyFN3<floatType, true>(myposition,
					cell2_pos_shared[offset], myf,
					cell2_forces_shared + offset);
		}
		__syncthreads();

		atomicAdd(cell2.forceX + idx, cell2_forces_shared[threadIdx.x].x);
		atomicAdd(cell2.forceY + idx, cell2_forces_shared[threadIdx.x].y);
		atomicAdd(cell2.forceZ + idx, cell2_forces_shared[threadIdx.x].z);
		__syncthreads();
	}
	if ((not NMisMultipleBlockSize) && (i > M)) {
		int idx = tile * block_size + threadIdx.x;
		if (idx < M) {
			cell2_pos_shared[threadIdx.x] = {cell2.posX[idx], cell2.posY[idx], cell2.posZ[idx]};
			cell2_forces_shared[threadIdx.x] = {0,0,0};
		}
		__syncthreads();

		const int size = block_size + M - i;
		for (int j = 0; j < size; ++j) {
			const int offset = (j + threadIdx.x) % size;
			myf = bodyBodyFN3<floatType>(myposition, cell2_pos_shared[offset],
					myf, cell2_forces_shared + offset);
		}
		__syncthreads();
		if (idx < M) {
			atomicAdd(cell2.forceX + idx, cell2_forces_shared[threadIdx.x].x);
			atomicAdd(cell2.forceY + idx, cell2_forces_shared[threadIdx.x].y);
			atomicAdd(cell2.forceZ + idx, cell2_forces_shared[threadIdx.x].z);
			__syncthreads();
		}
	}

	atomicAdd(cell1.forceX + tid, myf.x);
	atomicAdd(cell1.forceY + tid, myf.y);
	atomicAdd(cell1.forceZ + tid, myf.z);
}

template<typename floatType>
void CudaWrapper::SoAFunctorNoN3Wrapper(int N, floatType* posX, floatType* posY,
		floatType* posZ, floatType* forceX, floatType* forceY,
		floatType* forceZ, cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		SoAFunctorNoN3<floatType, 32> <<<numRequiredBlocks(N), 32, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 64:
		SoAFunctorNoN3<floatType, 64> <<<numRequiredBlocks(N), 64, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 96:
		SoAFunctorNoN3<floatType, 96> <<<numRequiredBlocks(N), 96, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 128:
		SoAFunctorNoN3<floatType, 128> <<<numRequiredBlocks(N), 128, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 256:
		SoAFunctorNoN3<floatType, 256> <<<numRequiredBlocks(N), 256, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 512:
		SoAFunctorNoN3<floatType, 512> <<<numRequiredBlocks(N), 512, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 1024:
		SoAFunctorNoN3<floatType, 1024> <<<numRequiredBlocks(N), 1024, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::SoAFunctorNoN3Wrapper<float>(int N, float* posX,
		float* posY, float* posZ, float* forceX, float* forceY, float* forceZ,
		cudaStream_t stream);
template void CudaWrapper::SoAFunctorNoN3Wrapper<double>(int N, double* posX,
		double* posY, double* posZ, double* forceX, double* forceY,
		double* forceZ, cudaStream_t stream);

template<typename floatType>
void CudaWrapper::SoAFunctorNoN3PairWrapper(int N, floatType* posX,
		floatType* posY, floatType* posZ, floatType* forceX, floatType* forceY,
		floatType* forceZ, int M, floatType* posX2, floatType* posY2,
		floatType* posZ2, cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		SoAFunctorNoN3Pair<floatType, 32> <<<numRequiredBlocks(N), 32, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ, M, posX2, posY2,
				posZ2);
		break;
	case 64:
		SoAFunctorNoN3Pair<floatType, 64> <<<numRequiredBlocks(N), 64, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ, M, posX2, posY2,
				posZ2);
		break;
	case 96:
		SoAFunctorNoN3Pair<floatType, 96> <<<numRequiredBlocks(N), 96, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ, M, posX2, posY2,
				posZ2);
		break;
	case 128:
		SoAFunctorNoN3Pair<floatType, 128> <<<numRequiredBlocks(N), 128, 0,
				stream>>>(N, posX, posY, posZ, forceX, forceY, forceZ, M, posX2,
				posY2, posZ2);
		break;
	case 256:
		SoAFunctorNoN3Pair<floatType, 256> <<<numRequiredBlocks(N), 256, 0,
				stream>>>(N, posX, posY, posZ, forceX, forceY, forceZ, M, posX2,
				posY2, posZ2);
		break;
	case 512:
		SoAFunctorNoN3Pair<floatType, 512> <<<numRequiredBlocks(N), 512, 0,
				stream>>>(N, posX, posY, posZ, forceX, forceY, forceZ, M, posX2,
				posY2, posZ2);
		break;
	case 1024:
		SoAFunctorNoN3Pair<floatType, 1024> <<<numRequiredBlocks(N), 1024, 0,
				stream>>>(N, posX, posY, posZ, forceX, forceY, forceZ, M, posX2,
				posY2, posZ2);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::SoAFunctorNoN3PairWrapper<float>(int N, float* posX,
		float* posY, float* posZ, float* forceX, float* forceY, float* forceZ,
		int M, float* posX2, float* posY2, float* posZ2, cudaStream_t stream);
template void CudaWrapper::SoAFunctorNoN3PairWrapper<double>(int N,
		double* posX, double* posY, double* posZ, double* forceX,
		double* forceY, double* forceZ, int M, double* posX2, double* posY2,
		double* posZ2, cudaStream_t stream);

template<typename floatType>
void CudaWrapper::SoAFunctorN3Wrapper(int N, floatType* posX, floatType* posY,
		floatType* posZ, floatType* forceX, floatType* forceY,
		floatType* forceZ, cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		SoAFunctorN3<floatType, 32> <<<numRequiredBlocks(N), 32, 0, stream>>>(N,
				posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 64:
		SoAFunctorN3<floatType, 64> <<<numRequiredBlocks(N), 64, 0, stream>>>(N,
				posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 128:
		SoAFunctorN3<floatType, 128> <<<numRequiredBlocks(N), 128, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 256:
		SoAFunctorN3<floatType, 256> <<<numRequiredBlocks(N), 256, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 512:
		SoAFunctorN3<floatType, 512> <<<numRequiredBlocks(N), 512, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	case 1024:
		SoAFunctorN3<floatType, 1024> <<<numRequiredBlocks(N), 1024, 0, stream>>>(
				N, posX, posY, posZ, forceX, forceY, forceZ);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::SoAFunctorN3Wrapper<float>(int N, float* posX,
		float* posY, float* posZ, float* forceX, float* forceY, float* forceZ,
		cudaStream_t stream);
template void CudaWrapper::SoAFunctorN3Wrapper<double>(int N, double* posX,
		double* posY, double* posZ, double* forceX, double* forceY,
		double* forceZ, cudaStream_t stream);

template<typename floatType>
void CudaWrapper::SoAFunctorN3PairWrapper(int N, floatType* posX,
		floatType* posY, floatType* posZ, floatType* forceX, floatType* forceY,
		floatType* forceZ, int M, floatType* posX2, floatType* posY2,
		floatType* posZ2, floatType* forceX2, floatType* forceY2,
		floatType* forceZ2, cudaStream_t stream) {
	soa<floatType> cell1;
	cell1.posX = posX;
	cell1.posY = posY;
	cell1.posZ = posZ;
	cell1.forceX = forceX;
	cell1.forceY = forceY;
	cell1.forceZ = forceZ;
	soa<floatType> cell2;
	cell2.posX = posX2;
	cell2.posY = posY2;
	cell2.posZ = posZ2;
	cell2.forceX = forceX2;
	cell2.forceY = forceY2;
	cell2.forceZ = forceZ2;
	switch (_num_threads) {
	case 32:
		SoAFunctorN3Pair<floatType, 32> <<<numRequiredBlocks(N), 32, 0, stream>>>(
				N, cell1, M, cell2);
		break;
	case 64:
		SoAFunctorN3Pair<floatType, 64> <<<numRequiredBlocks(N), 64, 0, stream>>>(
				N, cell1, M, cell2);
		break;
	case 128:
		SoAFunctorN3Pair<floatType, 128> <<<numRequiredBlocks(N), 128, 0, stream>>>(
				N, cell1, M, cell2);
		break;
	case 256:
		SoAFunctorN3Pair<floatType, 256> <<<numRequiredBlocks(N), 256, 0, stream>>>(
				N, cell1, M, cell2);
		break;
	case 512:
		SoAFunctorN3Pair<floatType, 512> <<<numRequiredBlocks(N), 512, 0, stream>>>(
				N, cell1, M, cell2);
		break;
	case 1024:
		SoAFunctorN3Pair<floatType, 1024> <<<numRequiredBlocks(N), 1024, 0,
				stream>>>(N, cell1, M, cell2);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();

}

template void CudaWrapper::SoAFunctorN3PairWrapper<float>(int N, float* posX,
		float* posY, float* posZ, float* forceX, float* forceY, float* forceZ,
		int M, float* posX2, float* posY2, float* posZ2, float* forceX2,
		float* forceY2, float* forceZ2, cudaStream_t stream);
template void CudaWrapper::SoAFunctorN3PairWrapper<double>(int N, double* posX,
		double* posY, double* posZ, double* forceX, double* forceY,
		double* forceZ, int M, double* posX2, double* posY2, double* posZ2,
		double* forceX2, double* forceY2, double* forceZ2, cudaStream_t stream);

template<typename floatType, int block_size>
__global__
void LinkedCellsTraversalNoN3(floatType* posX, floatType* posY, floatType* posZ,
		floatType* forceX, floatType* forceY, floatType* forceZ,
		unsigned int* cids, size_t* cellSizes, unsigned int offsets_size,
		int* offsets) {

	int own_cid = cids[blockIdx.x];
	__shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];

	typename vec3<floatType>::Type myposition = { getInfinity<floatType>(),
			getInfinity<floatType>(), getInfinity<floatType>() };
	typename vec3<floatType>::Type myf = { 0, 0, 0 };

	int index = cellSizes[blockIdx.x];
	if (threadIdx.x < cellSizes[own_cid + 1] - cellSizes[own_cid]) {
		myposition.x = posX[index];
		myposition.y = posY[index];
		myposition.z = posZ[index];
	}

	//other cells
	for (auto other_index = 0; other_index < offsets_size; ++other_index) {
		const int other_id = own_cid + offsets[other_index];
		const int cell2Start = cellSizes[other_id];
		const int sizeCell2 = cellSizes[other_id + 1] - cell2Start;

		cell2_pos_shared[threadIdx.x] = {posX[cell2Start+threadIdx.x], posY[cell2Start+threadIdx.x], posZ[cell2Start+threadIdx.x]};
		__syncthreads();
		for (int j = 0; j < sizeCell2; ++j) {
			myf = bodyBodyF<floatType>(myposition, cell2_pos_shared[j], myf);
		}
		__syncthreads();
	}

	if (threadIdx.x < cellSizes[own_cid + 1] - cellSizes[own_cid]) {
		atomicAdd(forceX + index, myf.x);
		atomicAdd(forceY + index, myf.y);
		atomicAdd(forceZ + index, myf.z);
	}
}

template<typename floatType>
void CudaWrapper::LinkedCellsTraversalNoN3Wrapper(floatType* posX,
		floatType* posY, floatType* posZ, floatType* forceX, floatType* forceY,
		floatType* forceZ, unsigned int cids_size, unsigned int* cids,
		unsigned int cellSizes_size, size_t* cellSizes,
		unsigned int offsets_size, int* offsets, cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		LinkedCellsTraversalNoN3<floatType, 32> <<<cids_size, 32, 0, stream>>>(
				posX, posY, posZ, forceX, forceY, forceZ, cids, cellSizes,
				offsets_size, offsets);
		break;
	case 64:
		LinkedCellsTraversalNoN3<floatType, 64> <<<cids_size, 64, 0, stream>>>(
				posX, posY, posZ, forceX, forceY, forceZ, cids, cellSizes,
				offsets_size, offsets);
		break;
	case 96:
		LinkedCellsTraversalNoN3<floatType, 96> <<<cids_size, 96, 0, stream>>>(
				posX, posY, posZ, forceX, forceY, forceZ, cids, cellSizes,
				offsets_size, offsets);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

template void CudaWrapper::LinkedCellsTraversalNoN3Wrapper<float>(float* posX,
		float* posY, float* posZ, float* forceX, float* forceY, float* forceZ,
		unsigned int cids_size, unsigned int* cids, unsigned int cellSizes_size,
		size_t* cellSizes, unsigned int offsets_size, int* offsets,
		cudaStream_t stream);

template void CudaWrapper::LinkedCellsTraversalNoN3Wrapper<double>(double* posX,
		double* posY, double* posZ, double* forceX, double* forceY,
		double* forceZ, unsigned int cids_size, unsigned int* cids,
		unsigned int cellSizes_size, size_t* cellSizes,
		unsigned int offsets_size, int* offsets, cudaStream_t stream);

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
