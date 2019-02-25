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
#include "autopas/utils/CudaDeviceVector.h"
#include "autopas/utils/CudaSoA.h"

__constant__ constants global_constants;

template<typename floating_precision = double>
class soa {
public:
	floating_precision* posX;
	floating_precision* posY;
	floating_precision* posZ;
	floating_precision* forceX;
	floating_precision* forceY;
	floating_precision* forceZ;
};

template<typename T> struct vec3 {
	typedef T Type;
};
template<> struct vec3<float> {
	typedef float3 Type;
};
template<> struct vec3<double> {
	typedef double3 Type;
};

template<typename floatType>
__device__
typename vec3<floatType>::Type bodyBodyF(typename vec3<floatType>::Type i,
		typename vec3<floatType>::Type j, typename vec3<floatType>::Type fi) {
	auto drx = i.x - j.x;
	auto dry = i.y - j.y;
	auto drz = i.z - j.z;

	auto dr2 = drx * drx + dry * dry + drz * drz;

	if (dr2 > global_constants.cutoffsquare | dr2 == 0.0) {
		return fi;
	}

	auto invdr2 = 1. / dr2;
	auto lj6 = global_constants.sigmasquare * invdr2;
	lj6 = lj6 * lj6 * lj6;
	auto lj12 = lj6 * lj6;
	auto lj12m6 = lj12 - lj6;
	auto fac = global_constants.epsilon24 * (lj12 + lj12m6) * invdr2;

	fi.x += drx * fac;
	fi.y += dry * fac;
	fi.z += drz * fac;

	return fi;
}

template<typename floatType>
__device__
typename vec3<floatType>::Type bodyBodyFN3(typename vec3<floatType>::Type i,
		typename vec3<floatType>::Type j, typename vec3<floatType>::Type fi,
		typename vec3<floatType>::Type* fj) {
	auto drx = i.x - j.x;
	auto dry = i.y - j.y;
	auto drz = i.z - j.z;

	auto dr2 = drx * drx + dry * dry + drz * drz;

	if (dr2 > global_constants.cutoffsquare) {
		return fi;
	}

	auto invdr2 = 1. / dr2;
	auto lj6 = global_constants.sigmasquare * invdr2;
	lj6 = lj6 * lj6 * lj6;
	auto lj12 = lj6 * lj6;
	auto lj12m6 = lj12 - lj6;
	auto fac = global_constants.epsilon24 * (lj12 + lj12m6) * invdr2;

	auto dfx = drx * fac;
	auto dfy = dry * fac;
	auto dfz = drz * fac;

	fi.x += dfx;
	fi.y += dfy;
	fi.z += dfz;

	atomicAdd(&(fj->x), -dfx);
	atomicAdd(&(fj->y), -dfy);
	atomicAdd(&(fj->z), -dfz);

	return fi;
}

__device__ double3 bodyBodyFcalcGlobals(double3 i, double3 j, double3 fi,
		double4 globals) {
	double drx = i.x - j.x;
	double dry = i.y - j.y;
	double drz = i.z - j.z;

	double dr2 = drx * drx + dry * dry + drz * drz;

	if (dr2 > global_constants.cutoffsquare | dr2 == 0.0) {
		return fi;
	}

	double invdr2 = 1. / dr2;
	double lj6 = global_constants.sigmasquare * invdr2;
	lj6 = lj6 * lj6 * lj6;
	double lj12 = lj6 * lj6;
	double lj12m6 = lj12 - lj6;
	double fac = global_constants.epsilon24 * (lj12 + lj12m6) * invdr2;

	const double fx = drx * fac;
	const double fy = dry * fac;
	const double fz = drz * fac;

	const double virialx = drx * fx;
	const double virialy = dry * fy;
	const double virialz = drz * fz;
	const double upot = (global_constants.epsilon24 * lj12m6
			+ global_constants.shift6);

	// these calculations assume that this functor is not called for halo cells!
	globals.w += upot;
	globals.x += virialx;
	globals.y += virialy;
	globals.z += virialz;

	fi.x += fx;
	fi.y += fy;
	fi.z += fz;

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
				myf = bodyBodyF<double>(myposition, block_pos[j], myf);
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
void SoAFunctorN3(int N, double* posX, double* posY, double* posZ,
		double* forceX, double* forceY, double* forceZ) {
	static_assert((block_size & (block_size - 1)) == 0, "block size must be power of 2");
	__shared__ typename vec3<floatType>::Type cell1_pos_shared[block_size];
	__shared__ typename vec3<floatType>::Type cell1_forces_shared[block_size];
	int tid = blockIdx.x * block_size + threadIdx.x;
	typename vec3<floatType>::Type myposition = { CUDART_INF, CUDART_INF,
			CUDART_INF };
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

		for (int j = threadIdx.x -1; j >= 0; --j) {
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
void SoAFunctorN3Pair(int N, soa<> cell1, int M, soa<> cell2) {
	static_assert((block_size & (block_size - 1)) == 0, "block size must be power of 2");
	__shared__ typename vec3<floatType>::Type cell2_pos_shared[block_size];
	__shared__ typename vec3<floatType>::Type cell2_forces_shared[block_size];
	int tid = blockIdx.x * block_size + threadIdx.x;
	typename vec3<floatType>::Type myposition = { CUDART_INF, CUDART_INF,
			CUDART_INF };
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
			myf = bodyBodyFN3<floatType>(myposition, cell2_pos_shared[offset],
					myf, cell2_forces_shared + offset);
		}
		__syncthreads();

		atomicAdd(cell2.forceX + idx, cell2_forces_shared[threadIdx.x].x);
		atomicAdd(cell2.forceY + idx, cell2_forces_shared[threadIdx.x].y);
		atomicAdd(cell2.forceZ + idx, cell2_forces_shared[threadIdx.x].z);
		__syncthreads();
	}
	if (not NMisMultipleBlockSize && i > M) {
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

void CudaWrapper::SoAFunctorNoN3Wrapper(int N, double* posX, double* posY,
		double* posZ, double* forceX, double* forceY, double* forceZ,
		cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		SoAFunctorNoN3<double, 32> <<<numRequiredBlocks(N), 32>>>(N, posX, posY,
				posZ, forceX, forceY, forceZ);
		break;
	case 64:
		SoAFunctorNoN3<double, 64> <<<numRequiredBlocks(N), 64>>>(N, posX, posY,
				posZ, forceX, forceY, forceZ);
		break;
	case 96:
		SoAFunctorNoN3<double, 96> <<<numRequiredBlocks(N), 96>>>(N, posX, posY,
				posZ, forceX, forceY, forceZ);
		break;
	case 128:
		SoAFunctorNoN3<double, 128> <<<numRequiredBlocks(N), 128>>>(N, posX,
				posY, posZ, forceX, forceY, forceZ);
		break;
	case 256:
		SoAFunctorNoN3<double, 256> <<<numRequiredBlocks(N), 256>>>(N, posX,
				posY, posZ, forceX, forceY, forceZ);
		break;
	case 512:
		SoAFunctorNoN3<double, 512> <<<numRequiredBlocks(N), 512>>>(N, posX,
				posY, posZ, forceX, forceY, forceZ);
		break;
	case 1024:
		SoAFunctorNoN3<double, 1024> <<<numRequiredBlocks(N), 1024>>>(N, posX,
				posY, posZ, forceX, forceY, forceZ);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

void CudaWrapper::SoAFunctorNoN3PairWrapper(int N, double* posX, double* posY,
		double* posZ, double* forceX, double* forceY, double* forceZ, int M,
		double* posX2, double* posY2, double* posZ2, cudaStream_t stream) {
	switch (_num_threads) {
	case 32:
		SoAFunctorNoN3Pair<double, 32> <<<numRequiredBlocks(N), 32>>>(N, posX,
				posY, posZ, forceX, forceY, forceZ, M, posX2, posY2, posZ2);
		break;
	case 64:
		SoAFunctorNoN3Pair<double, 64> <<<numRequiredBlocks(N), 64>>>(N, posX,
				posY, posZ, forceX, forceY, forceZ, M, posX2, posY2, posZ2);
		break;
	case 96:
		SoAFunctorNoN3Pair<double, 96> <<<numRequiredBlocks(N), 96>>>(N, posX,
				posY, posZ, forceX, forceY, forceZ, M, posX2, posY2, posZ2);
		break;
	case 128:
		SoAFunctorNoN3Pair<double, 128> <<<numRequiredBlocks(N), 128>>>(N, posX,
				posY, posZ, forceX, forceY, forceZ, M, posX2, posY2, posZ2);
		break;
	case 256:
		SoAFunctorNoN3Pair<double, 256> <<<numRequiredBlocks(N), 256>>>(N, posX,
				posY, posZ, forceX, forceY, forceZ, M, posX2, posY2, posZ2);
		break;
	case 512:
		SoAFunctorNoN3Pair<double, 512> <<<numRequiredBlocks(N), 512>>>(N, posX,
				posY, posZ, forceX, forceY, forceZ, M, posX2, posY2, posZ2);
		break;
	case 1024:
		SoAFunctorNoN3Pair<double, 1024> <<<numRequiredBlocks(N), 1024>>>(N,
				posX, posY, posZ, forceX, forceY, forceZ, M, posX2, posY2,
				posZ2);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

void CudaWrapper::SoAFunctorN3Wrapper(int N, double* posX, double* posY,
		double* posZ, double* forceX, double* forceY, double* forceZ) {

	switch (_num_threads) {
	case 32:
		SoAFunctorN3<double, 32> <<<numRequiredBlocks(N), 32>>>(N, posX, posY,
				posZ, forceX, forceY, forceZ);
		break;
	case 64:
		SoAFunctorN3<double, 64> <<<numRequiredBlocks(N), 64>>>(N, posX, posY,
				posZ, forceX, forceY, forceZ);
		break;
	case 128:
		SoAFunctorN3<double, 128> <<<numRequiredBlocks(N), 128>>>(N, posX, posY,
				posZ, forceX, forceY, forceZ);
		break;
	case 256:
		SoAFunctorN3<double, 256> <<<numRequiredBlocks(N), 256>>>(N, posX, posY,
				posZ, forceX, forceY, forceZ);
		break;
	case 512:
		SoAFunctorN3<double, 512> <<<numRequiredBlocks(N), 512>>>(N, posX, posY,
				posZ, forceX, forceY, forceZ);
		break;
	case 1024:
		SoAFunctorN3<double, 1024> <<<numRequiredBlocks(N), 1024>>>(N, posX, posY,
				posZ, forceX, forceY, forceZ);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

void CudaWrapper::SoAFunctorN3PairWrapper(int N, double* posX, double* posY,
		double* posZ, double* forceX, double* forceY, double* forceZ, int M,
		double* posX2, double* posY2, double* posZ2, double* forceX2,
		double* forceY2, double* forceZ2) {
	soa<double> cell1;
	cell1.posX = posX;
	cell1.posY = posY;
	cell1.posZ = posZ;
	cell1.forceX = forceX;
	cell1.forceY = forceY;
	cell1.forceZ = forceZ;
	soa<double> cell2;
	cell2.posX = posX2;
	cell2.posY = posY2;
	cell2.posZ = posZ2;
	cell2.forceX = forceX2;
	cell2.forceY = forceY2;
	cell2.forceZ = forceZ2;
	switch (_num_threads) {
	case 32:
		SoAFunctorN3Pair<double, 32> <<<numRequiredBlocks(N), 32>>>(N, cell1, M,
				cell2);
		break;
	case 64:
		SoAFunctorN3Pair<double, 64> <<<numRequiredBlocks(N), 64>>>(N, cell1, M,
				cell2);
		break;
	case 128:
		SoAFunctorN3Pair<double, 128> <<<numRequiredBlocks(N), 128>>>(N, cell1,
				M, cell2);
		break;
	case 256:
		SoAFunctorN3Pair<double, 256> <<<numRequiredBlocks(N), 256>>>(N, cell1,
				M, cell2);
		break;
	case 512:
		SoAFunctorN3Pair<double, 512> <<<numRequiredBlocks(N), 512>>>(N, cell1,
				M, cell2);
		break;
	case 1024:
		SoAFunctorN3Pair<double, 1024> <<<numRequiredBlocks(N), 1024>>>(N,
				cell1, M, cell2);
		break;
	default:
		autopas::utils::ExceptionHandler::exception(
				std::string("cuda Kernel size not available"));
		break;
	}
	autopas::utils::CudaExceptionHandler::checkLastCudaCall();
}

void CudaWrapper::loadConstants(double cutoffsquare, double epsilon24,
		double sigmasquare) {

	constants c;
	c.cutoffsquare = cutoffsquare;
	c.epsilon24 = epsilon24;
	c.sigmasquare = sigmasquare;

	cudaMemcpyToSymbol(global_constants, &c, sizeof(constants));
}
