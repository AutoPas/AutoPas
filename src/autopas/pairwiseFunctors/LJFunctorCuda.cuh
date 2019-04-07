/**
 * @file LJFunctorCuda.cuh
 *
 * @date 14.12.2018
 * @author jspahl
 */

#pragma once

#include "cuda_runtime.h"

namespace autopas {

template<typename floatType>
struct constants {
	floatType cutoffsquare;
	floatType epsilon24;
	floatType sigmasquare;
	floatType shift6;
	floatType _upotSum;
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

class CudaWrapper {
public:
	CudaWrapper() {
		_num_threads = 32;
	}
	virtual ~CudaWrapper() {

	}

	void setNumThreads(int num_threads) {
		_num_threads = num_threads;
	}

	template<typename floatType>
	void loadConstants(floatType cutoffsquare, floatType epsilon24,
			floatType sigmasquare);

	void AoSFunctorNoN3Wrapper(int N, double* particles);
	void AoSFunctorNoN3PairWrapper(int N, int M, double* particles1,
			double* particles2);

	template<typename floatType>
	void SoAFunctorNoN3Wrapper(int N, floatType* posX, floatType* posY,
			floatType* posZ, floatType* forceX, floatType* forceY,
			floatType* forceZ, cudaStream_t stream = 0);
	template<typename floatType>
	void SoAFunctorNoN3PairWrapper(int N, floatType* posX, floatType* posY,
			floatType* posZ, floatType* forceX, floatType* forceY,
			floatType* forceZ, int M, floatType* posX2, floatType* posY2,
			floatType* posZ2, cudaStream_t stream);

	template<typename floatType>
	void SoAFunctorN3Wrapper(int N, floatType* posX, floatType* posY,
			floatType* posZ, floatType* forceX, floatType* forceY,
			floatType* forceZ, cudaStream_t stream = 0);
	template<typename floatType>
	void SoAFunctorN3PairWrapper(int N, floatType* posX, floatType* posY,
			floatType* posZ, floatType* forceX, floatType* forceY,
			floatType* forceZ, int M, floatType* posX2, floatType* posY2,
			floatType* posZ2, floatType* forceX2, floatType* forceY2,
			floatType* forceZ2, cudaStream_t stream);

	template<typename floatType>
	void LinkedCellsTraversalNoN3Wrapper(floatType* posX, floatType* posY,
			floatType* posZ, floatType* forceX, floatType* forceY,
			floatType* forceZ, unsigned int cids_size, unsigned int* cids,
			unsigned int cellSizes_size, size_t* cellSizes,
			unsigned int offsets_size, int* offsets, cudaStream_t stream);

	template<typename floatType>
	void LinkedCellsTraversalN3Wrapper(floatType* posX, floatType* posY,
			floatType* posZ, floatType* forceX, floatType* forceY,
			floatType* forceZ, unsigned int cids_size, unsigned int* cids,
			unsigned int cellSizes_size, size_t* cellSizes,
			unsigned int offsets_size, int* offsets, cudaStream_t stream);

	template<typename floatType>
	void CellVerletTraversalNoN3Wrapper(floatType* posX, floatType* posY,
			floatType* posZ, floatType* forceX, floatType* forceY,
			floatType* forceZ, unsigned int ncells, unsigned int others_size,
			unsigned int* other_ids, cudaStream_t stream);

	template<typename floatType>
	void CellVerletTraversalN3Wrapper(floatType* posX, floatType* posY,
			floatType* posZ, floatType* forceX, floatType* forceY,
			floatType* forceZ, unsigned int ncells, unsigned int others_size,
			unsigned int* other_ids, cudaStream_t stream);

private:
	int numRequiredBlocks(int n) {
		return ((n - 1) / _num_threads) + 1;
	}

	int _num_threads;
};

} // namespace autopas

