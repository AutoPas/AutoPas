/**
 * @file LJFunctorCuda.cuh
 *
 * @date 14.12.2018
 * @author jspahl
 */

#pragma once

#include "cuda_runtime.h"
#include "autopas/pairwiseFunctors/FunctorCuda.cuh"

namespace autopas {

/**
 * Stores all pointers to the device Memory SoAs as needed by the LJ Functor
 * @tparam floatType of all vectors
 */
template<typename floatType>
class LJFunctorCudaSoA: public FunctorCudaSoA<floatType> {
public:
	/**
	 * Constructor for only positions
	 * @param size Number of particles
	 * @posX x positions of the particles
	 * @posY y positions of the particles
	 * @posZ z positions of the particles
	 */
	LJFunctorCudaSoA(unsigned int size, floatType* posX, floatType* posY,
			floatType* posZ) :
			_size(size), _posX(posX), _posY(posY), _posZ(posZ), _forceX(NULL), _forceY(
					NULL), _forceZ(NULL) {

	}

	/**
	 * Constructor for only positions
	 * @param size Number of particles
	 * @posX x positions of the particles
	 * @posY y positions of the particles
	 * @posZ z positions of the particles
	 * @forceX x forces of the particles
	 * @forceY y forces of the particles
	 * @forceZ z forces of the particles
	 */
	LJFunctorCudaSoA(unsigned int size, floatType* posX, floatType* posY,
			floatType* posZ, floatType* forceX, floatType* forceY,
			floatType* forceZ) :
			_size(size), _posX(posX), _posY(posY), _posZ(posZ), _forceX(forceX), _forceY(
					forceY), _forceZ(forceZ) {
	}

	/**
	 * CopyConstructor
	 * @param obj other object
	 */
	LJFunctorCudaSoA(const LJFunctorCudaSoA& obj) :
			_size(obj._size), _posX(obj._posX), _posY(obj._posY), _posZ(
					obj._posZ), _forceX(obj._forceX), _forceY(obj._forceY), _forceZ(
					obj._forceZ) {
	}

	unsigned int _size;
	floatType* _posX;
	floatType* _posY;
	floatType* _posZ;
	floatType* _forceX;
	floatType* _forceY;
	floatType* _forceZ;
};

/**
 * Stores all constants needed for the calculation
 * @tparam floatType of constants
 */
template<typename floatType>
class LJFunctorConstants: public FunctorCudaConstants<floatType> {
public:
	LJFunctorConstants() {

	}
	LJFunctorConstants(floatType csq, floatType ep24, floatType sqs,
			floatType sh6) :
			cutoffsquare(csq), epsilon24(ep24), sigmasquare(sqs), shift6(sh6) {
	}
	floatType cutoffsquare;
	floatType epsilon24;
	floatType sigmasquare;
	floatType shift6;
};

/**
 * Wraps vectors of size 3 with the required precision
 * @tparam T floating point Type
 */
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
class LJFunctorCudaWrapper: public CudaWrapperInterface<floatType> {
public:
	LJFunctorCudaWrapper() {
		_num_threads = 32;
	}
	virtual ~LJFunctorCudaWrapper() {

	}

	void setNumThreads(int num_threads) override {
		_num_threads = num_threads;
	}

	void loadConstants(FunctorCudaConstants<floatType>* constants) override;

	void SoAFunctorNoN3Wrapper(FunctorCudaSoA<floatType>* cell1Base,
			cudaStream_t stream = 0) override;
	void SoAFunctorNoN3PairWrapper(FunctorCudaSoA<floatType>* cell1Base,
			FunctorCudaSoA<floatType>* cell2Base, cudaStream_t stream) override;

	void SoAFunctorN3Wrapper(FunctorCudaSoA<floatType>* cell1Base,
			cudaStream_t stream = 0) override;
	void SoAFunctorN3PairWrapper(FunctorCudaSoA<floatType>* cell1Base,
			FunctorCudaSoA<floatType>* cell2Base, cudaStream_t stream) override;

	void LinkedCellsTraversalNoN3Wrapper(FunctorCudaSoA<floatType>* cell1Base,
			unsigned int reqThreads, unsigned int cids_size, unsigned int* cids,
			unsigned int cellSizes_size, size_t* cellSizes,
			unsigned int offsets_size, int* offsets, cudaStream_t stream)
					override;

	void LinkedCellsTraversalN3Wrapper(FunctorCudaSoA<floatType>* cell1Base,
			unsigned int reqThreads, unsigned int cids_size, unsigned int* cids,
			unsigned int cellSizes_size, size_t* cellSizes,
			unsigned int offsets_size, int* offsets, cudaStream_t stream)
					override;

private:
	int numRequiredBlocks(int n) {
		return ((n - 1) / _num_threads) + 1;
	}

	int _num_threads;
};

} // namespace autopas

