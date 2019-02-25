/**
 * @file LJFunctorCuda.cuh
 *
 * @date 14.12.2018
 * @author jspahl
 */

#pragma once

#include "cuda_runtime.h"

struct constants{
  double cutoffsquare;
  double epsilon24;
  double sigmasquare;
  double shift6;
  double _upotSum;
  double3 _virialSum;
  double3 _lowCorner;
  double3 _highCorner;
};

class CudaWrapper{
public:
	CudaWrapper() {
		_num_threads = 32;
	}
	virtual ~CudaWrapper(){

	}

	void setNumThreads(int num_threads){
		_num_threads=num_threads;
	}
	void loadConstants(double cutoffsquare, double epsilon24, double sigmasquare);

	void AoSFunctorNoN3Wrapper(int N, double* particles);
	void AoSFunctorNoN3PairWrapper(int N, int M, double* particles1, double* particles2);

	void SoAFunctorNoN3Wrapper(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ, cudaStream_t stream = 0);
	void SoAFunctorNoN3PairWrapper(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ,
			int M, double* posX2, double* posY2, double* posZ2, cudaStream_t stream = 0);

	void SoAFunctorN3Wrapper(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ, cudaStream_t stream = 0);
	void SoAFunctorN3PairWrapper(int N, double* posX, double* posY, double* posZ,
			double* forceX, double* forceY, double* forceZ, int M, double* posX2,
			double* posY2, double* posZ2, double* forceX2,
			double* forceY2, double* forceZ2, cudaStream_t stream = 0);
private:
	int numRequiredBlocks(int n){
		return ((n - 1) / _num_threads) + 1;

	}
private:
	int _num_threads;
};

