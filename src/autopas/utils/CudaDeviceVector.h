/**
 * @file CudaDeviceVector.h
 * @author jspahl
 * @date 2/4/18
 */

#pragma once

#include "cuda_runtime.h"

namespace autopas {

template <typename T>
class CudaDeviceVector{
	public:
	CudaDeviceVector(size_t max):_max_size(max), _size(0){
		cudaMalloc((void **)&_data, sizeof(T) * _max_size);
	}

	virtual ~CudaDeviceVector(){
		cudaFree(_data);
	}

	T* get(){
		return _data;
	}

	size_t size(){
		return _size;
	}

	cudaError_t copyHostToDevice(int n, T* hostData){
		_size = n;
		if (n > _max_size){
			_max_size = (n/32 + 1) * 32;

			cudaFree(_data);
			cudaMalloc((void **)&_data, sizeof(T) * _max_size);
		}
		return cudaMemcpy(_data, hostData, n, cudaMemcpyHostToDevice);
	}

	cudaError_t copyDeviceToHost(int n, T* hostData){
		return cudaMemcpy(hostData, _data, _size, cudaMemcpyDeviceToHost);
	}

	private:
	size_t _max_size;
	size_t _size;

	T* _data;
};

}  // namespace autopas
