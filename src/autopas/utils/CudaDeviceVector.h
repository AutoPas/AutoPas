/**
 * @file CudaDeviceVector.h
 * @author jspahl
 * @date 2/4/18
 */

#pragma once

#include "autopas/utils/CudaExceptionHandler.h"
#include "autopas/utils/ExceptionHandler.h"
#if defined(AUTOPAS_CUDA)
#include "cuda_runtime.h"
#endif

namespace autopas {
namespace utils {

template <typename T>
class CudaDeviceVector {
#if defined(AUTOPAS_CUDA)
 public:
  CudaDeviceVector() : CudaDeviceVector(32) {}

  CudaDeviceVector(size_t max) : _max_size(max), _size(0) {
    autopas::utils::CudaExceptionHandler::checkErrorCode(cudaMalloc((void**)&_data, sizeof(T) * _max_size));
  }

  CudaDeviceVector(const CudaDeviceVector<T>& obj) : _max_size(obj._max_size), _size(obj._size) {
    autopas::utils::CudaExceptionHandler::checkErrorCode(cudaMalloc((void**)&_data, sizeof(T) * _max_size));
    cudaMemcpy(_data, obj._data, _size * sizeof(T), cudaMemcpyDeviceToDevice);
  }

  virtual ~CudaDeviceVector() { autopas::utils::CudaExceptionHandler::checkErrorCode(cudaFree(_data)); }

  T* get() { return _data; }

  size_t size() { return _size; }

  void copyHostToDevice(const size_t n, T* hostData) {
    _size = n;
    if (n > _max_size) {
      _max_size = (n / 32 + 1) * 32;

      autopas::utils::CudaExceptionHandler::checkErrorCode(cudaFree(_data));
      autopas::utils::CudaExceptionHandler::checkErrorCode(cudaMalloc((void**)&_data, sizeof(T) * _max_size));
    }
    autopas::utils::CudaExceptionHandler::checkErrorCode(
        cudaMemcpy(_data, hostData, n * sizeof(T), cudaMemcpyHostToDevice));
  }

  void copyDeviceToHost(const size_t n, T* hostData) {
    autopas::utils::CudaExceptionHandler::checkErrorCode(
        cudaMemcpy(hostData, _data, n * sizeof(T), cudaMemcpyDeviceToHost));
  }

 private:
  size_t _max_size;
  size_t _size;

  T* _data;

#else
 public:
  CudaDeviceVector() {}
  CudaDeviceVector(size_t max) {}

  CudaDeviceVector(const CudaDeviceVector<T>& obj) {}
#endif
};

}  // namespace utils
}  // namespace autopas
