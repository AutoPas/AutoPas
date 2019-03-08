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

/**
 * This class handles storage of vectors on the gpu.
 * @tparam T types of the elements in the vector
 */
template <typename T>
class CudaDeviceVector {
#if defined(AUTOPAS_CUDA)
 public:
  CudaDeviceVector() : CudaDeviceVector(32) {}

  CudaDeviceVector(size_t max) : _max_size(max), _size(0), _padToMultiple(1024) {
    autopas::utils::CudaExceptionHandler::checkErrorCode(cudaMalloc((void**)&_data, sizeof(T) * _max_size));
  }

  CudaDeviceVector(const CudaDeviceVector<T>& obj) : _max_size(obj._max_size), _size(obj._size) {
    autopas::utils::CudaExceptionHandler::checkErrorCode(cudaMalloc((void**)&_data, sizeof(T) * _max_size));
    cudaMemcpy(_data, obj._data, _size * sizeof(T), cudaMemcpyDeviceToDevice);
  }

  virtual ~CudaDeviceVector() { autopas::utils::CudaExceptionHandler::checkErrorCode(cudaFree(_data)); }

  T* get() { return _data; }

  size_t size() { return _size; }

  void copyHostToDevice(const size_t n, T* hostData, const cudaStream_t stream = 0) {
    _size = n;
    if (n > _max_size) {
      _max_size = (n / _padToMultiple + 1) * _padToMultiple;

      autopas::utils::CudaExceptionHandler::checkErrorCode(cudaFree(_data));
      autopas::utils::CudaExceptionHandler::checkErrorCode(cudaMalloc((void**)&_data, sizeof(T) * _max_size));
    }
    autopas::utils::CudaExceptionHandler::checkErrorCode(
        cudaMemcpyAsync(_data, hostData, n * sizeof(T), cudaMemcpyHostToDevice, stream));
  }

  void copyDeviceToHost(const size_t n, T* hostData, const cudaStream_t stream = 0) {
    autopas::utils::CudaExceptionHandler::checkErrorCode(
        cudaMemcpyAsync(hostData, _data, n * sizeof(T), cudaMemcpyDeviceToHost, stream));
  }

 private:
  size_t _max_size;
  size_t _size;
  size_t _padToMultiple;

  T* _data;

#else
 public:
  /**
   * @brief Dummy default constructor.
   */
  CudaDeviceVector() {}
  /**
   * @brief Dummy construcor.
   * @param max initial maximal vector size
   */
  CudaDeviceVector(size_t max) {}

  /**
   * @brief Dummy Copy Constructor
   * @param obj other object
   */
  CudaDeviceVector(const CudaDeviceVector<T>& obj) {}
#endif
};

}  // namespace utils
}  // namespace autopas
