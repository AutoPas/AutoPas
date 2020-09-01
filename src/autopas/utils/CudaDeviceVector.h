/**
 * @file CudaDeviceVector.h
 * @author jspahl
 * @date 2/4/18
 */

#pragma once

#include "autopas/utils/CudaExceptionHandler.h"
#if defined(AUTOPAS_CUDA)
#include "cuda_runtime.h"
#endif

namespace autopas::utils {

/**
 * This class handles storage of vectors on the gpu.
 * @tparam T types of the elements in the vector
 */
template <typename T>
class CudaDeviceVector {
#if defined(AUTOPAS_CUDA)
 public:
  /**
   * @brief Construct vector without allocating memory
   */
  CudaDeviceVector() : _max_size(0), _size(0), _padToMultiple(1024) {}

  /**
   * @brief Construct vector with custom maximal length and allocate the Memory on the GPU
   * @param max size to be allocated
   */
  CudaDeviceVector(size_t max) : _max_size(max), _size(0), _padToMultiple(1024) {
    autopas::utils::CudaExceptionHandler::checkErrorCode(cudaMalloc((void **)&_data, sizeof(T) * _max_size));
  }

  /**
   * @brief Copy Constructor
   * @param obj other object
   */
  CudaDeviceVector(const CudaDeviceVector<T> &obj) : _max_size(obj._max_size), _size(obj._size) {
    if (obj._max_size > 0) {
      autopas::utils::CudaExceptionHandler::checkErrorCode(cudaMalloc((void **)&_data, sizeof(T) * _max_size));
      cudaMemcpy(_data, obj._data, _size * sizeof(T), cudaMemcpyDeviceToDevice);
    }
  }

  /**
   * @brief Frees Memory on GPU
   */
  virtual ~CudaDeviceVector() {
    if (_max_size > 0) autopas::utils::CudaExceptionHandler::checkErrorCode(cudaFree(_data));
  }

  /**
   * @brief return Pointer to the data in the GPU Memory
   * @return Pointer to the data in the GPU Memory
   */
  T *get() { return _data; }

  /**
   * @brief return number of elemnts in the vector
   * @return number of elemnts in the vector
   */
  size_t size() { return _size; }

  /**
   * @brief copy data from host to this vector on the device
   * @param n number of elements
   * @param hostData Pointer to host data
   * @param stream Cuda Stream to use for the copy
   */
  void copyHostToDevice(const size_t n, T *hostData, const cudaStream_t stream = 0) {
    _size = n;
    if (n > _max_size) {
      if (_max_size > 0) autopas::utils::CudaExceptionHandler::checkErrorCode(cudaFree(_data));
      _max_size = (n / _padToMultiple + 1) * _padToMultiple;

      autopas::utils::CudaExceptionHandler::checkErrorCode(cudaMalloc((void **)&_data, sizeof(T) * _max_size));
    }
    autopas::utils::CudaExceptionHandler::checkErrorCode(
        cudaMemcpyAsync(_data, hostData, n * sizeof(T), cudaMemcpyHostToDevice, stream));
  }

  /**
   * @brief copy data from this vectors device Memory back to the host
   * @param n number of elements
   * @param hostData Pointer to host data
   * @param stream Cuda Stream to use for the copy
   */
  void copyDeviceToHost(const size_t n, T *hostData, const cudaStream_t stream = 0) {
    autopas::utils::CudaExceptionHandler::checkErrorCode(
        cudaMemcpyAsync(hostData, _data, n * sizeof(T), cudaMemcpyDeviceToHost, stream));
  }

 private:
  size_t _max_size;
  size_t _size;
  size_t _padToMultiple;

  T *_data;

#else
 public:
  /**
   * @brief Dummy default constructor.
   */
  CudaDeviceVector() {}
  /**
   * @brief Dummy constructor.
   * @param max initial maximal vector size
   */
  CudaDeviceVector(size_t max) {}

  /**
   * @brief Dummy Copy Constructor
   * @param obj other object
   */
  CudaDeviceVector(const CudaDeviceVector<T> &obj) {}
#endif
};

}  // namespace autopas::utils
