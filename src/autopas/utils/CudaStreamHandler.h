/**
 * @file CudaStreamHandler.h
 * @author jspahl
 * @date 3/1/19
 */

#pragma once

#include "autopas/utils/CudaExceptionHandler.h"
#include "autopas/utils/ExceptionHandler.h"
#include "vector"
#if defined(AUTOPAS_CUDA)
#include "cuda_runtime.h"
#endif

namespace autopas::utils {

/**
 * Handles an Array of Cuda streams and provides different access Methods.
 */
class CudaStreamHandler {
#if defined(AUTOPAS_CUDA)
 public:
  /**
   * Creates a CudaStreamHandler object with nStreams different cuda streams
   * @param nStreams number of streams in Handler
   */
  CudaStreamHandler(size_t nStreams) : _streams(nStreams), random_index(0) {
    for (size_t i = 0; i < nStreams; ++i) {
      cudaStreamCreate(_streams.data() + i);
    }
  }

  /**
   * Deallocates all cuda streams
   */
  virtual ~CudaStreamHandler() {
    for (auto it : _streams) {
      cudaStreamDestroy(it);
    }
  }
  /**
   * Returns cuda stream at the given index
   * @param index
   * @return corresponding cuda stream
   */
  cudaStream_t &getStream(int index) { return _streams.at(index); }
  /**
   * Returns cuda stream at the hash value modulo the number of streams
   * @param hash
   * @return corresponding cuda stream
   */
  cudaStream_t &getStreambyHash(size_t hash) { return _streams[hash % _streams.size()]; }
  /**
   * Returns cuda stream at the least recently used index
   * @return ranodm cuda stream
   */
  cudaStream_t &getStreamRandom() { return _streams[(++random_index) % _streams.size()]; }

 private:
  std::vector<cudaStream_t> _streams;
  size_t random_index;
#else
 public:
  /**
   * Empty Dummy default Constructor in case there is no cuda support
   */
  CudaStreamHandler() = default;
#endif
};

}  // namespace autopas::utils
