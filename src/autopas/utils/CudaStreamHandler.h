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

namespace autopas {
namespace utils {

class CudaStreamHandler {
#if defined(AUTOPAS_CUDA)
 public:
  CudaStreamHandler(size_t nStreams) : _streams(nStreams), random_index(0) {
    for (int i = 0; i < nStreams; ++i) {
      cudaStreamCreate(_streams.data() + i);
    }
  }

  virtual ~CudaStreamHandler() {
    for (auto it : _streams) {
      cudaStreamDestroy(it);
    }
  }

  cudaStream_t& getStream(int index) { return _streams.at(index); }
  cudaStream_t& getStreambyHash(size_t hash) { return _streams[hash % _streams.size()]; }
  cudaStream_t& getStreamRandom() { return _streams[(++random_index) % _streams.size()]; }

 private:
  std::vector<cudaStream_t> _streams;
  size_t random_index;
#else
 public:
  CudaStreamHandler() {}
#endif
};

}  // namespace utils
}  // namespace autopas
