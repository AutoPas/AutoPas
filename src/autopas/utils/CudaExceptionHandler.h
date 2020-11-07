/**
 * @file CudaExceptionHandler.h
 * @author jspahl
 * @date 2/10/18
 */

#pragma once

#include "autopas/utils/ExceptionHandler.h"
#if defined(AUTOPAS_CUDA)
#include "cuda_runtime.h"
#endif

namespace autopas {
namespace utils {
/**
 * Handles exceptions in the context of cuda
 */
class CudaExceptionHandler {
 public:
  CudaExceptionHandler() = delete;

#if defined(AUTOPAS_CUDA)
  /**
   * checks Cuda Error Code and throws corresponding autpas exception
   * @param error Cuda Error Code
   */
  static void checkErrorCode(cudaError_t error) {
    if (error == cudaSuccess) return;
    std::string errorname = std::string(cudaGetErrorName(error));
    std::string errorstring = std::string(cudaGetErrorString(error));

    autopas::utils::ExceptionHandler::exception(std::string("cuda error: ") + errorname);
  }

  /**
   * checks cuda Error Code of the most recent cuda call and throws corresponding autpas exception
   */
  static void checkLastCudaCall() {
    cudaError error = cudaGetLastError();
    if (error == cudaSuccess) return;
    std::string errorname = std::string(cudaGetErrorName(error));
    std::string errorstring = std::string(cudaGetErrorString(error));

    autopas::utils::ExceptionHandler::exception(std::string("cuda error: ") + errorname);
  }
#endif
};

}  // namespace utils
}  // namespace autopas
