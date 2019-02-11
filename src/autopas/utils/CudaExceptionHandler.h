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

class CudaExceptionHandler {
 public:
  CudaExceptionHandler() = delete;

  static void checkErrorCode(cudaError_t error) {
#if defined(AUTOPAS_CUDA)
    if (error == 0) return;
    std::string errorname = std::string(cudaGetErrorName(error));
    std::string errorstring = std::string(cudaGetErrorString(error));
#endif

    autopas::utils::ExceptionHandler::exception(std::string("cuda error") + errorname);
  }
};

}  // namespace utils
}  // namespace autopas
