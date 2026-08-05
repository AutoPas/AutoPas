#include <sstream>

#include <cudawrappers/cu.hpp>

#include "NVIDIA.h"
#if defined(PMT_BUILD_NVML)
#include "nvml/NVML.h"
#endif
#if defined(PMT_BUILD_TEGRA)
#include "tegra/Tegra.h"
#endif

#if defined(PMT_BUILD_TEGRA)
#define checkCudaCall(val) __checkCudaCall((val), #val, __FILE__, __LINE__)

inline void __checkCudaCall(cudaError_t result, const char *const func,
                            const char *const file, int const line) {
  if (result != cudaSuccess) {
    std::stringstream error;
    error << "CUDA Error at " << file;
    error << ":" << line;
    error << " in function " << func;
    error << ": " << cudaGetErrorString(result);
    error << std::endl;
    throw std::runtime_error(error.str());
  }
}
#endif

namespace pmt::nvidia {

std::unique_ptr<PMT> NVIDIA::Create(int device_number) {
#if defined(PMT_BUILD_TEGRA)
  cu::init();
  cu::Device device(device_number);
  if (device.getAttribute(CU_DEVICE_ATTRIBUTE_INTEGRATED)) {
    return tegra::Tegra::Create();
  }
#endif
#if defined(PMT_BUILD_NVML)
  return nvml::NVML::Create(device_number);
#endif

  throw std::runtime_error("Neither Tegra nor NVML are available.");
}

}  //  end namespace pmt::nvidia