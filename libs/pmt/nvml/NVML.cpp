#include "NVML.h"
#include "NVMLImpl.h"

namespace pmt::nvml {

std::unique_ptr<NVML> NVML::Create(int device_number) {
  return std::unique_ptr<NVML>(new NVMLImpl(device_number));
}

}  // end namespace pmt::nvml