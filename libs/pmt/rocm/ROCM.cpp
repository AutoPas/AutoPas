#include "ROCM.h"
#include "ROCMImpl.h"

namespace pmt::rocm {

std::unique_ptr<ROCM> ROCM::Create(int device_number) {
  return std::make_unique<ROCMImpl>(device_number);
}

}  // end namespace pmt::rocm