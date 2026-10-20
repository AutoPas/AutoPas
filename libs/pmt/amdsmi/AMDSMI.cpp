#include "AMDSMI.h"
#include "AMDSMIImpl.h"

namespace pmt::amdsmi {

std::unique_ptr<AMDSMI> AMDSMI::Create(int device_number) {
  return std::unique_ptr<AMDSMI>(new AMDSMIImpl(device_number));
}

}  // end namespace pmt::amdsmi
