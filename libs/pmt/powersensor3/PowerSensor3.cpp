#include "PowerSensor3.h"
#include "PowerSensor3Impl.h"

namespace pmt::powersensor3 {

std::unique_ptr<PowerSensor3> PowerSensor3::Create(const char *device) {
  return std::unique_ptr<PowerSensor3>(new PowerSensor3Impl(device));
}

}  // end namespace pmt::powersensor3
