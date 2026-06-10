#ifndef PMT_POWERSENSOR3_H_
#define PMT_POWERSENSOR3_H_

#include <memory>
#include <string_view>

#include "common/PMT.h"

namespace pmt::powersensor3 {

class PowerSensor3 : public PMT {
 public:
  constexpr static inline std::string_view name = "powersensor3";
  static std::unique_ptr<PowerSensor3> Create(
      const char *device = default_device().c_str());
  static std::string default_device() { return "/dev/ttyACM0"; }
};

}  // end namespace pmt::powersensor3

#endif  // PMT_POWERSENSOR3_H_
