#ifndef PMT_POWERSENSOR2_H_
#define PMT_POWERSENSOR2_H_

#include <memory>
#include <string_view>

#include "common/PMT.h"

namespace pmt::powersensor2 {

class PowerSensor2 : public PMT {
 public:
  constexpr static inline std::string_view name = "powersensor2";
  static std::unique_ptr<PowerSensor2> Create(
      const char *device = default_device().c_str());
  static std::string default_device() { return "/dev/ttyACM0"; }
};

}  // end namespace pmt::powersensor2

#endif  // PMT_POWERSENSOR2_H_
