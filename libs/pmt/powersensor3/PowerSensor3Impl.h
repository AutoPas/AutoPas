#include <memory>
#include <string_view>

#include <PowerSensor.hpp>

#include "common/PMT.h"
#include "PowerSensor3.h"

namespace pmt::powersensor3 {

class PowerSensor3State {
 public:
  operator State();
  ::PowerSensor3::State state_;
  std::vector<std::string> pair_names_;
};

class PowerSensor3Impl : public PowerSensor3 {
 public:
  PowerSensor3Impl(const char *device);
  ~PowerSensor3Impl();

 private:
  State GetState() override { return GetPowerSensor3State(); }

  virtual const char *GetDumpFilename() override {
    return "/tmp/pmt_powersensor3.out";
  }

  PowerSensor3State state_previous_;
  PowerSensor3State GetPowerSensor3State();

  bool stopped_ = false;

  std::unique_ptr<::PowerSensor3::PowerSensor> powersensor_;
  std::vector<std::string> pair_names_;
};

}  // end namespace pmt::powersensor3