#include "Xilinx.h"

#include <istream>
#include <stdexcept>
#include <vector>

#include <errno.h>
#include <ext/alloc_traits.h>
#include <stdlib.h>

namespace anonymous {
float GetPower(std::string &filename) {
  // Open power file, e.g.
  // /sys/devices/pci0000:a0/0000:a0:03.1/0000:a1:00.0/hwmon/hwmon3/power1_input
  std::ifstream file(filename, std::ios::in | std::ios::binary);
  if (errno != 0) {
    throw std::runtime_error("Could not open: " + filename);
  }

  // This file has one line with instantenous power consumption in uW
  size_t power;
  file >> power;
  return power;
}

}  // namespace anonymous
namespace pmt::xilinx {

class XilinxState {
 public:
  operator State();
  Timestamp timestamp_;
  double watt_ = 0;
  double joules_ = 0;
};

class XilinxImpl : public Xilinx {
 public:
  XilinxImpl(const char *device);

 private:
  State GetState() override { return GetXilinxState(); }

  virtual const char *GetDumpFilename() override {
    return "/tmp/pmt_xilinx.out";
  }

  std::string filename_;

  XilinxState previous_state_;
  XilinxState GetXilinxState();
};

XilinxState::operator State() {
  State state;
  state.timestamp_ = timestamp_;
  state.name_[0] = "device";
  state.joules_[0] = joules_;
  state.watt_[0] = watt_;
  return state;
}

std::unique_ptr<Xilinx> Xilinx::Create(const char *device) {
  return std::make_unique<XilinxImpl>(device);
}

XilinxImpl::XilinxImpl(const char *device) {
  char *pmt_device = getenv("PMT_DEVICE");
  filename_ = pmt_device ? pmt_device : device;

  previous_state_ = GetXilinxState();
  previous_state_.joules_ = 0;
}

XilinxState XilinxImpl::GetXilinxState() {
  XilinxState state;
  state.timestamp_ = GetTime();
  state.watt_ = anonymous::GetPower(filename_) * 1e-6;
  state.joules_ = previous_state_.joules_;
  const float watts = (state.watt_ + previous_state_.watt_) / 2;
  const double duration = seconds(previous_state_.timestamp_, state.timestamp_);
  state.joules_ += watts * duration;
  previous_state_ = state;
  return state;
}

}  // end namespace pmt::xilinx
