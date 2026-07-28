#include <vector>

#include <ext/alloc_traits.h>

#include "PowerSensor2.h"

#include <PowerSensor.h>

namespace {
double Seconds(const PowerSensor::State &first,
               const PowerSensor::State &second) {
  return PowerSensor::seconds(first, second);
}

double Joules(const PowerSensor::State &first, const PowerSensor::State &second,
              int sensor_id) {
  return PowerSensor::Joules(first, second, sensor_id);
}

double Watt(const PowerSensor::State &first, const PowerSensor::State &second,
            int sensor_id) {
  return PowerSensor::Watt(first, second, sensor_id);
}
}  // namespace

namespace pmt::powersensor2 {

template <typename PowerSensor, typename PowerSensorState>
class PowerSensor2Impl : public PowerSensor2 {
 public:
  PowerSensor2Impl(const char *device)
      : powersensor_(std::make_unique<PowerSensor>(device)),
        first_state_(powersensor_->read()) {}

  State GetState() override {
    const PowerSensorState powersensor_state = powersensor_->read();
    State state;
    state.timestamp_ =
        Timestamp(std::chrono::microseconds(static_cast<long long>(
            ::Seconds(first_state_, powersensor_state) * 1e9)));
    state.name_[0] = "device";
    state.joules_[0] = ::Joules(first_state_, powersensor_state, -1);
    state.watt_[0] = ::Watt(first_state_, powersensor_state, -1);
    return state;
  }

 private:
  virtual const char *GetDumpFilename() override {
    return "/tmp/pmt_powersensor2.out";
  }

  std::unique_ptr<PowerSensor> powersensor_{};
  PowerSensorState first_state_{};
};

std::unique_ptr<PowerSensor2> PowerSensor2::Create(const char *device) {
  return std::make_unique<
      PowerSensor2Impl<PowerSensor::PowerSensor, PowerSensor::State>>(device);
}

}  // end namespace pmt::powersensor2
