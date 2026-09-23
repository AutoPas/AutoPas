#include "PowerSensor3Impl.h"

namespace pmt::powersensor3 {

PowerSensor3State::operator State() {
  State state(::PowerSensor3::MAX_PAIRS);
  state.timestamp_ = state_.timeAtRead;

  for (size_t i = 0; i < state.nr_measurements_; i++) {
    state.name_[i] = pair_names_[i];
    state.watt_[i] = state_.voltage[i] * state_.current[i];
    state.joules_[i] = state_.consumedEnergy[i];
  }

  return state;
}

PowerSensor3Impl::PowerSensor3Impl(const char *device)
    : powersensor_(std::make_unique<::PowerSensor3::PowerSensor>(device)) {
  for (unsigned pair_id = 0; pair_id < ::PowerSensor3::MAX_PAIRS; pair_id++) {
    pair_names_.push_back(powersensor_->getPairName(pair_id));
  }
}

PowerSensor3Impl::~PowerSensor3Impl() { stopped_ = true; }

PowerSensor3State PowerSensor3Impl::GetPowerSensor3State() {
  if (stopped_) {
    return state_previous_;
  }

  try {
    const PowerSensor3State state{powersensor_->read(), pair_names_};
    state_previous_ = state;
    return state;
  } catch (std::runtime_error &e) {
    return state_previous_;
  }
}
}  // end namespace pmt::powersensor3
