#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <cudawrappers/nvml.hpp>

#include <ext/alloc_traits.h>

#include "NVMLImpl.h"

namespace pmt::nvml {

NVMLState::operator State() {
  State state(measurements_.size());
  state.timestamp_ = timestamp_;
  state.joules_[0] = joules_ * 1e-3;
  for (size_t i = 0; i < measurements_.size(); i++) {
    state.name_[i] = measurements_[i].name;
    state.watt_[i] = measurements_[i].value * 1e-3;
  }
  return state;
}

NVMLImpl::NVMLImpl(int device_number) {
  const char *pmt_device = getenv("PMT_DEVICE");
  device_number = pmt_device ? atoi(pmt_device) : device_number;

  // Initialize CUDA
  cu::init();
  cu::Device device(device_number);

  // Initialize NVML
  context_ = std::make_unique<::nvml::Context>();
  device_ = std::make_unique<::nvml::Device>(*context_, device);

  // Check whether the CPU+GPU scope is supported (e.g. Grace Hopper)
#if not defined(PMT_NVML_LEGACY_MODE)
  nvmlFieldValue_t values[1];
  values[0].fieldId = kFieldIdPowerAverage;
  values[0].scopeId = 1;
  device_->getFieldValues(1, values);
  nr_scopes_ = 1 + (values[0].nvmlReturn == NVML_SUCCESS);
#endif
}

NVMLImpl::~NVMLImpl() { stopped_ = true; }

#if defined(PMT_NVML_LEGACY_MODE)
std::vector<NVMLMeasurement> NVMLImpl::GetMeasurements() {
  NVMLMeasurement measurement;
  measurement.name = "gpu_average";
  measurement.value = device_->getPower();
  measurement.timestamp = GetTime();
  return {measurement};
}
#else
std::vector<NVMLMeasurement> NVMLImpl::GetMeasurements() {
  const int nr_field_ids = 2;
  const int nr_measurements = nr_scopes_ * nr_field_ids;
  nvmlFieldValue_t values[nr_measurements];
  const unsigned int field_ids[] = {kFieldIdPowerInstant, kFieldIdPowerAverage};

  std::vector<NVMLMeasurement> measurements(nr_measurements);

  for (int i = 0; i < nr_measurements; i += nr_field_ids) {
    const unsigned int scopeId = i / nr_field_ids;
    values[i].fieldId = field_ids[0];
    values[i].scopeId = scopeId;
    values[i + 1].fieldId = field_ids[1];
    values[i + 1].scopeId = scopeId;
  }

  device_->getFieldValues(nr_measurements, values);

  const std::string scopeNames[] = {"gpu", "module"};
  const std::string suffixes[] = {"_instant", "_average"};

  for (int i = 0; i < nr_scopes_; ++i) {
    for (int j = 0; j < nr_field_ids; ++j) {
      int idx = nr_field_ids * i + j;
      measurements[idx].name = scopeNames[i] + suffixes[j];
      measurements[idx].value = values[idx].value.uiVal;
      measurements[idx].timestamp =
          Timestamp(std::chrono::microseconds(values[idx].timestamp));
    }
  }

  return measurements;
}
#endif

NVMLState NVMLImpl::GetNVMLState() {
  if (stopped_) {
    return state_previous_;
  }

  NVMLState state;
  try {
    state.measurements_ = GetMeasurements();

    // Default: use use the instantaneous GPU power
    // Grace Hopper: use the instantaneous module power
#if defined(PMT_NVML_LEGACY_MODE)
    state.watt_ = state.measurements_[0].value;
    state.timestamp_ = state.measurements_[0].timestamp;
#else
    const unsigned int measurement_id = nr_scopes_ == 1 ? 0 : 2;
    state.watt_ = state.measurements_[measurement_id].value;
    state.timestamp_ = state.measurements_[measurement_id].timestamp;
#endif

    // Set derived fields of state
    state.joules_ = state_previous_.joules_;
    const float watt = (state.watt_ + state_previous_.watt_) / 2;
    const double duration =
        seconds(state_previous_.timestamp_, state.timestamp_);
    state.joules_ += watt * duration;
    state.watt_ = watt;

    state_previous_ = state;
  } catch (std::runtime_error &e) {
    return state_previous_;
  }

  return state;
}
}  // end namespace pmt::nvml