#include <algorithm>
#include <exception>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <vector>

#include <ext/alloc_traits.h>

#include "Cray.h"
#include "FilenamesHelper.h"

namespace {
double GetPower(const std::string& filePath) {
  std::ifstream powerFile(filePath);
  if (!powerFile.is_open()) {
    throw std::runtime_error("Failed to open power file");
  }

  std::string line;
  if (!std::getline(powerFile, line)) {
    throw std::runtime_error("Failed to read power value");
  }

  double power;
  try {
    power = std::stod(line);
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to parse power value");
  }

  powerFile.close();
  return power;
}
}  // namespace

namespace pmt::cray {

class CrayImpl : public Cray {
 public:
  CrayImpl();

  State GetState() override;

 private:
  std::vector<std::string> filenames;
  std::string cray_pm_counters_path = "/sys/cray/pm_counters";

  Timestamp previous_timestamp_;
  std::vector<CrayMeasurement> previous_measurements_;

  // Mutex used to guard GetMeasurements()
  std::mutex mutex_;

  virtual const char* GetDumpFilename() override { return "/tmp/pmt_cray.out"; }

  unsigned int device_number_;

  std::vector<CrayMeasurement> GetMeasurements();
};

std::unique_ptr<Cray> Cray::Create() { return std::make_unique<CrayImpl>(); }

CrayImpl::CrayImpl() {
  filenames = filenames_helper::GetFilenames(cray_pm_counters_path);
  filenames = filenames_helper::ReorderFilenames(filenames);
#if defined(DEBUG)
  filenames_helper::PrintFilenames(filenames);
#endif

  previous_timestamp_ = GetTime();
  previous_measurements_ = GetMeasurements();
}

std::vector<CrayMeasurement> CrayImpl::GetMeasurements() {
  std::lock_guard<std::mutex> lock(mutex_);

  std::vector<CrayMeasurement> measurements;

  std::string file_path = "";

  for (const auto& filename : filenames) {
    file_path = cray_pm_counters_path + "/" + filename;
    CrayMeasurement measurement;
    measurement.name = filename;
    measurement.watt = GetPower(file_path);
    measurements.push_back(measurement);
  }

  return measurements;
}

State CrayImpl::GetState() {
  std::vector<CrayMeasurement> measurements = GetMeasurements();
  State state(measurements.size());
  state.timestamp_ = GetTime();

  const double duration = seconds(previous_timestamp_, state.timestamp_);

  for (size_t i = 0; i < measurements.size(); i++) {
    state.name_[i] = measurements[i].name;
    state.watt_[i] = measurements[i].watt;
    const double watt =
        (measurements[i].watt + previous_measurements_[i].watt) / 2;
    state.joules_[i] += watt * duration;
  }

  return state;
}

}  // end namespace pmt::cray
