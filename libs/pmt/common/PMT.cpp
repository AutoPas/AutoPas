#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>

#include <pmt.h>

#include "common/PMT.h"
namespace {
template <typename T>
bool isEqual(T x, T y) {
  return std::fabs(x - y) <= std::numeric_limits<T>::epsilon();
}

bool isNumber(const std::string &s) {
  return !s.empty() && std::all_of(s.begin(), s.end(),
                                   [](char c) { return std::isdigit(c); });
}
}  // end namespace

namespace pmt {

PMT::~PMT() {
};

double PMT::seconds(const Timestamp &timestamp) {
  return std::chrono::duration_cast<std::chrono::microseconds>(
             timestamp.time_since_epoch())
             .count() /
         1.e6;
}

double PMT::seconds(const Timestamp &first, const Timestamp &second) {
  return std::chrono::duration_cast<std::chrono::microseconds>((second - first))
             .count() /
         1.e6;
}

double PMT::seconds(const State &first, const State &second) {
  return seconds(first.timestamp_, second.timestamp_);
}

double PMT::joules(const State &first, const State &second) {
  return second.joules_[0] - first.joules_[0];
}

double PMT::watts(const State &first, const State &second) {
  return joules(first, second) / seconds(first, second);
}

unsigned int PMT::GetDumpInterval() {
  const char *dump_interval_ = std::getenv(kDumpIntervalVariable.c_str());
  if (dump_interval_) {
    unsigned int dump_interval = std::stoi(dump_interval_);
    assert(dump_interval > 0);
    return dump_interval;
  } else {
    return GetMeasurementInterval();
  }
}

std::string State::name(int i) {
  assert(i < nr_measurements_);
  return name_[i];
}

float State::joules(int i) {
  assert(i < nr_measurements_);
  return joules_[i];
}

float State::watts(int i) {
  assert(i < nr_measurements_);
  return watt_[i];
}

void PMT::StartDump(const char *filename) {
  const char *filename_ = std::getenv(kDumpFilenameVariable.c_str());
  if (filename_) {
    filename = filename_;
  }
  if (!filename) {
    filename = GetDumpFilename();
  }
  assert(filename);

  dump_file_ = std::make_unique<std::ofstream>(filename);
  Read();
}

void PMT::StopDump() { dump_file_.reset(); }

void PMT::DumpHeader(const State &state) {
  if (dump_file_ != nullptr) {
    std::unique_lock<std::mutex> lock(dump_file_mutex_);
    *dump_file_ << "timestamp";
    for (const std::string &name : state.name_) {
      *dump_file_ << " " << name;
    }
    *dump_file_ << std::endl;
  }
}

void PMT::Dump(const State &state) {
  if (dump_file_ != nullptr) {
    std::unique_lock<std::mutex> lock(dump_file_mutex_);
    *dump_file_ << std::fixed << std::setprecision(3)
                << seconds(state.timestamp_);
    for (float watt : state.watt_) {
      *dump_file_ << " " << watt;
    }
    *dump_file_ << std::endl;
  }
}

void PMT::Mark(const State &start, const State &current,
               const std::string &message) const {
  if (dump_file_ != nullptr) {
    std::unique_lock<std::mutex> lock(dump_file_mutex_);
    *dump_file_ << "M " << seconds(start, current) << " \"" << message << "\""
                << std::endl;
  }
}

void PMT::SetMeasurementInterval(unsigned int milliseconds) {
  if (milliseconds > 0) {
    measurement_interval_ = milliseconds;
  } else {
    unsigned int measurement_interval = 1;
    State state_first;
    State state_second;
    do {
      state_first = GetState();
      std::this_thread::sleep_for(
          std::chrono::milliseconds(measurement_interval));
      state_second = GetState();
      if (!isEqual(state_first.watt_[0], state_second.watt_[0])) {
        break;
      } else {
        measurement_interval *= 2;
      }
    } while (measurement_interval < 100);
    measurement_interval_ = measurement_interval < 10
                                ? measurement_interval
                                : round(measurement_interval / 10) * 10;
  }
}

Timestamp PMT::GetTime() { return std::chrono::system_clock::now(); }

State PMT::Read() {
  // Get the latest measurement
  state_latest_ = GetState();
  return state_latest_;
}

std::unique_ptr<PMT> Create(const std::string &name,
                            const std::string &argument) {
  if (argument.empty()) {
    // Create PMT instance without argument
#if defined(PMT_BUILD_CRAY)
    if (name == cray::Cray::name) {
      return cray::Cray::Create();
    }
#endif
#if defined(PMT_BUILD_LIKWID)
    if (name == likwid::Likwid::name) {
      return likwid::Likwid::Create();
    }
#endif
#if defined(PMT_BUILD_RAPL)
    if (name == rapl::Rapl::name) {
      return rapl::Rapl::Create();
    }
#endif
#if defined(PMT_BUILD_TEGRA)
    if (name == tegra::Tegra::name) {
      return tegra::Tegra::Create();
    }
#endif
#if defined(PMT_BUILD_NVML)
    if (name == nvml::NVML::name) {
      return nvml::NVML::Create();
    }
#endif
#if defined(PMT_BUILD_NVIDIA)
    if (name == nvidia::NVIDIA::name) {
      return nvidia::NVIDIA::Create();
    }
#endif
#if defined(PMT_BUILD_AMDSMI)
    if (name == amdsmi::AMDSMI::name) {
      return amdsmi::AMDSMI::Create();
    }
#endif
#if defined(PMT_BUILD_ROCM)
    if (name == rocm::ROCM::name) {
      return rocm::ROCM::Create();
    }
#endif
  } else if (isNumber(argument)) {
    // Create PMT instance with device number
    const int device_number = std::stoi(argument);

#if defined(PMT_BUILD_NVML)
    if (name == nvml::NVML::name) {
      return nvml::NVML::Create(device_number);
    }
#endif
#if defined(PMT_BUILD_NVIDIA)
    if (name == nvidia::NVIDIA::name) {
      return nvidia::NVIDIA::Create(device_number);
    }
#endif
#if defined(PMT_BUILD_AMDSMI)
    if (name == amdsmi::AMDSMI::name) {
      return amdsmi::AMDSMI::Create(device_number);
    }
#endif
#if defined(PMT_BUILD_ROCM)
    if (name == rocm::ROCM::name) {
      return rocm::ROCM::Create(device_number);
    }
#endif
  } else {
    // Create PMT instance with string argument

    if (argument.size() > 1) {
#if defined(PMT_BUILD_POWERSENSOR2)
      if (name == powersensor2::PowerSensor2::name) {
        return powersensor2::PowerSensor2::Create(argument.c_str());
      }
#endif
#if defined(PMT_BUILD_POWERSENSOR3)
      if (name == powersensor3::PowerSensor3::name) {
        return powersensor3::PowerSensor3::Create(argument.c_str());
      }
#endif
#if defined(PMT_BUILD_XILINX)
      if (name == xilinx::Xilinx::name) {
        return xilinx::Xilinx::Create(argument.c_str());
      }
#endif
    }
  }

#if defined(DEBUG)
  std::stringstream error;
  error << "Invalid or unavailable platform specified: " << name << std::endl;
  throw std::runtime_error(error.str());
#endif

  return Dummy::Create();
}
}  // end namespace pmt

std::ostream &operator<<(std::ostream &os, const pmt::Timestamp &timestamp) {
  os << pmt::PMT::seconds(timestamp);
  return os;
}
