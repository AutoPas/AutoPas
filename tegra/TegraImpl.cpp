#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

#include <ext/alloc_traits.h>
#include <signal.h>

#include "TegraImpl.h"

namespace detail {
bool FileExists(const std::string& name) {
  if (FILE* file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}

std::string Execute(const std::string& commandline) {
  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(commandline.c_str(), "r"),
                                                pclose);
  if (!pipe) {
    throw std::runtime_error("popen() failed!");
  }
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
    result += buffer.data();
  }
  return result;
}

std::vector<std::string> SplitString(const std::string& string,
                                     std::string delimiter = " ") {
  std::vector<std::string> substrings;
  int start = 0;
  int end = string.find_first_of(" \n");
  while (end != -1) {
    start = end + delimiter.size();
    end = string.find(delimiter, start);
    substrings.push_back(string.substr(start, end - start));
  }
  substrings.push_back(string.substr(start, end - start));
  return substrings;
}

std::string FindLogfile() {
  const std::string ps = Execute("ps aux | grep tegrastats");
  std::vector<std::string> substrings = SplitString(ps);
  bool have_tegrastats = false;
  int index_logfile = -1;
  for (size_t i = 0; i < substrings.size(); ++i) {
    if (substrings[i].compare("tegrastats") == 0 ||
        substrings[i].compare("/usr/bin/tegrastats") == 0) {
      have_tegrastats = true;
    }
    if (substrings[i].compare("--logfile") == 0) {
      index_logfile = i + 1;
      break;
    }
  }
  if (have_tegrastats && index_logfile > 0) {
    return substrings[index_logfile];
  } else {
    return "";
  }
}

std::string StartTegraStats(int interval) {
  char filename[32] = "/tmp/tegrastats-XXXXXX";
  if (mkstemp(filename) == -1) {
    throw std::runtime_error("Could not create temporary file");
  }
  const char* binary = "/usr/bin/tegrastats";
  std::stringstream commandline;
  commandline << binary << " --start --interval " << interval << " --logfile "
              << filename;
  std::string output = Execute(commandline.str());
  if (output.compare("") != 0) {
    throw std::runtime_error("Error starting tegrastats: " + output);
  }
  return filename;
}

void StopTegraStats(const std::string& logfile) {
  if (logfile.compare("") != 0) {
    const char* binary = "/usr/bin/tegrastats";
    std::stringstream commandline;
    commandline << binary << " --stop";
    Execute(commandline.str());
    std::remove(logfile.c_str());
  }
}

std::string ReadLastLine(const std::string& file_name) {
  std::ifstream file(file_name);
  if (!file.is_open()) {
    return "";
  }

  std::string last_line;
  std::string current_line;

  while (std::getline(file, current_line)) {
    if (!current_line.empty()) {
      last_line = current_line;
    }
  }

  file.close();
  return last_line;
}

std::vector<pmt::tegra::TegraMeasurement> ReadPowerMeasurements(
    const std::string& line) {
  std::vector<pmt::tegra::TegraMeasurement> measurements;

  std::vector<std::regex> regexes;
  regexes.emplace_back(
      "(\\w+) (\\d+)mW\\/(\\d+)mW");  // Jetson AGX Xavier and Jetson Orin Nano
  regexes.emplace_back("(POM_\\w+) (\\d+)\\/(\\d+)");  // Jetson Nano

  for (std::regex& regex : regexes) {
    std::smatch matches;
    std::string input = line;

    while (std::regex_search(input, matches, regex)) {
      const std::string name = matches[1].str();
      const int instantaneous_mw = atoi(matches[2].str().c_str());
      const int average_mw = atoi(matches[3].str().c_str());
      measurements.emplace_back(name, instantaneous_mw);
      input = matches.suffix();
    }
  }

  return measurements;
}

bool CheckSensors(const std::vector<std::string>& sensors,
                  std::vector<pmt::tegra::TegraMeasurement>& measurements) {
  if (sensors.size() != measurements.size()) {
    return false;
  }
  const size_t n = sensors.size();
  for (size_t i = 0; i < n; i++) {
    if (sensors[i].compare(measurements[i].first) != 0) {
      return false;
    }
  }
  return true;
}
}  // end namespace detail

namespace pmt::tegra {

void SignalCallbackHandler(int num) {
  const std::string logfile = detail::FindLogfile();
  detail::StopTegraStats(logfile);
  signal(SIGINT, SIG_DFL);
  raise(num);
}

TegraImpl::TegraImpl() {
  filename_ = detail::FindLogfile();
  if (!detail::FileExists(filename_)) {
    detail::StopTegraStats(filename_);
    filename_ = detail::StartTegraStats(measurement_interval_);
    started_tegrastats_ = true;
    signal(SIGINT, SignalCallbackHandler);
  }

  previous_state_ = GetTegraState();
  previous_state_.joules_ = 0;
}

TegraImpl::~TegraImpl() {
  if (started_tegrastats_) {
    detail::StopTegraStats(filename_);
    started_tegrastats_ = false;
  }
}

std::vector<TegraMeasurement> TegraImpl::GetMeasurements() {
  std::string line = detail::ReadLastLine(filename_);

  while (line.compare("") == 0) {
    std::this_thread::sleep_for(
        std::chrono::milliseconds(measurement_interval_));
    line = detail::ReadLastLine(filename_);
  }

  return detail::ReadPowerMeasurements(line);
};

const std::vector<std::string> sensors_agx_xavier{"GPU", "CPU",   "SOC",
                                                  "CV",  "VDDRQ", "SYS5V"};
const std::vector<std::string> sensors_agx_orin{
    "VDD_GPU_SOC", "VDD_CPU_CV", "VIN_SYS_5V0", "VDDQ_VDD2_1V8AO"};
const std::vector<std::string> sensors_jetson_nano{"POM_5V_IN", "POM_5V_GPU",
                                                   "POM_5V_CPU"};
const std::vector<std::string> sensors_jetson_orin_nano{
    "VDD_IN", "VDD_CPU_GPU_CV", "VDD_SOC"};

TegraState::operator State() {
  State state(1 + measurements.size());
  state.timestamp_ = timestamp_;
  state.name_[0] = "total";
  state.joules_[0] = joules_ * 1e-3;
  state.watt_[0] = watt_ * 1e-3;

  for (size_t i = 0; i < measurements.size(); i++) {
    const std::string name = measurements[i].first;
    const double watt = double(measurements[i].second) / 1000;
    state.name_[i + 1] = name;
    state.watt_[i + 1] = watt;
  }

  return state;
}

TegraState TegraImpl::GetTegraState() {
  TegraState state;
  state.timestamp_ = GetTime();
  state.measurements = GetMeasurements();

  // Compute total power consumption as sum of individual measurements.
  // Which individual measurements to use differs per platform.
  state.watt_ = 0;
  if (detail::CheckSensors(sensors_agx_xavier, state.measurements)) {
    // Jetson AGX Xavier: sum all sensors
    for (auto& measurement : state.measurements) {
      state.watt_ += measurement.second;
    }
  } else if (detail::CheckSensors(sensors_agx_orin, state.measurements)) {
    // Jetson AGX Orin: sum all sensors
    for (auto& measurement : state.measurements) {
      state.watt_ += measurement.second;
    }
  } else if (detail::CheckSensors(sensors_jetson_nano, state.measurements)) {
    // Jetson Nano: POM_5V_IN only
    state.watt_ += state.measurements[0].second;
  } else if (detail::CheckSensors(sensors_jetson_orin_nano,
                                  state.measurements)) {
    // Jetson Nano: VDD_IN only
    state.watt_ += state.measurements[0].second;
  }

  const float watt = (state.watt_ + previous_state_.watt_) / 2;
  const double duration = seconds(previous_state_.timestamp_, state.timestamp_);
  state.joules_ = previous_state_.joules_ + (watt * duration);

  previous_state_ = state;

  return state;
}

}  // end namespace pmt::tegra