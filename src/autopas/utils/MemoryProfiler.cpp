/**
 * @file MemoryProfiler.cpp
 * @author F. Gratl
 * @date 9/18/18
 */

#include "MemoryProfiler.h"

size_t autopas::memoryProfiler::currentMemoryUsage() {
#ifdef __linux__
  std::ifstream statusFile(statusFileName);
  if (statusFile.rdstate() != std::ifstream::goodbit) {
    // this seems non-critical so we don't throw an exception
    AutoPasLog(ERROR, "Error opening {}", statusFileName);
    return 0;
  }

  std::string line;
  while (getline(statusFile, line)) {
    auto needleStart = line.find("VmRSS");
    if (needleStart != std::string::npos) {
      auto needleEnd = line.find_first_of(" \n", needleStart);
      // parse whole text after colon. strtol ignores leading whitespaces and stops after number.
      return (size_t)std::strtol(line.substr(needleStart + 7, needleEnd).c_str(), nullptr, 10);
    }
  }

  AutoPasLog(ERROR, "Could not read memory usage!");
  return 0;
#elif defined __APPLE__
  mach_task_basic_info info{};
  mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info), &infoCount) !=
      KERN_SUCCESS) {
    AutoPasLog(ERROR, "Error querying the resident memory size", statusFileName);
    return 0;
  }
  // info.resident_size is the resident memory size in bytes, hence we divide by 1000 to return kilobytes
  return static_cast<size_t>(info.resident_size / 1000);
#endif
}
