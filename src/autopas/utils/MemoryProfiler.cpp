/**
 * @file MemoryProfiler.cpp
 * @author F. Gratl
 * @date 9/18/18
 */

#include "MemoryProfiler.h"

#include <sys/stat.h>

#include <fstream>

#include "autopas/utils/logging/Logger.h"

size_t autopas::memoryProfiler::currentMemoryUsage() {
  std::ifstream statusFile(statusFileName);
  if (statusFile.rdstate() != std::ifstream::goodbit) {
    // this seems non-critical so we don't throw an exception
    AutoPasLog(error, "Error opening {}", statusFileName);
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

  AutoPasLog(error, "Could not read memory usage!");
  return 0;
}
