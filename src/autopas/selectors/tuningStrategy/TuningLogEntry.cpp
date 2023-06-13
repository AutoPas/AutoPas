/**
 * @file TuningLogEntry.cpp
 * @author humig
 * @date 24.09.2021
 */

#include "TuningLogEntry.h"

namespace autopas::tuningLogEntry {

std::string writeEvidence(long time, size_t iteration, const autopas::Configuration &config) {
  return toString(std::string{"evidence"}, time, iteration, config);
}

std::tuple<long, size_t, Configuration> readEvidence(std::stringstream &str) {
  return fromString<long, size_t, Configuration>(str);
}

std::string writeTune(bool currentInvalid) { return toString(std::string{"tune"}, currentInvalid); }

bool readTune(std::stringstream &str) { return std::get<0>(fromString<bool>(str)); }

std::string writeReset(size_t iteration) { return toString(std::string{"reset"}, iteration); }

size_t readReset(std::stringstream &str) { return std::get<0>(fromString<size_t>(str)); }

std::string writeLiveInfo(const autopas::LiveInfo &liveInfo) { return toString(std::string{"liveInfo"}, liveInfo); }

autopas::LiveInfo readLiveInfo(std::stringstream &str) { return std::get<0>(fromString<LiveInfo>(str)); }

}  // namespace autopas::tuningLogEntry