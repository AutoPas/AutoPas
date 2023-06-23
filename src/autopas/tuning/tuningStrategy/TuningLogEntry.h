/**
 * @file TuningLogEntry.h
 * @author humig
 * @date 24.09.2021
 */

#pragma once

#include <ostream>
#include <sstream>
#include <tuple>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/LiveInfo.h"

namespace autopas::tuningLogEntry {
/**
 * Writes the given argument into the given ostream with following white space.
 * @tparam T type of the argument.
 * @param in The ostream.
 * @param val The argument to write.
 */
template <class T>
void toStringHelper(std::ostream &in, const T &val) {
  in << val << ' ';
}

/**
 * Writes multiple arguments into a string using their << operator.
 * @tparam Payload The types of the arguments.
 * @param payload The arguments to write.
 * @return The string representation of the arguments, with whitespaces between.
 */
template <class... Payload>
auto toString(const Payload &... payload) {
  std::stringstream stream;
  (toStringHelper(stream, payload), ...);
  return stream.str();
}

/**
 * Reads multiple values from a stringstream using their >> operator.
 * @tparam Payload The types of the valuaes to read.
 * @param stream The stringstream to read from.
 * @return The read values as a tuple.
 */
template <class... Payload>
std::tuple<Payload...> fromString(std::stringstream &stream) {
  std::tuple<Payload...> tuple{};
  ((stream >> std::get<Payload>(tuple)), ...);
  return tuple;
}

/**
 * Writes evidence to a string.
 * @param time The measured time.
 * @param iteration The iteration in was measured in.
 * @param config The configuation used.
 * @return The string with the evidence.
 */
std::string writeEvidence(long time, size_t iteration, const Configuration &config);

/**
 * Reads  the arguments of an evidence entry in the log file from a stringstream.
 * @param str The stream to read from.
 * @return The evidence as a tuple.
 */
std::tuple<long, size_t, Configuration> readEvidence(std::stringstream &str);

/**
 * Writes a tune entry for the log file into a string.
 * @param currentInvalid If current configuration was invalid.
 * @return The string with the tune entry.
 */
std::string writeTune(bool currentInvalid);

/**
 * Reads the arguments of a tune entry from the stringstream.
 * @param str The stringstream.
 * @return The currentInvalid bool.
 */
bool readTune(std::stringstream &str);

/**
 * Writes a reset entry in the log file to a string.
 * @param iteration The iteration it was performed on.
 * @return The string with the reset entry.
 */
std::string writeReset(size_t iteration);

/**
 * Reads the arguments of a reset entry in the log file from a string.
 * @param str The stringstream to read from.
 * @return The iteration the reset happened in.
 */
size_t readReset(std::stringstream &str);

/**
 * Writes a liveInfo entry in the log file to a string.
 * @param liveInfo The live info to write.
 * @return The string with the live info.
 */
std::string writeLiveInfo(const LiveInfo &liveInfo);

/**
 * Reads the arguments of a live info entry in the log file from a stringstream.
 * @param str The stringstream to read from.
 * @return The LiveInfo read from the stream.
 */
LiveInfo readLiveInfo(std::stringstream &str);
}  // namespace autopas::tuningLogEntry