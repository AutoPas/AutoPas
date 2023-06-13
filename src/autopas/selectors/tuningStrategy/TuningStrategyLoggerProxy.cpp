/**
 * @file TuningStrategyLoggerProxy.cpp
 * @author humig
 * @date 24.09.2021
 */

#include <chrono>

#include "TuningStrategyLoggerWrapper.h"
#include "autopas/utils/Timer.h"

namespace autopas {

// Contains some helpers to write and read the tuning log entries.
namespace tuningLogEntry {
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
std::string writeEvidence(long time, size_t iteration, const Configuration &config) {
  return toString(std::string{"evidence"}, time, iteration, config);
}

/**
 * Reads  the arguments of an evidence entry in the log file from a stringstream.
 * @param str The stream to read from.
 * @return The evidence as a tuple.
 */
std::tuple<long, size_t, Configuration> readEvidence(std::stringstream &str) {
  return fromString<long, size_t, Configuration>(str);
}

/**
 * Writes a tune entry for the log file into a string.
 * @param currentInvalid If current configuration was invalid.
 * @return The string with the tune entry.
 */
std::string writeTune(bool currentInvalid) { return toString(std::string{"tune"}, currentInvalid); }

/**
 * Reads the arguments of a tune entry from the stringstream.
 * @param str The stringstream.
 * @return The currentInvalid bool.
 */
bool readTune(std::stringstream &str) { return std::get<0>(fromString<bool>(str)); }

/**
 * Writes a reset entry in the log file to a string.
 * @param iteration The iteration it was performed on.
 * @return The string with the reset entry.
 */
std::string writeReset(size_t iteration) { return toString(std::string{"reset"}, iteration); }

/**
 * Reads the arguments of a reset entry in the log file from a string.
 * @param str The stringstream to read from.
 * @return The iteration the reset happened in.
 */
size_t readReset(std::stringstream &str) { return std::get<0>(fromString<size_t>(str)); }

/**
 * Writes a liveInfo entry in the log file to a string.
 * @param liveInfo The live info to write.
 * @return The string with the live info.
 */
std::string writeLiveInfo(const LiveInfo &liveInfo) { return toString(std::string{"liveInfo"}, liveInfo); }

/**
 * Reads the arguments of a live info entry in the log file from a stringstream.
 * @param str The stringstream to read from.
 * @return The LiveInfo read from the stream.
 */
LiveInfo readLiveInfo(std::stringstream &str) { return std::get<0>(fromString<LiveInfo>(str)); }
};  // namespace tuningLogEntry

TuningStrategyLoggerWrapper::TuningStrategyLoggerWrapper(std::unique_ptr<TuningStrategyInterface> actualTuningStrategy,
                                                         const std::string &outputSuffix)
    : _actualTuningStrategy(std::move(actualTuningStrategy)) {
  std::stringstream filename;
  filename << "tuningLog-";
  filename << utils::Timer::getDateStamp();
  filename << '-' << outputSuffix;
  filename << ".txt";
  _logOut.open(filename.str());
}

TuningStrategyLoggerWrapper::~TuningStrategyLoggerWrapper() { _logOut.flush(); }

void TuningStrategyLoggerWrapper::addEvidence(long time, size_t iteration) {
  _actualTuningStrategy->addEvidence(time, iteration);

  _logOut << tuningLogEntry::writeEvidence(time, iteration, _actualTuningStrategy->getCurrentConfiguration())
          << std::endl;
}

long TuningStrategyLoggerWrapper::getEvidence(Configuration configuration) const {
  return _actualTuningStrategy->getEvidence(configuration);
}

const Configuration &TuningStrategyLoggerWrapper::getCurrentConfiguration() const {
  return _actualTuningStrategy->getCurrentConfiguration();
}

bool TuningStrategyLoggerWrapper::tune(bool currentInvalid) {
  _logOut << tuningLogEntry::writeTune(currentInvalid) << std::endl;

  return _actualTuningStrategy->tune(currentInvalid);
}

void TuningStrategyLoggerWrapper::reset(size_t iteration) {
  _logOut << tuningLogEntry::writeReset(iteration) << std::endl;

  _actualTuningStrategy->reset(iteration);
}

bool TuningStrategyLoggerWrapper::needsLiveInfo() const { return true; }

void TuningStrategyLoggerWrapper::receiveLiveInfo(const LiveInfo &info) {
  _logOut << tuningLogEntry::writeLiveInfo(info) << std::endl;

  _actualTuningStrategy->receiveLiveInfo(info);
}

std::set<ContainerOption> TuningStrategyLoggerWrapper::getAllowedContainerOptions() const {
  return _actualTuningStrategy->getAllowedContainerOptions();
}

void TuningStrategyLoggerWrapper::removeN3Option(Newton3Option option) {
  _actualTuningStrategy->removeN3Option(option);
}
bool TuningStrategyLoggerWrapper::searchSpaceIsTrivial() const { return _actualTuningStrategy->searchSpaceIsTrivial(); }

bool TuningStrategyLoggerWrapper::searchSpaceIsEmpty() const { return _actualTuningStrategy->searchSpaceIsEmpty(); }

bool TuningStrategyLoggerWrapper::smoothedHomogeneityAndMaxDensityNeeded() const {
  return _actualTuningStrategy->smoothedHomogeneityAndMaxDensityNeeded();
}

TuningStrategyLogReplayer::TuningStrategyLogReplayer(std::string filename,
                                                     std::shared_ptr<TuningStrategyInterface> tuningStrategy)
    : _filename(std::move(filename)), _tuningStrategy(std::move(tuningStrategy)) {}

std::optional<Configuration> TuningStrategyLogReplayer::replay() {
  std::ifstream in{_filename};

  if (not in.is_open()) {
    AutoPasLog(ERROR, "Could not open file {}", _filename);
    AutoPasLog(ERROR, "Exiting!");
    exit(-1);
  }

  std::unordered_map<Configuration, std::pair<long, size_t>, ConfigHash> traversalTimes;

  std::optional<Configuration> bestConfiguration;
  bool evidenceAdded = false;
  while (not in.eof()) {
    std::string line;
    std::getline(in, line, '\n');

    std::stringstream stream{line};
    std::string type;
    std::getline(stream, type, ' ');

    if (type == "evidence") {
      evidenceAdded = true;
      const auto &[time, iteration, config] = tuningLogEntry::readEvidence(stream);
      // std::cout << "Add data for configuration " << config << std::endl;
      traversalTimes[config] = {time, iteration};
    } else if (type == "tune") {
      // Do nothing in former tune
    } else if (type == "liveInfo") {
      const auto &liveInfo = tuningLogEntry::readLiveInfo(stream);
      AutoPasLog(INFO, "\t{}", liveInfo.toString());
      _tuningStrategy->receiveLiveInfo(liveInfo);
    } else if (type == "reset" || in.eof()) {
      auto cont = not traversalTimes.empty();
      while (cont) {
        cont = _tuningStrategy->tune();
        try {
          const auto &[time, iteration] = traversalTimes.at(_tuningStrategy->getCurrentConfiguration());
          _tuningStrategy->addEvidence(time, iteration);
        } catch (std::out_of_range &e) {
          // std::cerr << "Could not find data for configuration " << _tuningStrategy->getCurrentConfiguration() <<
          // std::endl;
        }
      }

      if (evidenceAdded) {
        AutoPasLog(INFO, "Best Configuration found: {}", _tuningStrategy->getCurrentConfiguration().toShortString());
        bestConfiguration = _tuningStrategy->getCurrentConfiguration();
        evidenceAdded = false;
      }

      const auto &iteration = tuningLogEntry::readReset(stream);
      _tuningStrategy->reset(iteration);
      traversalTimes.clear();
    }
  }

  return bestConfiguration;
}

}  // namespace autopas
