/**
 * @file TuningStrategyLoggerProxy.cpp
 * @author humig
 * @date 24.09.2021
 */

#include "TuningStrategyLoggerProxy.h"

#include <chrono>

#include "autopas/utils/Timer.h"

namespace autopas {

namespace tuningLogEntry {
template <class T>
void toStringHelper(std::ostream &in, const T &val) {
  in << val << ' ';
}

template <class... Payload>
auto toString(const Payload &... payload) {
  std::stringstream stream;
  (toStringHelper(stream, payload), ...);
  return stream.str();
}

template <class... Payload>
std::tuple<Payload...> fromString(std::stringstream &stream) {
  std::tuple<Payload...> tuple{};
  ((stream >> std::get<Payload>(tuple)), ...);
  return tuple;
}

std::string writeEvidence(long time, size_t iteration, const Configuration &config) {
  return toString(std::string{"evidence"}, time, iteration, config);
}

std::tuple<long, size_t, Configuration> readEvidence(std::stringstream &str) {
  return fromString<long, size_t, Configuration>(str);
}

std::string writeTune(bool currentInvalid) { return toString(std::string{"tune"}, currentInvalid); }

bool readTune(std::stringstream &str) { return std::get<0>(fromString<bool>(str)); }

std::string writeReset(size_t iteration) { return toString(std::string{"reset"}, iteration); }

size_t readReset(std::stringstream &str) { return std::get<0>(fromString<size_t>(str)); }

std::string writeLiveInfo(const LiveInfo &liveInfo) { return toString(std::string{"liveInfo"}, liveInfo); }

LiveInfo readLiveInfo(std::stringstream &str) { return std::get<0>(fromString<LiveInfo>(str)); }
};  // namespace tuningLogEntry

TuningStrategyLoggerProxy::TuningStrategyLoggerProxy(std::unique_ptr<TuningStrategyInterface> actualTuningStrategy,
                                                     const std::string &outputSuffix)
    : _actualTuningStrategy(std::move(actualTuningStrategy)) {
  std::stringstream filename;
  filename << "tuningLog-";
  filename << utils::Timer::getDateStamp();
  filename << '-' << outputSuffix;
  filename << ".txt";
  _logOut.open(filename.str());
}

TuningStrategyLoggerProxy::~TuningStrategyLoggerProxy() { _logOut.flush(); }

void TuningStrategyLoggerProxy::addEvidence(long time, size_t iteration) {
  _actualTuningStrategy->addEvidence(time, iteration);

  _logOut << tuningLogEntry::writeEvidence(time, iteration, _actualTuningStrategy->getCurrentConfiguration())
          << std::endl;
}

long TuningStrategyLoggerProxy::getEvidence(Configuration configuration) const {
  return _actualTuningStrategy->getEvidence(configuration);
}

const Configuration &TuningStrategyLoggerProxy::getCurrentConfiguration() const {
  return _actualTuningStrategy->getCurrentConfiguration();
}

bool TuningStrategyLoggerProxy::tune(bool currentInvalid) {
  _logOut << tuningLogEntry::writeTune(currentInvalid) << std::endl;

  return _actualTuningStrategy->tune(currentInvalid);
}

void TuningStrategyLoggerProxy::reset(size_t iteration) {
  _logOut << tuningLogEntry::writeReset(iteration) << std::endl;

  _actualTuningStrategy->reset(iteration);
}

bool TuningStrategyLoggerProxy::needsLiveInfo() const { return _actualTuningStrategy->needsLiveInfo(); }

void TuningStrategyLoggerProxy::receiveLiveInfo(LiveInfo info) {
  _logOut << tuningLogEntry::writeLiveInfo(info) << std::endl;

  _actualTuningStrategy->receiveLiveInfo(info);
}

std::set<ContainerOption> TuningStrategyLoggerProxy::getAllowedContainerOptions() const {
  return _actualTuningStrategy->getAllowedContainerOptions();
}

void TuningStrategyLoggerProxy::removeN3Option(Newton3Option option) { _actualTuningStrategy->removeN3Option(option); }
bool TuningStrategyLoggerProxy::searchSpaceIsTrivial() const { return _actualTuningStrategy->searchSpaceIsTrivial(); }

bool TuningStrategyLoggerProxy::searchSpaceIsEmpty() const { return _actualTuningStrategy->searchSpaceIsEmpty(); }

TuningStrategyLogReplayer::TuningStrategyLogReplayer(std::string filename,
                                                     std::unique_ptr<TuningStrategyInterface> tuningStrategy)
    : _filename(std::move(filename)), _tuningStrategy(std::move(tuningStrategy)) {}

void TuningStrategyLogReplayer::replay() const {
  std::ifstream in{_filename};

  if (not in.is_open()) {
    std::cerr << "Could not open file " << _filename << std::endl;
    std::cerr << "Exiting!" << std::endl;
    exit(-1);
  }

  std::unordered_map<Configuration, std::pair<long, size_t>, ConfigHash> traversalTimes;

  while (not in.eof()) {
    std::string line;
    std::getline(in, line, '\n');

    std::stringstream stream{line};
    std::string type;
    std::getline(stream, type, ' ');

    if (type == "evidence") {
      const auto &[time, iteration, config] = tuningLogEntry::readEvidence(stream);
      traversalTimes[config] = {time, iteration};
    } else if (type == "tune") {
      // Do nothing in former tune
    } else if (type == "liveInfo") {
      const auto &liveInfo = tuningLogEntry::readLiveInfo(stream);
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

      // std::cout << "Best Configuration found: " << _tuningStrategy->getCurrentConfiguration() << std::endl;

      const auto &iteration = tuningLogEntry::readReset(stream);
      _tuningStrategy->reset(iteration);
      traversalTimes.clear();
    }
  }
}

}  // namespace autopas
