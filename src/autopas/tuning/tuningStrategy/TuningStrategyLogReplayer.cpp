/**
 * @file TuningStrategyLogReplayer.cpp
 * @author humig
 * @date 24.09.2021
 */

#include "TuningStrategyLogReplayer.h"

#include <fstream>

#include "autopas/tuning/tuningStrategy/TuningLogEntry.h"

namespace autopas {
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