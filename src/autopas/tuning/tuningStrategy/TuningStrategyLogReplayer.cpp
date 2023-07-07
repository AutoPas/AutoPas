/**
 * @file TuningStrategyLogReplayer.cpp
 * @author humig
 * @date 24.09.2021
 */

#include "TuningStrategyLogReplayer.h"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <set>
#include <vector>

#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/searchSpace/Evidence.h"
#include "autopas/tuning/searchSpace/EvidenceCollection.h"
#include "autopas/tuning/tuningStrategy/TuningLogEntry.h"

namespace autopas {
TuningStrategyLogReplayer::TuningStrategyLogReplayer(std::string filename,
                                                     std::shared_ptr<TuningStrategyInterface> tuningStrategy,
                                                     const std::set<Configuration> &searchSpace)
    : _filename(std::move(filename)), _tuningStrategy(std::move(tuningStrategy)), _searchSpace(searchSpace) {}

std::optional<Configuration> TuningStrategyLogReplayer::replay() {
  std::ifstream in{_filename};

  if (not in.is_open()) {
    AutoPasLog(ERROR, "Could not open file {}", _filename);
    AutoPasLog(ERROR, "Exiting!");
    exit(-1);
  }

  std::unordered_map<Configuration, std::pair<long, size_t>, ConfigHash> traversalTimes;
  std::vector<Configuration> configQueue{};
  configQueue.reserve(_searchSpace.size());
  std::copy(_searchSpace.rbegin(), _searchSpace.rend(), std::back_inserter(configQueue));
  EvidenceCollection evidenceCollection;
  size_t tuningPhase{0};

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
      while (not configQueue.empty()) {
        _tuningStrategy->optimizeSuggestions(configQueue, evidenceCollection);
        try {
          const auto &[time, iteration] = traversalTimes.at(configQueue.back());
          Evidence evidence{iteration, tuningPhase, time};
          evidenceCollection.addEvidence(configQueue.back(), evidence);
          _tuningStrategy->addEvidence(configQueue.back(), evidence);
        } catch (std::out_of_range &e) {
          // std::cerr << "Could not find data for configuration " << _tuningStrategy->getCurrentConfiguration() <<
          // std::endl;
        }
      }

      if (evidenceAdded) {
        const auto [conf, evidence] = evidenceCollection.getLatestOptimalConfiguration();
        AutoPasLog(INFO, "Best Configuration found: {}", conf.toShortString());
        bestConfiguration = conf;
        evidenceAdded = false;
      }

      const auto &iteration = tuningLogEntry::readReset(stream);
      ++tuningPhase;
      _tuningStrategy->reset(iteration, tuningPhase, configQueue, evidenceCollection);
      traversalTimes.clear();
    }
  }

  return bestConfiguration;
}
}  // namespace autopas