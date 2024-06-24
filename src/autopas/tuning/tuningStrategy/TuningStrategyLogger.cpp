/**
 * @file TuningStrategyLogger.cpp
 * @author humig
 * @date 24.09.2021
 */

#include "TuningStrategyLogger.h"

#include "autopas/tuning/tuningStrategy/TuningLogEntry.h"
#include "autopas/utils/Timer.h"

namespace autopas {

TuningStrategyLogger::TuningStrategyLogger(const std::string &outputSuffix) {
  std::stringstream filename;
  filename << "tuningLog_";
  filename << outputSuffix;
  filename << (outputSuffix.empty() or outputSuffix.back() == '_' ? "" : "_");
  filename << utils::Timer::getDateStamp();
  filename << ".txt";
  _logOut.open(filename.str());
}

TuningStrategyLogger::~TuningStrategyLogger() { _logOut.flush(); }

void TuningStrategyLogger::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  _logOut << tuningLogEntry::writeEvidence(evidence.value, evidence.iteration, configuration) << std::endl;
}

void TuningStrategyLogger::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                               const EvidenceCollection &evidenceCollection) {
  _logOut << tuningLogEntry::writeTune() << std::endl;
}

void TuningStrategyLogger::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                                 const autopas::EvidenceCollection &evidenceCollection) {
  _logOut << tuningLogEntry::writeReset(iteration) << std::endl;
}

bool TuningStrategyLogger::needsLiveInfo() const {
  // we want to log this therefore request live info
  return true;
}

void TuningStrategyLogger::receiveLiveInfo(const LiveInfo &info) {
  _logOut << tuningLogEntry::writeLiveInfo(info) << std::endl;
}

TuningStrategyOption TuningStrategyLogger::getOptionType() const { return TuningStrategyOption::tuningStrategyLogger; }
}  // namespace autopas
