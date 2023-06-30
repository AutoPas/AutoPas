/**
 * @file TuningStrategyLoggerWrapper.cpp
 * @author humig
 * @date 24.09.2021
 */

#include "TuningStrategyLoggerWrapper.h"

#include "autopas/tuning/tuningStrategy/TuningLogEntry.h"
#include "autopas/utils/Timer.h"

namespace autopas {

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

void TuningStrategyLoggerWrapper::addEvidence(const Configuration &configuration, const Evidence &evidence) {
  _logOut << tuningLogEntry::writeEvidence(evidence.value, evidence.iteration, configuration) << std::endl;

  _actualTuningStrategy->addEvidence(configuration, evidence);
}

void TuningStrategyLoggerWrapper::optimizeSuggestions(std::vector<Configuration> &configQueue,
                                                      const EvidenceCollection &evidenceCollection) {
  _logOut << tuningLogEntry::writeTune() << std::endl;

  return _actualTuningStrategy->optimizeSuggestions(configQueue, evidenceCollection);
}

void TuningStrategyLoggerWrapper::reset(size_t iteration, size_t tuningPhase, std::vector<Configuration> &configQueue,
                                        const autopas::EvidenceCollection &evidenceCollection) {
  _logOut << tuningLogEntry::writeReset(iteration) << std::endl;

  _actualTuningStrategy->reset(iteration, tuningPhase, configQueue, evidenceCollection);
}

bool TuningStrategyLoggerWrapper::needsLiveInfo() const { return true; }

void TuningStrategyLoggerWrapper::receiveLiveInfo(const LiveInfo &info) {
  _logOut << tuningLogEntry::writeLiveInfo(info) << std::endl;

  _actualTuningStrategy->receiveLiveInfo(info);
}

bool TuningStrategyLoggerWrapper::smoothedHomogeneityAndMaxDensityNeeded() const {
  return _actualTuningStrategy->smoothedHomogeneityAndMaxDensityNeeded();
}
}  // namespace autopas
