/**
 * @file TuningStrategyLoggerWrapper.cpp
 * @author humig
 * @date 24.09.2021
 */

#include "TuningStrategyLoggerWrapper.h"

#include "autopas/selectors/tuningStrategy/TuningLogEntry.h"
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
}  // namespace autopas
