/**
 * @file MPIParallelizedStrategy.h
 * @author W. Thieme
 * @date 27.05.2020
 */

#pragma once

#include "TuningStrategyInterface.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/utils/WrapMPI.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"

namespace autopas {

/**
 * Wrapper for other tuning strategies which handles MPI communication between processes.
 * The splitting of the search space is not currently handled by this class, but by AutoPasConfigurationCommunicator.
 * This mainly overwrites the tune() method to globally compare the best configuration.
 */
class MPIParallelizedStrategy : public TuningStrategyInterface {
 public:
  /**
   * Constructor for the wrapper. Assumes that the tuningStrategy has already been constructed with the appropriate
   * search space.
   * @param tuningStrategy The underlying tuning strategy which tunes locally on it's node.
   * @param comm The communicator holding all ranks which participate in this tuning strategy
   */
  MPIParallelizedStrategy(std::unique_ptr<TuningStrategyInterface> tuningStrategy, const AutoPas_MPI_Comm comm)
      : _tuningStrategy(std::move(tuningStrategy)),
        _comm(comm) {}

  void addEvidence(long time, size_t iteration) override {
    if (time < _localOptimalTime || _localOptimalTime == -1) {
      _localOptimalTime = time;
    }
    _tuningStrategy->addEvidence(time, iteration);
  }

  const Configuration &getCurrentConfiguration() const override {
    // cellSizeFactor == -1 iff a config is invalid
    if (_optimalConfiguration.cellSizeFactor == -1) {
      return _tuningStrategy->getCurrentConfiguration();
    } else {
      return _optimalConfiguration;
    }
  }

  bool tune(bool currentInvalid) override;

  void reset(size_t iteration) override {
    _tuningStrategy->reset(iteration);
    _optimalConfiguration = Configuration();
    _allGlobalConfigurationsTested = false;
    _allLocalConfigurationsTested = false;
    _localOptimalTime = -1;
  }

  std::set<ContainerOption> getAllowedContainerOptions() const override {
    return _tuningStrategy->getAllowedContainerOptions();
  }

  void removeN3Option(Newton3Option badN3Option) override { _tuningStrategy->removeN3Option(badN3Option); }

  bool searchSpaceIsTrivial() const override { return _tuningStrategy->searchSpaceIsTrivial(); }

  bool searchSpaceIsEmpty() const override { return _tuningStrategy->searchSpaceIsEmpty(); }

  /**
   * Getter for the internal tuningStrategy
   * @return the tuning strategy which was provided in the constructor
   */
  inline const TuningStrategyInterface &getTuningStrategy() const { return *_tuningStrategy; }

 private:
  /**
   * The tuning strategy tuning locally
   */
  std::unique_ptr<TuningStrategyInterface> _tuningStrategy;

  AutoPas_MPI_Comm _comm;
  AutoPas_MPI_Request _request{AUTOPAS_MPI_REQUEST_NULL};

  /**
   * The globally optimal configuration.
   * Usually holds a value that is not in a given rank's search space
   */
   Configuration _optimalConfiguration{Configuration()};

   bool _allLocalConfigurationsTested{false};
   bool _allGlobalConfigurationsTested{false};

   long _localOptimalTime{-1};
};

bool MPIParallelizedStrategy::tune(bool currentInvalid) {
  if (not _allLocalConfigurationsTested) {
    _allLocalConfigurationsTested = not _tuningStrategy->tune(currentInvalid);
  }

  // Wait for the Iallreduce from the last tuning step to finish
  // Make all ranks ready for global tuning simultaneously
  AutoPas_MPI_Wait(&_request, AUTOPAS_MPI_STATUS_IGNORE);
  if (_allGlobalConfigurationsTested) {
    _optimalConfiguration = AutoPasConfigurationCommunicator::optimizeConfiguration(
        _comm, _tuningStrategy->getCurrentConfiguration(), _localOptimalTime);

    return false;
  }

  AutoPas_MPI_Iallreduce(&_allLocalConfigurationsTested, &_allGlobalConfigurationsTested, 1, AUTOPAS_MPI_CXX_BOOL,
                         AUTOPAS_MPI_LAND, _comm, &_request);

  return true;
}


}  // namespace autopas
