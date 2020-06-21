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

  inline void addEvidence(long time, size_t iteration) override { _tuningStrategy->addEvidence(time, iteration); }

  [[nodiscard]] inline long getEvidence(Configuration configuration) const override {
    return _tuningStrategy->getEvidence(configuration);
  }

  [[nodiscard]] const Configuration &getCurrentConfiguration() const override {
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
  }

  [[nodiscard]] std::set<ContainerOption> getAllowedContainerOptions() const override {
    return _tuningStrategy->getAllowedContainerOptions();
  }

  inline void removeN3Option(Newton3Option badN3Option) override { _tuningStrategy->removeN3Option(badN3Option); }

  [[nodiscard]] inline bool searchSpaceIsTrivial() const override { return _tuningStrategy->searchSpaceIsTrivial(); }

  [[nodiscard]] inline bool searchSpaceIsEmpty() const override { return _tuningStrategy->searchSpaceIsEmpty(); }

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
};

bool MPIParallelizedStrategy::tune(bool currentInvalid) {
  int rank;
  AutoPas_MPI_Comm_rank(_comm, &rank);

  if (not _allLocalConfigurationsTested) {
    _allLocalConfigurationsTested = not _tuningStrategy->tune(currentInvalid);
    if (currentInvalid) {
      return true;
    }
  }

  // Wait for the Iallreduce from the last tuning step to finish
  // Make all ranks ready for global tuning simultaneously
  AutoPas_MPI_Wait(&_request, AUTOPAS_MPI_STATUS_IGNORE);
  if (_allGlobalConfigurationsTested) {
    auto config = _tuningStrategy->getCurrentConfiguration();
    _optimalConfiguration = utils::AutoPasConfigurationCommunicator::optimizeConfiguration(
        _comm, config, _tuningStrategy->getEvidence(config));

    return false;
  }

  AutoPas_MPI_Iallreduce(&_allLocalConfigurationsTested, &_allGlobalConfigurationsTested, 1, AUTOPAS_MPI_CXX_BOOL,
                         AUTOPAS_MPI_LAND, _comm, &_request);

  return true;
}


}  // namespace autopas
