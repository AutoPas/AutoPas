/**
 * @file MPIParallelizedStrategy.h
 * @author W. Thieme
 * @date 27.05.2020
*/

#pragma once

#include "TuningStrategyInterface.h"
#include "autopas/utils/WrapMPI.h"
#include "autopas/options/TuningStrategyOption.h"

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
    _tuningStrategy->addEvidence(time, iteration);
  }

  const Configuration &getCurrentConfiguration() const override {
    return _tuningStrategy->getCurrentConfiguration();
  }

  bool tune(bool currentInvalid) override {
   return  _tuningStrategy->tune(currentInvalid);
  }

  void reset(size_t iteration) override {
    _tuningStrategy->reset(iteration);
  }

  std::set<ContainerOption> getAllowedContainerOptions() const override {
    return _tuningStrategy->getAllowedContainerOptions();
  }

  void removeN3Option(Newton3Option badN3Option) override {
    _tuningStrategy->removeN3Option(badN3Option);
  }

  bool searchSpaceIsTrivial() const override {
    return _tuningStrategy->searchSpaceIsTrivial();
  }

  bool searchSpaceIsEmpty () const override {
    return _tuningStrategy->searchSpaceIsEmpty();
  }

  /**
   * Getter for the internal tuningStrategy
   * @return the tuning strategy which was provided in the constructor
   */
  inline const TuningStrategyInterface &getTuningStrategy() const {
    return *_tuningStrategy;
  }

private:

  /**
   * The tuning strategy tuning locally
   */
  std::unique_ptr<TuningStrategyInterface> _tuningStrategy;

  /**
   * The mpi communicator holding all ranks that participate in this tuning strategy
   */
  AutoPas_MPI_Comm _comm;

};

}