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

class MPIParallelizedStrategy : public TuningStrategyInterface {
public:

  MPIParallelizedStrategy(const std::set<ContainerOption> &allowedContainerOptions,
                          const std::set<double> &allowedCellSizeFactors,
                          const std::set<TraversalOption> &allowedTraversalOptions,
                          const std::set<DataLayoutOption> &allowedDataLayoutOptions,
                          const std::set<Newton3Option> &allowedNewton3Options,
                          TuningStrategyInterface tuningStrategy,
                          AutoPas_MPI_Comm comm)
    : _tuningStrategy(std::make_unique<TuningStrategyInterface>(tuningStrategy)) {

  }

  void addEvidence(long time, size_t iteration) {
    _tuningStrategy->addEvidence(time, iteration);
  }

  const Configuration &getCurrentConfiguration() const {
    _tuningStrategy->getCurrentConfiguration();
  }

  bool tune(bool currentInvalid) {
    _tuningStrategy->tune(currentInvalid);
  }

  void reset(size_t iteration) {
    _tuningStrategy->reset(iteration);
  }

  std::set<ContainerOption> getAllowedContainerOptions() const {
    _tuningStrategy->getAllowedContainerOptions();
  }

  void removeN3Option(Newton3Option badN3Option) {
    _tuningStrategy->removeN3Option(badN3Option);
  }

  bool searchSpaceIsTrivial() const {
    _tuningStrategy->searchSpaceIsTrivial();
  }

  bool searchSpaceIsEmpty () const {
    _tuningStrategy->searchSpaceIsEmpty();
  }

  inline const TuningStrategyInterface &getTuningStrategy() {
    return *_tuningStrategy;
  }

private:

  std::unique_ptr<TuningStrategyInterface> _tuningStrategy;

};

}