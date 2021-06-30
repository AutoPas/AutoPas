/**
 * @file MPIParallelizedStrategy.h
 * @author W. Thieme
 * @date 27.05.2020
 */

#pragma once

#include "TuningStrategyInterface.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/utils/AutoPasConfigurationCommunicator.h"
#include "autopas/utils/ConfigurationAndRankIteratorHandler.h"
#include "autopas/utils/ExceptionHandler.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/WrapMPI.h"

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
   * If the underlying tuning strategy fails for some reason, uses the fallback-options to keep global search going.
   * @param tuningStrategy The underlying tuning strategy which tunes locally on its node.
   * @param comm The communicator holding all ranks which participate in this tuning strategy
   * @param fallbackContainers
   * @param fallbackTraversals
   * @param fallbackLoadEstimators
   * @param fallbackDataLayouts
   * @param fallbackNewton3s
   */
  MPIParallelizedStrategy(std::unique_ptr<TuningStrategyInterface> tuningStrategy, const AutoPas_MPI_Comm comm,
                          std::set<autopas::ContainerOption> fallbackContainers,
                          std::set<autopas::TraversalOption> fallbackTraversals,
                          std::set<autopas::LoadEstimatorOption> fallbackLoadEstimators,
                          std::set<autopas::DataLayoutOption> fallbackDataLayouts,
                          std::set<autopas::Newton3Option> fallbackNewton3s)
      : _tuningStrategy(std::move(tuningStrategy)),
        _comm(comm),
        _fallbackContainers(std::move(fallbackContainers)),
        _fallbackTraversalOptions(std::move(fallbackTraversals)),
        _fallbackLoadEstimators(std::move(fallbackLoadEstimators)),
        _fallbackDataLayouts(std::move(fallbackDataLayouts)),
        _fallbackNewton3s(std::move(fallbackNewton3s)) {}

  inline void addEvidence(long time, size_t iteration) override { _tuningStrategy->addEvidence(time, iteration); }

  [[nodiscard]] inline long getEvidence(Configuration configuration) const override {
    return _tuningStrategy->getEvidence(configuration);
  }

  [[nodiscard]] const Configuration &getCurrentConfiguration() const override {
    if (not _optimalConfiguration.hasValidValues()) {
      return _tuningStrategy->getCurrentConfiguration();
    } else {
      return _optimalConfiguration;
    }
  }

  bool tune(bool currentInvalid) override;

  void reset(size_t iteration) override {
    _optimalConfiguration = Configuration();
    _allGlobalConfigurationsTested = false;
    _allLocalConfigurationsTested = false;
    _strategyStillWorking = true;
    if (_configIterator != nullptr) {
      _configIterator.reset();
    }
    try {
      _tuningStrategy->reset(iteration);
    } catch (utils::ExceptionHandler::AutoPasException &exception) {
      AutoPasLog(warn,
                 "MPIParallelizedStrategy: Underlying strategy failed (with error: {}). Reverting to fallback-mode.",
                 exception.what());
      setupFallbackOptions();
    }
  }

  /**
   * Reset all internal parameters to the beginning of a new tuning cycle.
   * Also distribute ranks in buckets for MPI tuning
   * @tparam Particle
   * @param iteration Gives the current iteration to the tuning strategy.
   * @param container container of current simulation
   * @param MPITuningMaxDifferenceForBucket For MPI-tuning: Maximum of the relative difference in the comparison metric
   * for two ranks which exchange their tuning information.
   * @param MPITuningWeightForMaxDensity For MPI-tuning: Weight for maxDensity in the calculation for bucket
   * distribution.
   */
  template <class Particle>
  void reset(size_t iteration, std::shared_ptr<autopas::ParticleContainerInterface<Particle>> container,
             double MPITuningMaxDifferenceForBucket, double MPITuningWeightForMaxDensity) {
    _optimalConfiguration = Configuration();
    _allGlobalConfigurationsTested = false;
    _allLocalConfigurationsTested = false;
    _strategyStillWorking = true;
    if (_configIterator != nullptr) {
      _configIterator.reset();
    }
    autopas::utils::AutoPasConfigurationCommunicator::distributeRanksInBuckets<Particle>(
        _comm, &_bucket, container, MPITuningMaxDifferenceForBucket, MPITuningWeightForMaxDensity);
    AutoPasLog(debug, "finished bucket distribution");
    try {
      _tuningStrategy->reset(iteration);
    } catch (utils::ExceptionHandler::AutoPasException &exception) {
      AutoPasLog(warn,
                 "MPIParallelizedStrategy: Underlying strategy failed (with error: {}). Reverting to fallback-mode.",
                 exception.what());
      setupFallbackOptions();
    }
  }

  [[nodiscard]] std::set<ContainerOption> getAllowedContainerOptions() const override {
    return _tuningStrategy->getAllowedContainerOptions();
  }

  inline void removeN3Option(Newton3Option badN3Option) override { _tuningStrategy->removeN3Option(badN3Option); }

  [[nodiscard]] inline bool searchSpaceIsTrivial() const override { return _tuningStrategy->searchSpaceIsTrivial(); }

  [[nodiscard]] inline bool searchSpaceIsEmpty() const override { return _tuningStrategy->searchSpaceIsEmpty(); }

 private:
  /**
   * Essentially sets up a full-search-esque configuration selection, in case the underlying search strategy fails.
   * Usually the strategy fails, because it has no valid configurations when many ranks are present.
   */
  void setupFallbackOptions();

  /**
   * Uses _configIterator to get a new config from the global search space.
   */
  void nextFallbackConfig();

  /**
   * The tuning strategy tuning locally
   */
  std::unique_ptr<TuningStrategyInterface> _tuningStrategy;

  AutoPas_MPI_Comm _comm;
  AutoPas_MPI_Comm _bucket{AUTOPAS_MPI_COMM_NULL};
  AutoPas_MPI_Request _request{AUTOPAS_MPI_REQUEST_NULL};

  /**
   * The globally optimal configuration.
   * Usually holds a value that is not in a given rank's search space
   */
  Configuration _optimalConfiguration{Configuration()};

  bool _allLocalConfigurationsTested{false};
  bool _allGlobalConfigurationsTested{false};
  bool _strategyStillWorking{true};

  // fallback configurations, in case the underlying search strategy fails.
  const std::set<ContainerOption> _fallbackContainers;
  // fallback options cannot deal with continuous cellSizeFactors.
  const autopas::NumberSetFinite<double> _fallbackCellSizeFactor{1};
  const std::set<TraversalOption> _fallbackTraversalOptions;
  const std::set<LoadEstimatorOption> _fallbackLoadEstimators;
  const std::set<DataLayoutOption> _fallbackDataLayouts;
  const std::set<Newton3Option> _fallbackNewton3s;
  int _numFallbackConfigs{-1};
  std::unique_ptr<utils::ConfigurationAndRankIteratorHandler> _configIterator{nullptr};
};
}  // namespace autopas
