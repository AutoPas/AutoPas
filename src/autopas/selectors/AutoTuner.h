/**
 * @file AutoTuner.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>
#include <memory>
#include <set>
#include "OptimumSelector.h"
#include "autopas/autopasIncludes.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/selectors/Configuration.h"
#include "autopas/selectors/ContainerSelector.h"
#include "autopas/selectors/TraversalSelector.h"

namespace autopas {

/**
 * Possible choices for the auto tuner.
 * @todo: implement more options and then use this enum! :D
 */
enum TuningStrategy {
  /**
   * Test all allowed configurations and select the best.
   */
  timeMeasuring
};

/**
 * Automated tuner for optimal iteration performance.
 *
 * This class offers an interface to the iteratePairwise method.
 * Internally it chooses the best container, traversal etc for the current simulation.
 *
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle, class ParticleCell>
class AutoTuner {
 public:
  /**
   * Constructor for the AutoTuner that generates all configurations from the given options.
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff  Cutoff radius to be used in this container.
   * @param verletSkin Length added to the cutoff for the verlet lists' skin.
   * @param verletRebuildFrequency Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   * @param allowedContainerOptions Vector of container types AutoPas can choose from.
   * @param allowedTraversalOptions Vector of traversals AutoPas can choose from.
   * @param allowedDataLayoutOptions Vector of data layouts AutoPas can choose from.
   * @param allowedNewton3Options Vector of newton 3 options AutoPas can choose from.
   * @param selectorStrategy Strategy for the configuration selection.
   * @param tuningInterval Number of timesteps after which the auto-tuner shall reevaluate all selections.
   * @param maxSamples Number of samples that shall be collected for each combination.
   */
  AutoTuner(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff, double verletSkin,
            unsigned int verletRebuildFrequency, const std::vector<ContainerOption> &allowedContainerOptions,
            std::vector<TraversalOption> allowedTraversalOptions,
            const std::vector<DataLayoutOption> &allowedDataLayoutOptions,
            const std::vector<Newton3Option> &allowedNewton3Options, SelectorStrategy selectorStrategy,
            unsigned int tuningInterval, unsigned int maxSamples)
      : _tuningInterval(tuningInterval),
        _iterationsSinceTuning(tuningInterval),  // init to max so that tuning happens in first iteration
        _containerSelector(boxMin, boxMax, cutoff, verletSkin, verletRebuildFrequency),
        _maxSamples(maxSamples),
        _numSamples(maxSamples),
        _selectorStrategy(selectorStrategy),
        _allowedConfigurations(),
        _currentConfig(),
        _traversalTimes() {
    //@TODO needed until all containers support propper traversals
    allowedTraversalOptions.push_back(TraversalOption::dummyTraversal);

    // needs to be sorted for intersection with applicable traversals
    std::sort(allowedTraversalOptions.begin(), allowedTraversalOptions.end());

    // generate all potential configs
    for (auto &containerOption : allowedContainerOptions) {
      _containerSelector.selectContainer(containerOption);
      _traversalSelectors.emplace(containerOption,
                                  _containerSelector.getCurrentContainer()->generateTraversalSelector());

      // get all traversals of the container and restrict them to the allowed ones
      auto allContainerTraversals = _containerSelector.getCurrentContainer()->getAllTraversals();
      std::vector<TraversalOption> allowedAndApplicable;
      std::sort(allContainerTraversals.begin(), allContainerTraversals.end());
      std::set_intersection(allowedTraversalOptions.begin(), allowedTraversalOptions.end(),
                            allContainerTraversals.begin(), allContainerTraversals.end(),
                            std::back_inserter(allowedAndApplicable));

      for (auto &traversalOption : allowedAndApplicable) {
        for (auto &dataLayoutOption : allowedDataLayoutOptions) {
          for (auto &newton3Option : allowedNewton3Options) {
            _allowedConfigurations.emplace(containerOption, traversalOption, dataLayoutOption, newton3Option);
          }
        }
      }
    }

    if (_allowedConfigurations.empty()) {
      autopas::utils::ExceptionHandler::exception("AutoTuner: No valid configurations could be created.");
    }

    _currentConfig = _allowedConfigurations.begin();
    _containerSelector.selectContainer(_currentConfig->_container);
  }

  /**
   * Constructor for the AutoTuner that only contains the given configurations.
   * This constructor assumes only valid configurations are passed! Mainly for easier unit testing.
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff  Cutoff radius to be used in this container.
   * @param verletSkin Length added to the cutoff for the verlet lists' skin.
   * @param verletRebuildFrequency Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   * @param allowedConfigurations Set of configurations AutoPas can choose from.
   * @param selectorStrategy Strategy for the configuration selection.
   * @param tuningInterval Number of timesteps after which the auto-tuner shall reevaluate all selections.
   * @param maxSamples Number of samples that shall be collected for each combination.
   */
  AutoTuner(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff, double verletSkin,
            unsigned int verletRebuildFrequency, const std::set<Configuration> &allowedConfigurations,
            SelectorStrategy selectorStrategy, unsigned int tuningInterval, unsigned int maxSamples)
      : _tuningInterval(tuningInterval),
        _iterationsSinceTuning(tuningInterval),  // init to max so that tuning happens in first iteration
        _containerSelector(boxMin, boxMax, cutoff, verletSkin, verletRebuildFrequency),
        _maxSamples(maxSamples),
        _numSamples(maxSamples),
        _selectorStrategy(selectorStrategy),
        _allowedConfigurations(allowedConfigurations),
        _currentConfig(),
        _traversalTimes() {
    if (_allowedConfigurations.empty()) {
      autopas::utils::ExceptionHandler::exception("AutoTuner: No configurations passed.");
    }

    for (auto &conf : allowedConfigurations) {
      if (_traversalSelectors.find(conf._container) == _traversalSelectors.end()) {
        _containerSelector.selectContainer(conf._container);
        _traversalSelectors.emplace(conf._container,
                                    _containerSelector.getCurrentContainer()->generateTraversalSelector());
      }
    }

    _currentConfig = _allowedConfigurations.begin();
    _containerSelector.selectContainer(_currentConfig->_container);
  }

  /**
   * Getter for the current container.
   * @return Smart pointer to the current container.
   */
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> getContainer() {
    return _containerSelector.getCurrentContainer();
  }

  /**
   * Check if a configuration is applicable to the current domain with the given functor.
   * @tparam PairwiseFunctor
   * @param conf
   * @param pairwiseFunctor
   * @return
   */
  template <class PairwiseFunctor>
  bool configApplicable(const Configuration &conf, PairwiseFunctor &pairwiseFunctor);

  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @tparam PairwiseFunctor
   * @param f Functor that describes the pair-potential.
   * @return true if this was a tuning iteration.
   */
  template <class PairwiseFunctor>
  bool iteratePairwise(PairwiseFunctor *f);

  /**
   * Returns whether the configuration will be changed in the next iteration.
   * This does does not necessarily mean that the container will change.
   *
   * @return True if the next iteratePairwise() call uses a different configuration. False otherwise.
   */
  bool willRebuild() {
    if (_allowedConfigurations.size() == 1) {
      return false;
    }

    return _iterationsSinceTuning >= _tuningInterval and _numSamples >= _maxSamples;
  }

  /**
   * Save the runtime of a given traversal if the functor is relevant for tuning.
   * The time argument is a long because std::chrono::duration::count returns a long
   * @param pairwiseFunctor
   * @param time
   */
  template <class PairwiseFunctor>
  void addTimeMeasurement(PairwiseFunctor &pairwiseFunctor, long time) {
    if (pairwiseFunctor.isRelevantForTuning()) {
      _traversalTimes[*_currentConfig].push_back(time);
    }
  }

  /**
   * Get the currently selected configuration.
   * @return
   */
  const autopas::Configuration getCurrentConfig() const;

  /**
   * Get the set of all allowed configurations.
   * @return
   */
  const std::set<Configuration> &getAllowedConfigurations() const;

 private:
  void selectOptimalConfiguration();

  template <class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3, bool inTuningPhase>
  void iteratePairwiseTemplateHelper(PairwiseFunctor *f);

  /**
   * Tune available algorithm configurations.
   *
   * It is assumed this function is only called for relevant functors and that at least two configurations are allowed.
   * When in tuning phase selects next config to test. At the end of the tuning phase select optimum.
   * The function returns true if the selected config is not yet the optimum but something that should be sampled.
   *
   * @tparam PairwiseFunctor
   * @param pairwiseFunctor
   * @return true iff still in tuning phase.
   */
  template <class PairwiseFunctor>
  bool tune(PairwiseFunctor &pairwiseFunctor);

  unsigned int _tuningInterval, _iterationsSinceTuning;
  ContainerSelector<Particle, ParticleCell> _containerSelector;

  /**
   * One selector / factory per possible container because every container might have a different number of cells.
   */
  std::map<ContainerOption, TraversalSelector<ParticleCell>> _traversalSelectors;

  /**
   * How many times each configuration should be tested.
   */
  const size_t _maxSamples;
  /**
   * How many times this configurations has already been tested.
   * Initialize with max value to start tuning at start of simulation.
   */
  size_t _numSamples;

  SelectorStrategy _selectorStrategy;

  std::set<Configuration> _allowedConfigurations;
  std::set<Configuration>::iterator _currentConfig;
  std::unordered_map<Configuration, std::vector<size_t>, ConfigHash> _traversalTimes;
};

template <class Particle, class ParticleCell>
template <class PairwiseFunctor>
bool AutoTuner<Particle, ParticleCell>::iteratePairwise(PairwiseFunctor *f) {
  bool isTuning = false;
  // tune if :
  // - more than one config exists
  // - currently in tuning phase
  // - functor is relevant
  if (_allowedConfigurations.size() > 1 and _iterationsSinceTuning >= _tuningInterval and f->isRelevantForTuning()) {
    isTuning = tune<PairwiseFunctor>(*f);
    if (not isTuning) {
      _iterationsSinceTuning = 0;
    }
  }

  // large case differentiation for data layout and newton 3
  switch (_currentConfig->_dataLayout) {
    case DataLayoutOption::aos: {
      if (_currentConfig->_newton3 == Newton3Option::enabled) {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ true,
                                        /*tuning*/ true>(f);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ true,
                                        /*tuning*/ false>(f);
        }
      } else {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ false,
                                        /*tuning*/ true>(f);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ false,
                                        /*tuning*/ false>(f);
        }
      }
      break;
    }
    case DataLayoutOption::soa: {
      if (_currentConfig->_newton3 == Newton3Option::enabled) {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ true,
                                        /*tuning*/ true>(f);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ true,
                                        /*tuning*/ false>(f);
        }
      } else {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ false,
                                        /*tuning*/ true>(f);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ false,
                                        /*tuning*/ false>(f);
        }
      }
      break;
    }
#if defined(AUTOPAS_CUDA)
    case DataLayoutOption::cuda: {
      if (_currentConfig->_newton3 == Newton3Option::enabled) {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::cuda, /*Newton3*/ true,
                                        /*tuning*/ true>(f);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::cuda, /*Newton3*/ true,
                                        /*tuning*/ false>(f);
        }
      } else {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::cuda, /*Newton3*/ false,
                                        /*tuning*/ true>(f);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::cuda, /*Newton3*/ false,
                                        /*tuning*/ false>(f);
        }
      }
      break;
    }
#endif
    default:
      utils::ExceptionHandler::exception("AutoTuner: Unknown data layout : {}", _currentConfig->_dataLayout);
  }

  if (f->isRelevantForTuning()) {
    ++_numSamples;
    ++_iterationsSinceTuning;
  }
  return isTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3, bool inTuningPhase>
void AutoTuner<Particle, ParticleCell>::iteratePairwiseTemplateHelper(PairwiseFunctor *f) {
  auto container = getContainer();
  AutoPasLog(debug, "Iterating with configuration: {}", _currentConfig->toString());

  auto traversal =
      _traversalSelectors[_currentConfig->_container]
          .template generateTraversal<PairwiseFunctor, DataLayout, useNewton3>(_currentConfig->_traversal, *f);

  // if tuning execute with time measurements
  if (inTuningPhase) {
    auto start = std::chrono::high_resolution_clock::now();
    // @todo remove useNewton3 in iteratePairwise by introducing traversals for DS and VL
    switch (DataLayout) {
      case DataLayoutOption::aos:
        withStaticContainerType(container,
                                [&](auto container) { container->iteratePairwiseAoS(f, traversal.get(), useNewton3); });
        break;
      case DataLayoutOption::soa:
        withStaticContainerType(container,
                                [&](auto container) { container->iteratePairwiseSoA(f, traversal.get(), useNewton3); });
        break;
      case DataLayoutOption::cuda:
        withStaticContainerType(
            container, [&](auto container) { container->iteratePairwiseSoACuda(f, traversal.get(), useNewton3); });
        break;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    AutoPasLog(debug, "IteratePairwise took {} nanoseconds", runtime);
    //    _containerSelector.addTimeMeasurement(container->getContainerType(), runtime);
    //    traversalSelector.addTimeMeasurement(*f, traversal->getTraversalType(), runtime);
    addTimeMeasurement(*f, runtime);
  } else {
    switch (DataLayout) {
      case DataLayoutOption::aos:
        withStaticContainerType(container,
                                [&](auto container) { container->iteratePairwiseAoS(f, traversal.get(), useNewton3); });
        break;
      case DataLayoutOption::soa:
        withStaticContainerType(container,
                                [&](auto container) { container->iteratePairwiseSoA(f, traversal.get(), useNewton3); });
        break;
      case DataLayoutOption::cuda:
        withStaticContainerType(
            container, [&](auto container) { container->iteratePairwiseSoACuda(f, traversal.get(), useNewton3); });
        break;
    }
  }
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor>
bool AutoTuner<Particle, ParticleCell>::tune(PairwiseFunctor &pairwiseFunctor) {
  bool stillTuning = true;

  // need more samples; keep current config
  if (_numSamples < _maxSamples) {
    return stillTuning;
  }

  // first tuning iteration -> reset to first config
  if (_iterationsSinceTuning == _tuningInterval) {
    _currentConfig = _allowedConfigurations.begin();
    _numSamples = 0;
  } else {  // enough samples -> next config
    ++_currentConfig;
    _numSamples = 0;
  }

  // repeat as long as traversals are not applicable or we run out of configs
  while (_currentConfig != _allowedConfigurations.end()) {
    // check if newton3 works with this functor and remove config if not
    if ((_currentConfig->_newton3 == Newton3Option::enabled and not pairwiseFunctor.allowsNewton3()) or
        (_currentConfig->_newton3 == Newton3Option::disabled and not pairwiseFunctor.allowsNonNewton3())) {
      AutoPasLog(warn, "Configuration with newton 3 {} called with a functor that does not support this!",
                 utils::StringUtils::to_string(_currentConfig->_newton3));

      _currentConfig = _allowedConfigurations.erase(_currentConfig);
    } else {
      if (configApplicable(*_currentConfig, pairwiseFunctor)) {
        // we found a valid config!
        break;
      } else {
        ++_currentConfig;
      }
    }
  }

  // reached end of tuning phase
  // either wait until last config has enough samples or last config is not applicable
  if (_currentConfig == _allowedConfigurations.end()) {
    if (_traversalTimes.empty()) {
      autopas::utils::ExceptionHandler::exception("AutoTuner: No applicable configurations found!");
    }

    // sets _currentConfig
    selectOptimalConfiguration();
    stillTuning = false;
  }

  _containerSelector.selectContainer(_currentConfig->_container);
  return stillTuning;
}

template <class Particle, class ParticleCell>
void AutoTuner<Particle, ParticleCell>::selectOptimalConfiguration() {
  if (_allowedConfigurations.size() == 1) {
    return;
  }

  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception("AutoTuner: Trying to determine fastest configuration before measuring!");
  }

  long optimalTraversalTime = std::numeric_limits<long>::max();
  Configuration optimalConfiguration;
  // reduce sample values
  for (auto &configAndTimes : _traversalTimes) {
    long value = OptimumSelector::optimumValue(configAndTimes.second, _selectorStrategy);

    // save all values for debugging purposes
    if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
      configAndTimes.second.push_back(value);
    }

    if (value < optimalTraversalTime) {
      optimalTraversalTime = value;
      optimalConfiguration = configAndTimes.first;
    }
  }

  // print all configs, times and their reduced values
  if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
    std::stringstream ss;
    // print all configs
    for (auto &p : _traversalTimes) {
      ss << std::endl << p.first.toString() << " : [";
      // print all timings
      for (size_t i = 0; i < p.second.size() - 1; ++i) {
        ss << " " << p.second[i];
      }
      ss << " ] ";
      ss << "Reduced value: " << p.second[p.second.size() - 1];
    }
    AutoPasLog(debug, "Collected times: {}", ss.str());
  }

  _currentConfig = _allowedConfigurations.find(optimalConfiguration);
  // sanity check
  if (_currentConfig == _allowedConfigurations.end()) {
    autopas::utils::ExceptionHandler::exception(
        "AutoTuner: Optimal configuration not found in list of configurations!");
  }

  // measurements are not needed anymore
  _traversalTimes.clear();

  AutoPasLog(debug, "Selected Configuration {}", _currentConfig->toString());
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor>
bool AutoTuner<Particle, ParticleCell>::configApplicable(const Configuration &conf, PairwiseFunctor &pairwiseFunctor) {
  bool traversalApplicable = false;

  switch (conf._dataLayout) {
    case DataLayoutOption::aos: {
      switch (conf._newton3) {
        case Newton3Option::enabled: {
          traversalApplicable = _traversalSelectors[conf._container]
                                    .template generateTraversal<PairwiseFunctor, DataLayoutOption::aos, true>(
                                        conf._traversal, pairwiseFunctor)
                                    ->isApplicable();
          break;
        }
        case Newton3Option::disabled: {
          traversalApplicable = _traversalSelectors[conf._container]
                                    .template generateTraversal<PairwiseFunctor, DataLayoutOption::aos, false>(
                                        conf._traversal, pairwiseFunctor)
                                    ->isApplicable();
          break;
        }
      }
      break;
    }
    case DataLayoutOption::soa: {
      switch (conf._newton3) {
        case Newton3Option::enabled: {
          traversalApplicable = _traversalSelectors[conf._container]
                                    .template generateTraversal<PairwiseFunctor, DataLayoutOption::soa, true>(
                                        conf._traversal, pairwiseFunctor)
                                    ->isApplicable();
          break;
        }
        case Newton3Option::disabled: {
          traversalApplicable = _traversalSelectors[conf._container]
                                    .template generateTraversal<PairwiseFunctor, DataLayoutOption::soa, false>(
                                        conf._traversal, pairwiseFunctor)
                                    ->isApplicable();
          break;
        }
      }
      break;
    }
    case DataLayoutOption::cuda: {
      switch (conf._newton3) {
        case Newton3Option::enabled: {
          traversalApplicable = _traversalSelectors[conf._container]
                                    .template generateTraversal<PairwiseFunctor, DataLayoutOption::cuda, true>(
                                        conf._traversal, pairwiseFunctor)
                                    ->isApplicable();
          break;
        }
        case Newton3Option::disabled: {
          traversalApplicable = _traversalSelectors[conf._container]
                                    .template generateTraversal<PairwiseFunctor, DataLayoutOption::cuda, false>(
                                        conf._traversal, pairwiseFunctor)
                                    ->isApplicable();
          break;
        }
      }
      break;
    }
  }

  return traversalApplicable;
}

template <class Particle, class ParticleCell>
const autopas::Configuration AutoTuner<Particle, ParticleCell>::getCurrentConfig() const {
  return *_currentConfig;
}

template <class Particle, class ParticleCell>
const std::set<Configuration> &AutoTuner<Particle, ParticleCell>::getAllowedConfigurations() const {
  return _allowedConfigurations;
}
}  // namespace autopas
