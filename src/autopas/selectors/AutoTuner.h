/**
 * @file AutoTuner.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>
#include <memory>
#include <set>
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
   * Constructor for the AutoTuner
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff  Cutoff radius to be used in this container.
   * @param verletSkin Length added to the cutoff for the verlet lists' skin.
   * @param verletRebuildFrequency Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   * @param allowedContainerOptions Vector of container types AutoPas can choose from.
   * @param allowedTraversalOptions Vector of traversals AutoPas can choose from.
   * @param selectorStrategy Strategy for the configuration selection.
   * @param tuningInterval Number of timesteps after which the auto-tuner shall reevaluate all selections.
   * @param maxSamples Number of samples that shall be collected for each combination.
   */
  AutoTuner(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff, double verletSkin,
            unsigned int verletRebuildFrequency, std::vector<ContainerOption> allowedContainerOptions,
            std::vector<TraversalOption> allowedTraversalOptions,
            std::vector<DataLayoutOption> allowedDataLayoutOptions, std::vector<Newton3Option> allowedNewton3Options,
            SelectorStrategy selectorStrategy, unsigned int tuningInterval, unsigned int maxSamples)
      : _tuningInterval(tuningInterval),
        _iterationsSinceTuning(tuningInterval),  // init to max so that tuning happens in first iteration
        _containerSelector(boxMin, boxMax, cutoff, verletSkin, verletRebuildFrequency, allowedContainerOptions,
                           allowedTraversalOptions),
        _allowedTraversalOptions(allowedTraversalOptions),
        _maxSamples(maxSamples),
        _numSamples(maxSamples),
        _selectorStrategy(selectorStrategy),
        _allowedConfigurations(),
        _currentConfig(),
        _traversalTimes() {
    for (auto &containerOption : allowedContainerOptions) {
      for (auto &traversalOption : allowedTraversalOptions) {
        // ToDo: throw out traversal options not applicable to the current container
        for (auto &dataLayoutOption : allowedDataLayoutOptions) {
          for (auto &newton3Option : allowedNewton3Options) {
            _allowedConfigurations.emplace(containerOption, traversalOption, dataLayoutOption, newton3Option);
          }
        }
      }
    }
    _currentConfig = _allowedConfigurations.end();
  }

  /**
   * Getter for the optimal container.
   * Also checks if the container was already encountered and if not creates a new traversal selector for it.
   * @return Smartpointer to the optimal container.
   */
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> getContainer() {
    auto container = _containerSelector.getCurrentContainer();
    // if the container is new create a new traversal selector for it
    // @TODO probably not necessary anymore -> remove
    //    if (_traversalSelectors.find(container->getContainerType()) == _traversalSelectors.end()) {
    //      _traversalSelectors.insert(std::make_pair(container->getContainerType(),
    //                                                container->generateTraversalSelector(_allowedTraversalOptions)));
    //    }
    return container;
  };

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
   * Returns whether the container or the traversal will be changed in the next iteration.
   * @return True if the container will be rebuild on the next iteratePairwise() call. False otherwise.
   */
  bool willRebuild() {
    if (_iterationsSinceTuning >= _tuningInterval) {
      if (_numSamples < _maxSamples) {
        return false;
      } else {
        return (_currentConfig->_container != std::next(_currentConfig)->_container) or
               (_currentConfig->_traversal != std::next(_currentConfig)->_traversal) or
               (_currentConfig == _allowedConfigurations.end());
      }
    } else {
      return false;
    }
  }

  /**
   * Save the runtime of a given traversal if the functor is relevant for tuning.
   * The time argument is a long because std::chrono::duration::count returns a long
   * @param pairwiseFunctor
   * @param traversal
   * @param time
   */
  template <class PairwiseFunctor>
  void addTimeMeasurement(PairwiseFunctor &pairwiseFunctor, long time) {
    if (pairwiseFunctor.isRelevantForTuning()) {
      _traversalTimes[*_currentConfig].push_back(time);
    }
  }

 private:
  void selectOptimalConfiguration();

  template <class PairwiseFunctor, bool useSoA, bool useNewton3, bool inTuningPhase>
  void iteratePairwiseTemplateHelper(PairwiseFunctor *f);

  template <class PairwiseFunctor>
  bool tune(PairwiseFunctor &pairwiseFunctor);

  unsigned int _tuningInterval, _iterationsSinceTuning;
  ContainerSelector<Particle, ParticleCell> _containerSelector;
  std::vector<TraversalOption> _allowedTraversalOptions;
  /**
   * One selector per possible container.
   */
  //  std::map<ContainerOption, TraversalSelector<ParticleCell>> _traversalSelectors;
  TraversalSelector<ParticleCell> _traversalSelectors;

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
  // check if currently in tuning phase, execute iteration and take time measurement if necessary
  if (_iterationsSinceTuning >= _tuningInterval) {
    isTuning = tune<PairwiseFunctor>(*f);
  }

  // large case differentiation for data layout and newton 3
  switch (_currentConfig->_dataLayout) {
    case DataLayoutOption::aos: {
      if (_currentConfig->_newton3) {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, false, true, true>(f);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, false, true, false>(f);
        }
      } else {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, false, false, true>(f);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, false, false, false>(f);
        }
      }
      break;
    }
    case DataLayoutOption::soa: {
      if (_currentConfig->_newton3) {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, true, true, true>(f);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, true, true, false>(f);
        }
      } else {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, true, false, true>(f);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, true, false, false>(f);
        }
      }
      break;
    }
    default:
      utils::ExceptionHandler::exception("AutoTuner: Unknown data layout : {}", _currentConfig->_dataLayout);
  }
  return isTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3, bool inTuningPhase>
void AutoTuner<Particle, ParticleCell>::iteratePairwiseTemplateHelper(PairwiseFunctor *f) {
  auto container = getContainer();
  AutoPasLog(debug, "Using container {}", utils::StringUtils::to_string(container->getContainerType()));

  auto traversal = _traversalSelectors.template getCurrentTraversal<PairwiseFunctor, useSoA, useNewton3>(*f);

  // if tuning execute with time measurements
  if (inTuningPhase) {
    auto start = std::chrono::high_resolution_clock::now();
    // @todo remove useNewton3 in iteratePairwise by introducing traversals for DS and VL
    if (useSoA) {
      withStaticContainerType(container,
                              [&](auto container) { container->iteratePairwiseSoA(f, traversal.get(), useNewton3); });
    } else {
      withStaticContainerType(container,
                              [&](auto container) { container->iteratePairwiseAoS(f, traversal.get(), useNewton3); });
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    AutoPasLog(debug, "IteratePairwise took {} nanoseconds", runtime);
    //    _containerSelector.addTimeMeasurement(container->getContainerType(), runtime);
    //    traversalSelector.addTimeMeasurement(*f, traversal->getTraversalType(), runtime);
    addTimeMeasurement(*f, runtime);
  } else {
    if (useSoA) {
      withStaticContainerType(container,
                              [&](auto container) { container->iteratePairwiseSoA(f, traversal.get(), useNewton3); });
    } else {
      withStaticContainerType(container,
                              [&](auto container) { container->iteratePairwiseAoS(f, traversal.get(), useNewton3); });
    }
  }
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor>
bool AutoTuner<Particle, ParticleCell>::tune(PairwiseFunctor &pairwiseFunctor) {
  // check for sane newton3 settings and fix them if needed
  if ((_currentConfig->_newton3 == Newton3Option::enabled and not pairwiseFunctor.allowsNewton3()) or
      (_currentConfig->_newton3 == Newton3Option::disabled and not pairwiseFunctor.allowsNonNewton3())) {
    AutoPasLog(warn, "Configuration with newton 3 {} called with a functor that does not support this!",
               utils::StringUtils::to_string(_currentConfig->_newton3));

    Configuration modifiedCurrentConfig = *_currentConfig;
    // choose the other option
    modifiedCurrentConfig._newton3 = Newton3Option((static_cast<int>(_currentConfig->_newton3) + 1) % 2);

    // if modified config is equal to next delete current. Else insert modified config.
    // the disabled case should come after the enabled.
    int searchDirection = _currentConfig->_newton3 == Newton3Option::enabled ? 1 : -1;
    if (modifiedCurrentConfig == *(std::next(_currentConfig, searchDirection))) {
      AutoPasLog(warn, "Newton 3 {} automatically for this config!",
                 utils::StringUtils::to_string(modifiedCurrentConfig._newton3));
      if (pairwiseFunctor.isRelevantForTuning()) {
        _currentConfig = _allowedConfigurations.erase(_currentConfig);
      }
    } else {
      _currentConfig = _allowedConfigurations.insert(_currentConfig, modifiedCurrentConfig);
    }
  }

  bool stillTuning = true;

  // if beginning of tuning phase
  if (_iterationsSinceTuning == _tuningInterval) {
    _currentConfig = _allowedConfigurations.begin();
    _numSamples = 0;
  } else {
    if (_numSamples >= _maxSamples) {
      _numSamples = 0;
      ++_currentConfig;
    }
  }

  // reached end of tuning phase
  if (_currentConfig == _allowedConfigurations.end()) {
    selectOptimalConfiguration();
    _iterationsSinceTuning = 0;
    stillTuning = false;
  }

  _containerSelector.selectContainer(_currentConfig->_container);
  _traversalSelectors.selectTraversal(_currentConfig->_traversal);

  ++_iterationsSinceTuning;
  ++_numSamples;
  return stillTuning;
}

template <class Particle, class ParticleCell>
void AutoTuner<Particle, ParticleCell>::selectOptimalConfiguration() {
  // Time measure strategy
  if (_traversalTimes.empty()) {
    utils::ExceptionHandler::exception("AutoTuner: Trying to determine fastest configuration before measuring!");
  }

  long optimalTraversalTime = std::numeric_limits<long>::max();
  Configuration optimalConfiguration;
  // reduce sample values
  for (auto &configAndTimes : _traversalTimes) {
    long value = 0;
    switch (_selectorStrategy) {
      case SelectorStrategy::fastestAbs: {
        value = *std::min_element(configAndTimes.second.begin(), configAndTimes.second.end());
        break;
      }
      case SelectorStrategy::fastestMean: {
        value = std::accumulate(configAndTimes.second.begin(), configAndTimes.second.end(), 0l);
        value /= configAndTimes.second.size();
        break;
      }
      case SelectorStrategy::fastestMedian: {
        std::sort(configAndTimes.second.begin(), configAndTimes.second.end());
        value = configAndTimes.second[configAndTimes.second.size() / 2];
        break;
      }
      default:
        utils::ExceptionHandler::exception("AutoTuner: Unknown selector strategy {}", _selectorStrategy);
    }
    // save all values for debugging purposes
    if (spdlog::get("AutoPasLog")->level() <= spdlog::level::debug) {
      configAndTimes.second.push_back(value);
    }

    if (value < optimalTraversalTime) {
      optimalTraversalTime = value;
      optimalConfiguration = configAndTimes.first;
    }
  }

  _currentConfig = _allowedConfigurations.find(optimalConfiguration);

  // print all configs, times and their reduced values
  AutoPasLog(debug, "Collected times: {}", [&]() -> std::string {
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
    return ss.str();
  }());  // () to evaluate the function here.

  // measurements are not needed anymore
  _traversalTimes.clear();

  AutoPasLog(debug, "Selected Configuration {}", _currentConfig->toString());
}
}  // namespace autopas
