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
#include "autopas/selectors/OptimumSelector.h"
#include "autopas/selectors/TraversalSelector.h"
#include "autopas/selectors/tuningStrategy/TuningStrategyInterface.h"

namespace autopas {

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
   * @param cutoff Cutoff radius to be used in this container.
   * @param cellSizeFactor Cell size factor to be used in this container (only relevant for LinkedCells).
   * @param verletSkin Length added to the cutoff for the Verlet lists' skin.
   * @param verletRebuildFrequency Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   * @param tuningStrategy Object implementing the modelling and exploration of a search space.
   * @param selectorStrategy Strategy for the configuration selection.
   * @param tuningInterval Number of time steps after which the auto-tuner shall reevaluate all selections.
   * @param maxSamples Number of samples that shall be collected for each combination.
   */
  AutoTuner(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff, double cellSizeFactor,
            double verletSkin, unsigned int verletRebuildFrequency, TuningStrategyInterface *tuningStrategy,
            SelectorStrategyOption selectorStrategy, unsigned int tuningInterval, unsigned int maxSamples)
      : _selectorStrategy(selectorStrategy),
        _tuningStrategy(tuningStrategy),
        _tuningInterval(tuningInterval),
        _iterationsSinceTuning(tuningInterval),  // init to max so that tuning happens in first iteration
        _containerSelector(boxMin, boxMax, cutoff, cellSizeFactor, verletSkin, verletRebuildFrequency),
        _maxSamples(maxSamples),
        _samples(maxSamples) {
    if (_tuningStrategy->searchSpaceEmpty()) {
      autopas::utils::ExceptionHandler::exception("AutoTuner: Passed tuning strategy has an empty search space.");
    }

    // collect all potentially needed traversal selector infos
    for (auto &containerOption : _tuningStrategy->getAllowedContainerOptions()) {
      _containerSelector.selectContainer(containerOption);
      _traversalSelectorInfos.emplace(containerOption,
                                      _containerSelector.getCurrentContainer()->getTraversalSelectorInfo());
    }

    _containerSelector.selectContainer(tuningStrategy->getCurrentConfiguration()._container);
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
    if (_tuningStrategy->searchSpaceOneOption()) {
      return false;
    }

    return _iterationsSinceTuning >= _tuningInterval and _samples.size() >= _maxSamples;
  }

  /**
   * Save the runtime of a given traversal if the functor is relevant for tuning.
   *
   * Samples are collected and reduced to one single value according to _selectorStrategy. Only then the value is passed
   * on to the tuning strategy. This function expects that samples of the same configuration are taken consecutively.
   * The time argument is a long because std::chrono::duration::count returns a long.
   *
   * @param pairwiseFunctor
   * @param time
   */
  template <class PairwiseFunctor>
  void addTimeMeasurement(PairwiseFunctor &pairwiseFunctor, long time) {
    if (pairwiseFunctor.isRelevantForTuning()) {
      if (_samples.size() < _maxSamples) {
        AutoPasLog(trace, "Adding sample.");
        _samples.push_back(time);
        // if this was the last sample:
        if (_samples.size() == _maxSamples) {
          auto reducedValue = OptimumSelector::optimumValue(_samples, _selectorStrategy);
          _tuningStrategy->addEvidence(time);

          // print config, times and reduced value
          if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
            std::stringstream ss;
            // print all configs
            ss << std::endl << _tuningStrategy->getCurrentConfiguration().toString() << " : [";
            // print all timings
            for (auto &sample : _samples) {
              ss << " " << sample;
            }
            ss << " ] ";
            ss << "Reduced value: " << reducedValue;
            AutoPasLog(debug, "Collected times for  {}", ss.str());
          }
        }
      }
    } else {
      AutoPasLog(trace, "Skipping adding of time measurement because functor is not marked relevant.");
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

  SelectorStrategyOption _selectorStrategy;
  TuningStrategyInterface *_tuningStrategy;
  unsigned int _tuningInterval, _iterationsSinceTuning;
  ContainerSelector<Particle, ParticleCell> _containerSelector;

  /**
   * One selector / factory per possible container because every container might have a different number of cells.
   */
  std::map<ContainerOption, TraversalSelectorInfo<ParticleCell>> _traversalSelectorInfos;

  /**
   * How many times each configuration should be tested.
   */
  const size_t _maxSamples;
  /**
   * How many times this configurations has already been tested.
   * Initialize with max value to start tuning at start of simulation.
   */
  std::vector<size_t> _samples;
};

template <class Particle, class ParticleCell>
template <class PairwiseFunctor>
bool AutoTuner<Particle, ParticleCell>::iteratePairwise(PairwiseFunctor *f) {
  bool isTuning = false;
  // tune if :
  // - more than one config exists
  // - currently in tuning phase
  // - functor is relevant
  if ((not _tuningStrategy->searchSpaceOneOption()) and _iterationsSinceTuning >= _tuningInterval and
      f->isRelevantForTuning()) {
    isTuning = tune<PairwiseFunctor>(*f);
    if (not isTuning) {
      _iterationsSinceTuning = 0;
    }
  }

  // large case differentiation for data layout and newton 3
  switch (_tuningStrategy->getCurrentConfiguration()._dataLayout) {
    case DataLayoutOption::aos: {
      if (_tuningStrategy->getCurrentConfiguration()._newton3 == Newton3Option::enabled) {
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
      if (_tuningStrategy->getCurrentConfiguration()._newton3 == Newton3Option::enabled) {
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
      if (_tuningStrategy->getCurrentConfiguration()._newton3 == Newton3Option::enabled) {
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
      utils::ExceptionHandler::exception("AutoTuner: Unknown data layout : {}",
                                         _tuningStrategy->getCurrentConfiguration()._dataLayout);
  }

  if (f->isRelevantForTuning()) {
    ++_iterationsSinceTuning;
  }
  return isTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3, bool inTuningPhase>
void AutoTuner<Particle, ParticleCell>::iteratePairwiseTemplateHelper(PairwiseFunctor *f) {
  auto container = getContainer();
  AutoPasLog(debug, "Iterating with configuration: {}", _tuningStrategy->getCurrentConfiguration().toString());

  auto traversal = TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayout, useNewton3>(
      _tuningStrategy->getCurrentConfiguration()._traversal, *f,
      _traversalSelectorInfos[_tuningStrategy->getCurrentConfiguration()._container]);

  // if tuning execute with time measurements
  if (inTuningPhase) {
    auto start = std::chrono::high_resolution_clock::now();
    // @todo remove useNewton3 in iteratePairwise by introducing traversals for DS and VL

    f->initTraversal();
    withStaticContainerType(container,
                            [&](auto container) { container->iteratePairwise(f, traversal.get(), useNewton3); });
    f->endTraversal(useNewton3);

    auto stop = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    AutoPasLog(debug, "IteratePairwise took {} nanoseconds", runtime);
    //    _containerSelector.addTimeMeasurement(container->getContainerType(), runtime);
    //    traversalSelector.addTimeMeasurement(*f, traversal->getTraversalType(), runtime);
    addTimeMeasurement(*f, runtime);
  } else {
    f->initTraversal();
    withStaticContainerType(container,
                            [&](auto container) { container->iteratePairwise(f, traversal.get(), useNewton3); });
    f->endTraversal(useNewton3);
  }
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor>
bool AutoTuner<Particle, ParticleCell>::tune(PairwiseFunctor &pairwiseFunctor) {
  bool stillTuning = true;

  // need more samples; keep current config
  if (_samples.size() < _maxSamples) {
    return stillTuning;
  }

  // first tuning iteration -> reset to first config
  if (_iterationsSinceTuning == _tuningInterval) {
    _tuningStrategy->reset();
  } else {  // enough samples -> next config
    stillTuning = _tuningStrategy->tune();
  }
  // @TODO move to addMeasurement?
  _samples.clear();

  // repeat as long as traversals are not applicable or we run out of configs
  while (true) {
    // check if newton3 works with this functor and remove config if not
    if ((_tuningStrategy->getCurrentConfiguration()._newton3 == Newton3Option::enabled and
         not pairwiseFunctor.allowsNewton3()) or
        (_tuningStrategy->getCurrentConfiguration()._newton3 == Newton3Option::disabled and
         not pairwiseFunctor.allowsNonNewton3())) {
      AutoPasLog(warn, "Configuration with newton 3 {} called with a functor that does not support this!",
                 utils::StringUtils::to_string(_tuningStrategy->getCurrentConfiguration()._newton3));

      _tuningStrategy->removeN3Option(_tuningStrategy->getCurrentConfiguration()._newton3);
    } else {
      if (configApplicable(_tuningStrategy->getCurrentConfiguration(), pairwiseFunctor)) {
        // we found a valid config!
        break;
      } else {
        stillTuning = _tuningStrategy->tune();
      }
    }
  }

  _containerSelector.selectContainer(_tuningStrategy->getCurrentConfiguration()._container);
  return stillTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor>
bool AutoTuner<Particle, ParticleCell>::configApplicable(const Configuration &conf, PairwiseFunctor &pairwiseFunctor) {
  bool traversalApplicable = false;

  switch (conf._dataLayout) {
    case DataLayoutOption::aos: {
      switch (conf._newton3) {
        case Newton3Option::enabled: {
          traversalApplicable =
              TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::aos, true>(
                  conf._traversal, pairwiseFunctor, _traversalSelectorInfos[conf._container])
                  ->isApplicable();
          break;
        }
        case Newton3Option::disabled: {
          traversalApplicable =
              TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::aos,
                                                                          false>(
                  conf._traversal, pairwiseFunctor, _traversalSelectorInfos[conf._container])
                  ->isApplicable();
          break;
        }
      }
      break;
    }
    case DataLayoutOption::soa: {
      switch (conf._newton3) {
        case Newton3Option::enabled: {
          traversalApplicable =
              TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::soa, true>(
                  conf._traversal, pairwiseFunctor, _traversalSelectorInfos[conf._container])
                  ->isApplicable();
          break;
        }
        case Newton3Option::disabled: {
          traversalApplicable =
              TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::soa,
                                                                          false>(
                  conf._traversal, pairwiseFunctor, _traversalSelectorInfos[conf._container])
                  ->isApplicable();
          break;
        }
      }
      break;
    }
    case DataLayoutOption::cuda: {
      switch (conf._newton3) {
        case Newton3Option::enabled: {
          traversalApplicable =
              TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::cuda,
                                                                          true>(
                  conf._traversal, pairwiseFunctor, _traversalSelectorInfos[conf._container])
                  ->isApplicable();
          break;
        }
        case Newton3Option::disabled: {
          traversalApplicable =
              TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayoutOption::cuda,
                                                                          false>(
                  conf._traversal, pairwiseFunctor, _traversalSelectorInfos[conf._container])
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
  return _tuningStrategy->getCurrentConfiguration();
}
}  // namespace autopas
