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
#include "autopas/utils/ArrayUtils.h"

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
   * @param verletSkin Length added to the cutoff for the Verlet lists' skin.
   * @param verletClusterSize Number of particles in a cluster to use in verlet list.
   * @param tuningStrategy Object implementing the modelling and exploration of a search space.
   * @param selectorStrategy Strategy for the configuration selection.
   * @param tuningInterval Number of time steps after which the auto-tuner shall reevaluate all selections.
   * @param maxSamples Number of samples that shall be collected for each combination.
   */
  AutoTuner(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff, double verletSkin,
            unsigned int verletClusterSize, std::unique_ptr<TuningStrategyInterface> tuningStrategy,
            SelectorStrategyOption selectorStrategy, unsigned int tuningInterval, unsigned int maxSamples)
      : _selectorStrategy(selectorStrategy),
        _tuningStrategy(std::move(tuningStrategy)),
        _tuningInterval(tuningInterval),
        _iterationsSinceTuning(tuningInterval),  // init to max so that tuning happens in first iteration
        _containerSelector(boxMin, boxMax, cutoff),
        _verletSkin(verletSkin),
        _verletClusterSize(verletClusterSize),
        _maxSamples(maxSamples),
        _samples(maxSamples) {
    if (_tuningStrategy->searchSpaceIsEmpty()) {
      autopas::utils::ExceptionHandler::exception("AutoTuner: Passed tuning strategy has an empty search space.");
    }

    selectCurrentContainer();
  }

  /**
   * Move assignment operator
   * @param other
   * @return
   */
  AutoTuner &operator=(AutoTuner &&other) noexcept {
    _tuningStrategy = std::move(other._tuningStrategy);
    return *this;
  }

  /**
   * Getter for the current container.
   * @return Smart pointer to the current container.
   */
  std::shared_ptr<autopas::ParticleContainer<ParticleCell>> getContainer() {
    return _containerSelector.getCurrentContainer();
  }

  /**
   * Getter for the current container.
   * @return Smart pointer to the current container.
   */
  std::shared_ptr<const autopas::ParticleContainer<ParticleCell>> getContainer() const {
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
   * @param doListRebuild Indicates whether or not the verlet lists should be rebuild.
   * @return true if this was a tuning iteration.
   */
  template <class PairwiseFunctor>
  bool iteratePairwise(PairwiseFunctor *f, bool doListRebuild);

  /**
   * Returns whether the configuration will be changed in the next iteration.
   * This does does not necessarily mean that the container will change.
   *
   * @return True if the next iteratePairwise() call uses a different configuration. False otherwise.
   */
  bool willRebuild() {
    if (_tuningStrategy->searchSpaceIsTrivial()) {
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
          _tuningStrategy->addEvidence(reducedValue);

          // print config, times and reduced value
          if (autopas::Logger::get()->level() <= autopas::Logger::LogLevel::debug) {
            std::ostringstream ss;
            // print all configs
            ss << std::endl << _tuningStrategy->getCurrentConfiguration().toString() << " : [ ";
            // print all timings
            ss << utils::ArrayUtils::to_string(_samples, " ");
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
  /**
   * Initialize the container specified by the TuningStrategy.
   */
  void selectCurrentContainer();

  template <class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3, bool inTuningPhase>
  void iteratePairwiseTemplateHelper(PairwiseFunctor *f, bool doListRebuild);

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
  std::unique_ptr<TuningStrategyInterface> _tuningStrategy;
  unsigned int _tuningInterval, _iterationsSinceTuning;
  ContainerSelector<Particle, ParticleCell> _containerSelector;
  double _verletSkin;
  unsigned int _verletClusterSize;

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
void AutoTuner<Particle, ParticleCell>::selectCurrentContainer() {
  auto conf = _tuningStrategy->getCurrentConfiguration();
  _containerSelector.selectContainer(conf.container,
                                     ContainerSelectorInfo(conf.cellSizeFactor, _verletSkin, _verletClusterSize));
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor>
bool AutoTuner<Particle, ParticleCell>::iteratePairwise(PairwiseFunctor *f, bool doListRebuild) {
  bool isTuning = false;
  // tune if :
  // - more than one config exists
  // - currently in tuning phase
  // - functor is relevant
  if ((not _tuningStrategy->searchSpaceIsTrivial()) and _iterationsSinceTuning >= _tuningInterval and
      f->isRelevantForTuning()) {
    isTuning = tune<PairwiseFunctor>(*f);
    if (not isTuning) {
      _iterationsSinceTuning = 0;
    }
  }

  // large case differentiation for data layout and newton 3
  switch (_tuningStrategy->getCurrentConfiguration().dataLayout) {
    case DataLayoutOption::aos: {
      if (_tuningStrategy->getCurrentConfiguration().newton3 == Newton3Option::enabled) {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ true,
                                        /*tuning*/ true>(f, doListRebuild);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ true,
                                        /*tuning*/ false>(f, doListRebuild);
        }
      } else {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ false,
                                        /*tuning*/ true>(f, doListRebuild);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::aos, /*Newton3*/ false,
                                        /*tuning*/ false>(f, doListRebuild);
        }
      }
      break;
    }
    case DataLayoutOption::soa: {
      if (_tuningStrategy->getCurrentConfiguration().newton3 == Newton3Option::enabled) {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ true,
                                        /*tuning*/ true>(f, doListRebuild);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ true,
                                        /*tuning*/ false>(f, doListRebuild);
        }
      } else {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ false,
                                        /*tuning*/ true>(f, doListRebuild);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::soa, /*Newton3*/ false,
                                        /*tuning*/ false>(f, doListRebuild);
        }
      }
      break;
    }
#if defined(AUTOPAS_CUDA)
    case DataLayoutOption::cuda: {
      if (_tuningStrategy->getCurrentConfiguration().newton3 == Newton3Option::enabled) {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::cuda, /*Newton3*/ true,
                                        /*tuning*/ true>(f, doListRebuild);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::cuda, /*Newton3*/ true,
                                        /*tuning*/ false>(f, doListRebuild);
        }
      } else {
        if (isTuning) {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::cuda, /*Newton3*/ false,
                                        /*tuning*/ true>(f, doListRebuild);
        } else {
          iteratePairwiseTemplateHelper<PairwiseFunctor, DataLayoutOption::cuda, /*Newton3*/ false,
                                        /*tuning*/ false>(f, doListRebuild);
        }
      }
      break;
    }
#endif
    default:
      utils::ExceptionHandler::exception("AutoTuner: Unknown data layout : {}",
                                         _tuningStrategy->getCurrentConfiguration().dataLayout);
  }

  if (f->isRelevantForTuning()) {
    ++_iterationsSinceTuning;
  }
  return isTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor, DataLayoutOption DataLayout, bool useNewton3, bool inTuningPhase>
void AutoTuner<Particle, ParticleCell>::iteratePairwiseTemplateHelper(PairwiseFunctor *f, bool doListRebuild) {
  auto containerPtr = getContainer();
  AutoPasLog(debug, "Iterating with configuration: {}", _tuningStrategy->getCurrentConfiguration().toString());

  auto traversal = TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor, DataLayout, useNewton3>(
      _tuningStrategy->getCurrentConfiguration().traversal, *f, containerPtr->getTraversalSelectorInfo());

  if (not traversal->isApplicable()) {
    autopas::utils::ExceptionHandler::exception(
        "Error: Trying to execute a traversal that is not applicable. This normally happens only if the search space "
        "is trivial, but no traversals are applicable. Config: {}",
        _tuningStrategy->getCurrentConfiguration().toString());
  }

  // if tuning execute with time measurements
  if (inTuningPhase) {
    auto start = std::chrono::high_resolution_clock::now();

    f->initTraversal();
    if (doListRebuild) {
      containerPtr->rebuildNeighborLists(traversal.get());
    }
    containerPtr->iteratePairwise(traversal.get());
    f->endTraversal(useNewton3);

    auto stop = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    AutoPasLog(debug, "IteratePairwise took {} nanoseconds", runtime);
    addTimeMeasurement(*f, runtime);
  } else {
    f->initTraversal();
    if (doListRebuild) {
      containerPtr->rebuildNeighborLists(traversal.get());
    }
    containerPtr->iteratePairwise(traversal.get());
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

  // repeat as long as traversals are not applicable or we run out of configs
  while (true) {
    // check if newton3 works with this functor and remove config if not
    if ((_tuningStrategy->getCurrentConfiguration().newton3 == Newton3Option::enabled and
         not pairwiseFunctor.allowsNewton3()) or
        (_tuningStrategy->getCurrentConfiguration().newton3 == Newton3Option::disabled and
         not pairwiseFunctor.allowsNonNewton3())) {
      AutoPasLog(warn, "Configuration with newton 3 {} called with a functor that does not support this!",
                 utils::StringUtils::to_string(_tuningStrategy->getCurrentConfiguration().newton3));

      _tuningStrategy->removeN3Option(_tuningStrategy->getCurrentConfiguration().newton3);
    } else {
      if (configApplicable(_tuningStrategy->getCurrentConfiguration(), pairwiseFunctor)) {
        // we found a valid config!
        break;
      } else {
        AutoPasLog(debug, "Skip not applicable configuration {}",
                   _tuningStrategy->getCurrentConfiguration().toString());
        stillTuning = _tuningStrategy->tune(true);
      }
    }
  }
  // samples should only be cleared if we are still tuning, see `if (_samples.size() < _maxSamples)` from before.
  if (stillTuning) {
    // samples are no longer needed. Delete them here so willRebuild() works as expected.
    _samples.clear();
  }

  selectCurrentContainer();
  return stillTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor>
bool AutoTuner<Particle, ParticleCell>::configApplicable(const Configuration &conf, PairwiseFunctor &pairwiseFunctor) {
  auto allContainerTraversals = compatibleTraversals::allCompatibleTraversals(conf.container);
  if (allContainerTraversals.find(conf.traversal) == allContainerTraversals.end()) {
    // container and traversal mismatch
    return false;
  }

  _containerSelector.selectContainer(conf.container, ContainerSelectorInfo(conf.cellSizeFactor, _verletSkin));
  auto traversalInfo = _containerSelector.getCurrentContainer()->getTraversalSelectorInfo();

  return TraversalSelector<ParticleCell>::template generateTraversal<PairwiseFunctor>(
             conf.traversal, pairwiseFunctor, traversalInfo, conf.dataLayout, conf.newton3)
      ->isApplicable();
}

template <class Particle, class ParticleCell>
const autopas::Configuration AutoTuner<Particle, ParticleCell>::getCurrentConfig() const {
  return _tuningStrategy->getCurrentConfiguration();
}
}  // namespace autopas
