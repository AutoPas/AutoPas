/**
 * @file AutoTuner.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>
#include <memory>
#include "autopas/autopasIncludes.h"
#include "autopas/options/DataLayoutOptions.h"
#include "autopas/options/TraversalOptions.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/selectors/ContainerSelector.h"
#include "autopas/selectors/TraversalSelector.h"

namespace autopas {

/**
 * Provides a way to iterate over the possible choices of data layouts.
 */
static std::vector<DataLayoutOption> allDataLayoutOptions = {DataLayoutOption::aos, DataLayoutOption::soa};

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
   * @param containerSelectorStrategy Strategy for the container selector.
   * @param traversalSelectorStrategy Strategy for the traversal selector.
   * @param tuningInterval Number of timesteps after which the auto-tuner shall reevaluate all selections.
   * @param maxSamples Number of samples that shall be collected for each combination.
   */
  AutoTuner(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff, double verletSkin,
            unsigned int verletRebuildFrequency, std::vector<ContainerOptions> allowedContainerOptions,
            std::vector<TraversalOptions> allowedTraversalOptions, SelectorStrategy containerSelectorStrategy,
            SelectorStrategy traversalSelectorStrategy, unsigned int tuningInterval, unsigned int maxSamples)
      : _tuningInterval(tuningInterval),
        _iterationsSinceTuning(tuningInterval),  // init to max so that tuning happens in first iteration
        _containerSelector(boxMin, boxMax, cutoff, verletSkin, verletRebuildFrequency, allowedContainerOptions,
                           allowedTraversalOptions),
        _allowedTraversalOptions(allowedTraversalOptions),
        _maxSamples(maxSamples),
        _numSamples(maxSamples),
        _containerSelectorStrategy(containerSelectorStrategy),
        _traversalSelectorStrategy(traversalSelectorStrategy) {}

  /**
   * Getter for the optimal container.
   * Also checks if the container was already encountered and if not creates a new traversal selector for it.
   * @return Smartpointer to the optimal container.
   */
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> getContainer() {
    auto container = _containerSelector.getOptimalContainer();
    // if the container is new create a new traversal selector for it
    if (_traversalSelectors.find(container->getContainerType()) == _traversalSelectors.end()) {
      _traversalSelectors.insert(std::make_pair(container->getContainerType(),
                                                container->generateTraversalSelector(_allowedTraversalOptions)));
    }
    return container;
  };

  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @tparam PairwiseFunctor
   * @param f Functor that describes the pair-potential.
   * @param dataLayoutOption if true SoA data structure is used otherwise AoS.
   * @return true if this was a tuning iteration.
   */
  template <class PairwiseFunctor>
  bool iteratePairwise(PairwiseFunctor *f, DataLayoutOption dataLayoutOption);

 private:
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  bool iteratePairwiseTemplateHelper(PairwiseFunctor *f);

  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  bool tune(PairwiseFunctor &pairwiseFunctor);

  unsigned int _tuningInterval, _iterationsSinceTuning;
  ContainerSelector<Particle, ParticleCell> _containerSelector;
  std::vector<TraversalOptions> _allowedTraversalOptions;
  /**
   * One selector per possible container.
   */
  std::map<ContainerOptions, TraversalSelector<ParticleCell>> _traversalSelectors;

  /**
   * How many times each configuration should be tested.
   */
  const size_t _maxSamples;
  /**
   * How many times this configurations has already been tested.
   * Initialize with max value to start tuning at start of simulation.
   */
  size_t _numSamples;

  // TODO: _containerSelectorStrategy currently unused
  SelectorStrategy _containerSelectorStrategy;
  SelectorStrategy _traversalSelectorStrategy;
};

template <class Particle, class ParticleCell>
template <class PairwiseFunctor>
bool AutoTuner<Particle, ParticleCell>::iteratePairwise(PairwiseFunctor *f, DataLayoutOption dataLayoutOption) {
  bool newton3Allowed = f->allowsNewton3();
  bool nonNewton3Allowed = f->allowsNonNewton3();
  bool useNewton3 = false;
  if (newton3Allowed and nonNewton3Allowed) {
    // @todo auto-tune (far off future)
  } else if (not newton3Allowed and not nonNewton3Allowed) {
    utils::ExceptionHandler::exception("Functor::iteratePairwise() Functor neither allows Newton 3 nor not.");
  } else {
    useNewton3 = newton3Allowed;
  }

  bool isTuning = false;

  switch (dataLayoutOption) {
    case DataLayoutOption::aos: {
      if (useNewton3) {
        isTuning = iteratePairwiseTemplateHelper<PairwiseFunctor, false, true>(f);
      } else {
        isTuning = iteratePairwiseTemplateHelper<PairwiseFunctor, false, false>(f);
      }
      break;
    }
    case DataLayoutOption::soa: {
      if (useNewton3) {
        isTuning = iteratePairwiseTemplateHelper<PairwiseFunctor, true, true>(f);
      } else {
        isTuning = iteratePairwiseTemplateHelper<PairwiseFunctor, true, false>(f);
      }
      break;
    }
  }
  return isTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
bool AutoTuner<Particle, ParticleCell>::iteratePairwiseTemplateHelper(PairwiseFunctor *f) {
  bool isTuning = false;
  // large case differentiation for data layout and newton 3
  // check if currently in tuning phase, execute iteration and take time measurement if necessary
  if (_iterationsSinceTuning >= _tuningInterval) {
    if (_numSamples < _maxSamples) {
      isTuning = true;
    } else {
      _numSamples = 0;
      isTuning = tune<PairwiseFunctor, useSoA, useNewton3>(*f);
    }
  }

  auto container = getContainer();
  AutoPasLog(debug, "Using container {}", utils::StringUtils::to_string(container->getContainerType()));

  TraversalSelector<ParticleCell> &traversalSelector = _traversalSelectors[container->getContainerType()];
  auto traversal = traversalSelector.template getOptimalTraversal<PairwiseFunctor, useSoA, useNewton3>(*f);

  // if tuning execute with time measurements
  if (isTuning) {
    auto start = std::chrono::high_resolution_clock::now();
    // @todo remove useNewton3 in iteratePairwise by introducing traversals for DS and VL
    if (useSoA) {
      WithStaticContainerType(container, container->iteratePairwiseSoA(f, traversal.get(), useNewton3););
    } else {
      WithStaticContainerType(container, container->iteratePairwiseAoS(f, traversal.get(), useNewton3););
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    AutoPasLog(debug, "IteratePairwiese took {} nanoseconds", runtime);
    _containerSelector.addTimeMeasurement(container->getContainerType(), runtime);
    traversalSelector.addTimeMeasurement(*f, traversal->getTraversalType(), runtime);
  } else {
    if (useSoA) {
      WithStaticContainerType(container, container->iteratePairwiseSoA(f, traversal.get(), useNewton3););
    } else {
      WithStaticContainerType(container, container->iteratePairwiseAoS(f, traversal.get(), useNewton3););
    }
  }
  ++_iterationsSinceTuning;
  ++_numSamples;
  return isTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
bool AutoTuner<Particle, ParticleCell>::tune(PairwiseFunctor &pairwiseFunctor) {
  auto traversal = _traversalSelectors[getContainer()->getContainerType()]
                       .template selectNextTraversal<PairwiseFunctor, useSoA, useNewton3>(pairwiseFunctor);
  // if there is no next traversal take next container
  if (traversal) {
    return true;
  } else {
    auto container = _containerSelector.selectNextContainer();

    // if there is no next container everything is tested and the optimum can be chosen
    if (container) {
      _traversalSelectors[getContainer()->getContainerType()]
          .template selectNextTraversal<PairwiseFunctor, useSoA, useNewton3>(pairwiseFunctor);
      return true;
    } else {
      _containerSelector.selectOptimalContainer();
      _traversalSelectors[getContainer()->getContainerType()]
          .template selectOptimalTraversal<PairwiseFunctor, useSoA, useNewton3>(_traversalSelectorStrategy,
                                                                                pairwiseFunctor);
      _iterationsSinceTuning = 0;
      return false;
    }
  }
}
}  // namespace autopas
