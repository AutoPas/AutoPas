/**
 * @file AutoTuner.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>
#include <memory>
#include "autopas/autopasIncludes.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/selectors/ContainerSelector.h"
#include "autopas/selectors/TraversalSelector.h"

namespace autopas {

/**
 * Possible Choices for the particle data layout.
 */
enum DataLayoutOption { aos, soa };

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
   * @param tuningInterval Number of timesteps after which the auto-tuner shall reevaluate all selections.
   */
  AutoTuner(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff, double verletSkin,
            unsigned int verletRebuildFrequency, std::vector<ContainerOptions> allowedContainerOptions,
            std::vector<TraversalOptions> allowedTraversalOptions, unsigned int tuningInterval)
      : _tuningInterval(tuningInterval),
        _iterationsSinceTuning(tuningInterval),  // init to max so that tuning happens in first iteration
        _containerSelector(boxMin, boxMax, cutoff, verletSkin, verletRebuildFrequency, allowedContainerOptions,
                           allowedTraversalOptions),
        _allowedTraversalOptions(allowedTraversalOptions) {}

  /**
   * Getter for the optimal container.
   * Also checks if the container was already encountered and if not creates a new traversal selector for it
   * @return Smartpointer to the optimal container.
   */
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> getContainer() {
    auto container = _containerSelector.getOptimalContainer();
    // if the container is new create a new traversal selector for it
    if (_traversalSelectors.find(container->getContainerType()) == _traversalSelectors.end()) {
      // @todo as soon as all containers work with traversalselectors this if can be removed
      //      if (container->getContainerType() == ContainerOptions::linkedCells) {
      _traversalSelectors.insert(std::make_pair(container->getContainerType(),
                                                container->generateTraversalSelector(_allowedTraversalOptions)));
      // @todo think about how to handle traversal selectors for other containers (dedicated selector types?)
      //      } else {
      //      }
    }
    return container;
  };

  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @tparam ParticleFunctor
   * @param f Functor that describes the pair-potential
   * @param dataLayoutOption if true SoA data structure is used otherwise AoS
   * @return true if this was a tuning iteration
   */
  template <class ParticleFunctor>
  bool iteratePairwise(ParticleFunctor *f, DataLayoutOption dataLayoutOption);

 private:
  template <class ParticleFunctor, bool useSoA, bool useNewton3>
  bool iteratePairwiseTemplateHelper(ParticleFunctor *f);

  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  bool tune(PairwiseFunctor &pairwiseFunctor);

  unsigned int _tuningInterval, _iterationsSinceTuning;
  ContainerSelector<Particle, ParticleCell> _containerSelector;
  bool _isTuningTraversals = false;
  bool _isTuningContainers = false;
  std::vector<TraversalOptions> _allowedTraversalOptions;
  /**
   * One selector per possible container.
   */
  std::map<ContainerOptions, TraversalSelector<ParticleCell>> _traversalSelectors;

  /**
   * How many times one configurations should be tested.
   */
  const size_t _maxSamples = 3;
  /**
   * How many times this configurations has already been tested.
   * Initialize with max value to start tuning at start of simulation.
   */
  size_t _numSamples = _maxSamples;
};

template <class Particle, class ParticleCell>
template <class ParticleFunctor>
bool AutoTuner<Particle, ParticleCell>::iteratePairwise(ParticleFunctor *f, DataLayoutOption dataLayoutOption) {
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
    case DataLayoutOption::aos : {
      if (useNewton3) {
        isTuning = iteratePairwiseTemplateHelper<ParticleFunctor, false, true>(f);
      } else {
        isTuning = iteratePairwiseTemplateHelper<ParticleFunctor, false, false>(f);
      }
      break;
    }
    case DataLayoutOption::soa : {
      if (useNewton3) {
        isTuning = iteratePairwiseTemplateHelper<ParticleFunctor, true, true>(f);
      } else {
        isTuning = iteratePairwiseTemplateHelper<ParticleFunctor, true, false>(f);
      }
      break;
    }
  }
  return isTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
bool AutoTuner<Particle, ParticleCell>::iteratePairwiseTemplateHelper(PairwiseFunctor *f) {
  // @todo: WANT one single iteratePairwise(CellFunctor) for containers
  // @todo: CellFunctor for iteration should be build here using selectors for SoA and N3
  bool isTuning = false;
  // large case differentiation for data layout and newton 3
  // check if currently in tuning phase, execute iteration and take time measurement if necessary
  if (_iterationsSinceTuning >= _tuningInterval) {
    if (_numSamples < _maxSamples) {
      isTuning = true;
    } else {
      _numSamples = 0;
      if (useNewton3) {
        isTuning = tune<PairwiseFunctor, useSoA, useNewton3>(*f);
      }
    }
  }

  auto container = getContainer();
  AutoPasLog(debug, "Using container {}", container->getContainerType());

  TraversalSelector<ParticleCell> &traversalSelector = _traversalSelectors[container->getContainerType()];
  auto traversal = traversalSelector.template getOptimalTraversal<PairwiseFunctor, useSoA, useNewton3>(*f);

  if (isTuning) {
    auto start = std::chrono::high_resolution_clock::now();
    // @todo remove useNewton3 in iteratePairwise by introducing traversals for DS and VL
    WithStaticContainerType(container, container->iteratePairwiseSoA(f, traversal.get(), useNewton3););
    auto stop = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
    AutoPasLog(debug, "IteratePairwiese took {} nanoseconds", runtime);
    _containerSelector.addTimeMeasurement(container->getContainerType(), runtime);
    traversalSelector.addTimeMeasurement(*f, traversal->getTraversalType(), runtime);
  } else {
    WithStaticContainerType(container, container->iteratePairwiseSoA(f, traversal.get(), useNewton3););
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
          .template selectOptimalTraversal<PairwiseFunctor, useSoA, useNewton3>(pairwiseFunctor);
      _iterationsSinceTuning = 0;
      return false;
    }
  }
}
}  // namespace autopas