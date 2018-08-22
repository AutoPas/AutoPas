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
      if (container->getContainerType() == ContainerOptions::linkedCells) {
        _traversalSelectors.insert(std::make_pair(container->getContainerType(),
                                                  container->generateTraversalSelector(_allowedTraversalOptions)));
        // @todo think about how to handle traversal selectors for other containers (dedicated selector types?)
        //      } else {
      }
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
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  bool tune(PairwiseFunctor &pairwiseFunctor);

  unsigned int _tuningInterval, _iterationsSinceTuning;
  ContainerSelector<Particle, ParticleCell> _containerSelector;
  bool _isTuningTraversals = false;
  bool _isTuningContainers = false;
  std::vector<TraversalOptions> _allowedTraversalOptions;
  // one selector per possible container
  std::map<ContainerOptions, TraversalSelector<ParticleCell>> _traversalSelectors;
};

template <class Particle, class ParticleCell>
template <class ParticleFunctor>
bool AutoTuner<Particle, ParticleCell>::iteratePairwise(ParticleFunctor *f, DataLayoutOption dataLayoutOption) {
  bool newton3Allowed = f->allowsNewton3();
  bool nonNewton3Allowed = f->allowsNonNewton3();
  bool useNewton3 = false;
  if (newton3Allowed and nonNewton3Allowed) {
    /// @todo auto-tune (far off future)
  } else if (not newton3Allowed and not nonNewton3Allowed) {
    /// @todo throw exception
  } else {
    useNewton3 = newton3Allowed;
  }

  /// @todo: WANT one single iteratePairwise(CellFunctor) for containers
  /// @todo: CellFunctor for iteration should be build here using selectors for SoA and N3
  bool isTuning = false;
  switch (dataLayoutOption) {
    case autopas::soa: {
      if (useNewton3) {
        if (_iterationsSinceTuning >= _tuningInterval) {
          isTuning = tune<ParticleFunctor, true, true>(*f);
        }
      } else {
        if (_iterationsSinceTuning >= _tuningInterval) {
          isTuning = tune<ParticleFunctor, true, false>(*f);
        }
      }
      auto container = getContainer();
      AutoPasLogger->debug("AutoTuner: Using container {}", container->getContainerType());

      std::unique_ptr<CellPairTraversal<ParticleCell>> traversal;
      if (container->getContainerType() == ContainerOptions::linkedCells) {
        TraversalSelector<ParticleCell> &traversalSelector = _traversalSelectors[container->getContainerType()];
        if (useNewton3) {
          traversal = traversalSelector.template getOptimalTraversal<ParticleFunctor, true, true>(*f);
        } else {
          traversal = traversalSelector.template getOptimalTraversal<ParticleFunctor, true, false>(*f);
        }

        auto start = std::chrono::high_resolution_clock::now();
        // @todo remove useNewton3 in iteratePairwise by introducing traversals for DS and VL
        WithStaticContainerType(container, container->iteratePairwiseSoA(f, traversal.get(), useNewton3););
        auto stop = std::chrono::high_resolution_clock::now();
        auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
        _containerSelector.addTimeMeasurement(container->getContainerType(), runtime);
        traversalSelector.addTimeMeasurement(*f, traversal->getTraversalType(), runtime);
      } else {
        // dummy traversal
        traversal = std::unique_ptr<C08Traversal<ParticleCell, ParticleFunctor, false, false>>(
            new C08Traversal<ParticleCell, ParticleFunctor, false, false>({0, 0, 0}, f));
        auto start = std::chrono::high_resolution_clock::now();
        WithStaticContainerType(container, container->iteratePairwiseSoA(f, traversal.get(), useNewton3););
        auto stop = std::chrono::high_resolution_clock::now();
        auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
        _containerSelector.addTimeMeasurement(container->getContainerType(), runtime);
      }
      break;
    }
    case autopas::aos: {
      if (useNewton3) {
        if (_iterationsSinceTuning >= _tuningInterval) {
          isTuning = tune<ParticleFunctor, false, true>(*f);
        }
      } else {
        if (_iterationsSinceTuning >= _tuningInterval) {
          isTuning = tune<ParticleFunctor, false, false>(*f);
        }
      }
      auto container = getContainer();
      AutoPasLogger->debug("AutoTuner: Using container {}", container->getContainerType());

      std::unique_ptr<CellPairTraversal<ParticleCell>> traversal;
      if (container->getContainerType() == ContainerOptions::linkedCells) {
        TraversalSelector<ParticleCell> &traversalSelector = _traversalSelectors[container->getContainerType()];
        if (useNewton3) {
          traversal = traversalSelector.template getOptimalTraversal<ParticleFunctor, false, true>(*f);
        } else {
          traversal = traversalSelector.template getOptimalTraversal<ParticleFunctor, false, false>(*f);
        }

        auto start = std::chrono::high_resolution_clock::now();
        WithStaticContainerType(container, container->iteratePairwiseAoS(f, traversal.get(), useNewton3););
        auto stop = std::chrono::high_resolution_clock::now();
        auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
        _containerSelector.addTimeMeasurement(container->getContainerType(), runtime);
        traversalSelector.addTimeMeasurement(*f, traversal->getTraversalType(), runtime);
      } else {
        // dummy traversal
        traversal = std::unique_ptr<C08Traversal<ParticleCell, ParticleFunctor, false, false>>(
            new C08Traversal<ParticleCell, ParticleFunctor, false, false>({0, 0, 0}, f));
        auto start = std::chrono::high_resolution_clock::now();
        WithStaticContainerType(container, container->iteratePairwiseAoS(f, traversal.get(), useNewton3););
        auto stop = std::chrono::high_resolution_clock::now();
        auto runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count();
        _containerSelector.addTimeMeasurement(container->getContainerType(), runtime);
      }
      break;
    }
  }
  ++_iterationsSinceTuning;
  return isTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
bool AutoTuner<Particle, ParticleCell>::tune(PairwiseFunctor &pairwiseFunctor) {
  if (not _isTuningTraversals) {
    _isTuningContainers = _containerSelector.tune();
  };

  // only tune traversals for containers that actually have a tuner
  if (_traversalSelectors.find(getContainer()->getContainerType()) != _traversalSelectors.end()) {
    _isTuningTraversals =
        _traversalSelectors[getContainer()->getContainerType()].template tune<PairwiseFunctor, useSoA, useNewton3>(
            pairwiseFunctor);
  } else {
    _isTuningTraversals = false;
  }

  // reset counter only if all tuning phases are done
  if (not _isTuningTraversals and not _isTuningContainers) {
    _iterationsSinceTuning = 0;
    return false;
  }
  return true;
}
}  // namespace autopas