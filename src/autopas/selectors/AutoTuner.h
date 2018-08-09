/**
 * @file AutoTuner.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>
#include "autopas/AutoPas.h"
#include "autopas/autopasIncludes.h"
#include "autopas/pairwiseFunctors/Functor.h"
#include "autopas/selectors/ContainerSelector.h"

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
                           allowedTraversalOptions) {}

  /**
   * Getter for the optimal container.
   * @return Smartpointer to the optimal container.
   */
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> getContainer() {
    return _containerSelector.getOptimalContainer();
  };

  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @tparam ParticleFunctor
   * @param f Functor that describes the pair-potential
   * @param dataLayoutOption if true SoA data structure is used otherwise AoS
   */
  template <class ParticleFunctor>
  void iteratePairwise(ParticleFunctor *f, DataLayoutOption dataLayoutOption);

 private:
  template <class PairwiseFunctor, bool useSoA, bool useNewton3>
  void tune(PairwiseFunctor &pairwiseFunctor);

  unsigned int _tuningInterval, _iterationsSinceTuning;
  ContainerSelector<Particle, ParticleCell> _containerSelector;
  bool _isTuningTraversals = false;
  bool _isTuningContainers = false;
};
template <class Particle, class ParticleCell>
template <class ParticleFunctor>
void AutoTuner<Particle, ParticleCell>::iteratePairwise(ParticleFunctor *f, DataLayoutOption dataLayoutOption) {
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
  switch (dataLayoutOption) {
    case autopas::soa: {
      if (useNewton3) {
        if (_iterationsSinceTuning >= _tuningInterval) {
          tune<ParticleFunctor, true, true>(*f);
        }
      } else {
        if (_iterationsSinceTuning >= _tuningInterval) {
          tune<ParticleFunctor, true, false>(*f);
        }
      }
      auto container = _containerSelector.getOptimalContainer();
      WithStaticContainerType(container, container->iteratePairwiseSoA(f, useNewton3););
      break;
    }
    case autopas::aos: {
      if (useNewton3) {
        if (_iterationsSinceTuning >= _tuningInterval) {
          tune<ParticleFunctor, false, true>(*f);
        }
      } else {
        if (_iterationsSinceTuning >= _tuningInterval) {
          tune<ParticleFunctor, false, false>(*f);
        }
      }
      auto container = _containerSelector.getOptimalContainer();
      WithStaticContainerType(container, container->iteratePairwiseAoS(f, useNewton3););
      break;
    }
  }
  ++_iterationsSinceTuning;
}

template <class Particle, class ParticleCell>
template <class PairwiseFunctor, bool useSoA, bool useNewton3>
void AutoTuner<Particle, ParticleCell>::tune(PairwiseFunctor &pairwiseFunctor) {
  if (not _isTuningTraversals) {
    _isTuningContainers = _containerSelector.tune();
  };

  _isTuningTraversals =
      _containerSelector.getOptimalContainer()->template tuneTraversal<PairwiseFunctor, useSoA, useNewton3>(
          pairwiseFunctor);

  // reset counter only if all tuning phases are done
  if (not _isTuningTraversals and not _isTuningContainers) _iterationsSinceTuning = 0;
}
}  // namespace autopas