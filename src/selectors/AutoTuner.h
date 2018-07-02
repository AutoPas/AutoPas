/**
 * @file AutoTuner.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <AutoPas.h>
#include <autopasIncludes.h>
#include <pairwiseFunctors/Functor.h>
#include <array>
#include "ContainerSelector.h"
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
   * Constructor for the AutoTuner
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff  Cutoff radius to be used in this container.
   * @param allowedContainerOptions Vector of container types AutoPas can choose from.
   * @param allowedTraversalOptions Vector of traversals AutoPas can choose from.
   * @param tuningInterval Number of timesteps after which the auto-tuner shall reevaluate all selections.
   */
  AutoTuner(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff,
            std::vector<ContainerOptions> allowedContainerOptions,
            std::vector<TraversalOptions> allowedTraversalOptions, unsigned int tuningInterval)
      : _tuningInterval(tuningInterval),
        _iterationsSinceTuning(tuningInterval),  // init to max so that tuning happens in first iteration
        _containerSelector(boxMin, boxMax, cutoff, allowedContainerOptions, allowedTraversalOptions) {}

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
   * @param useSoA if true SoA data structure is used otherwise AoS
   */
  template <class ParticleFunctor>
  void iteratePairwise(ParticleFunctor *f, bool useSoA);

 private:
  template <class CellFunctor>
  void tune(CellFunctor &cellFunctor);

  unsigned int _tuningInterval, _iterationsSinceTuning;
  ContainerSelector<Particle, ParticleCell> _containerSelector;
};
template <class Particle, class ParticleCell>
template <class ParticleFunctor>
void AutoTuner<Particle, ParticleCell>::iteratePairwise(ParticleFunctor *f, bool useSoA) {
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

  auto container = _containerSelector.getOptimalContainer();

  /// @todo: WANT one single iteratePairwise(CellFunctor) for containers
  /// @todo: CellFunctor for iteration should be build here using selectors for SoA and N3
  if (useSoA) {
    if (useNewton3) {
      //      CellFunctor<Particle, ParticleCell, ParticleFunctor, true, true> cellFunctor(f);
      //      if (_iterationsSinceTuning == _tuningInterval)
      //        tune(cellFunctor);
      //      _containerSelector.getOptimalContainer()->iteratePairwiseSoA(cellFunctor);
      WithStaticContainerType(container, container->iteratePairwiseSoA(f, useNewton3););
    } else {
      //      CellFunctor<Particle, ParticleCell, ParticleFunctor, true, false> cellFunctor(f);
      //      if (_iterationsSinceTuning == _tuningInterval)
      //        tune(cellFunctor);
      //      _containerSelector.getOptimalContainer()->iteratePairwiseSoA(cellFunctor);
      WithStaticContainerType(container, container->iteratePairwiseSoA(f, useNewton3););
    }
  } else {
    if (useNewton3) {
      //      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, true> cellFunctor(f);
      //      if (_iterationsSinceTuning == _tuningInterval)
      //        tune(cellFunctor);
      //      _containerSelector.getOptimalContainer()->iteratePairwiseAoS(cellFunctor);
      WithStaticContainerType(container, container->iteratePairwiseAoS(f, useNewton3););
    } else {
      //      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, false> cellFunctor(f);
      //      if (_iterationsSinceTuning == _tuningInterval)
      //        tune(cellFunctor);
      //      _containerSelector.getOptimalContainer()->iteratePairwiseAoS(cellFunctor);
      WithStaticContainerType(container, container->iteratePairwiseAoS(f, useNewton3););
    }
  }

  ++_iterationsSinceTuning;
}
template <class Particle, class ParticleCell>
template <class CellFunctor>
void AutoTuner<Particle, ParticleCell>::tune(CellFunctor &cellFunctor) {
  _containerSelector.tune();
  //  _containerSelector.getOptimalContainer()->tuneTraversal(cellFunctor);
  _iterationsSinceTuning = 0;
}
}  // namespace autopas