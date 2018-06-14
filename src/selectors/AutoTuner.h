/**
 * AutoTuner.h
 *
 *  Created on: 6/11/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include <array>
#include <pairwiseFunctors/Functor.h>
#include <AutoPas.h>
#include "ContainerSelector.h"
namespace autopas {

template<class Particle, class ParticleCell>
class AutoTuner {
 public:
  AutoTuner(std::array<double, 3> boxMin,
            std::array<double, 3> boxMax,
            double cutoff,
            unsigned int tuningInterval,
            std::vector<ContainerOptions> allowedContainerOptions,
            std::vector<TraversalOptions> allowedTraversalOptions
  ) : tuningInterval(tuningInterval),
      iterationsSinceTuning(tuningInterval),    // init to max so that tuning happens in first iteration
      _containerSelector(boxMin,
                         boxMax,
                         cutoff,
                         allowedContainerOptions,
                         allowedTraversalOptions) {
  }

  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> getContainer() {
    return _containerSelector.getOptimalContainer();
  };

  template<class ParticleFunctor>
  void iteratePairwise(ParticleFunctor *f,
                       bool useSoA);

 private:
  template<class CellFunctor>
  void tune(CellFunctor &cellFunctor);

  unsigned int tuningInterval, iterationsSinceTuning;
  ContainerSelector<Particle, ParticleCell> _containerSelector;
};
template<class Particle, class ParticleCell>
template<class ParticleFunctor>
void AutoTuner<Particle, ParticleCell>::iteratePairwise(ParticleFunctor *f,
                                                        bool useSoA) {

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
  if (useSoA) {
    if (useNewton3) {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, true, true> cellFunctor(f);
//      if (iterationsSinceTuning == tuningInterval)
//        tune(cellFunctor);
//      _containerSelector.getOptimalContainer()->iteratePairwiseSoA(cellFunctor);
      _containerSelector.getOptimalContainer()->iteratePairwiseSoA(f, useNewton3);
    } else {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, true, false> cellFunctor(f);
//      if (iterationsSinceTuning == tuningInterval)
//        tune(cellFunctor);
//      _containerSelector.getOptimalContainer()->iteratePairwiseSoA(cellFunctor);
      _containerSelector.getOptimalContainer()->iteratePairwiseSoA(f, useNewton3);
    }
  } else {
    if (useNewton3) {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, true> cellFunctor(f);
//      if (iterationsSinceTuning == tuningInterval)
//        tune(cellFunctor);
//      _containerSelector.getOptimalContainer()->iteratePairwiseAoS(cellFunctor);
      _containerSelector.getOptimalContainer()->iteratePairwiseAoS(f, useNewton3);
    } else {
      CellFunctor<Particle, ParticleCell, ParticleFunctor, false, false> cellFunctor(f);
//      if (iterationsSinceTuning == tuningInterval)
//        tune(cellFunctor);
//      _containerSelector.getOptimalContainer()->iteratePairwiseAoS(cellFunctor);
      _containerSelector.getOptimalContainer()->iteratePairwiseAoS(f, useNewton3);
    }

  }

  ++iterationsSinceTuning;
}
template<class Particle, class ParticleCell>
template<class CellFunctor>
void AutoTuner<Particle, ParticleCell>::tune(CellFunctor &cellFunctor) {
  _containerSelector.tune();
//  _containerSelector.getOptimalContainer()->tuneTraversal(cellFunctor);
  iterationsSinceTuning = 0;
}
}