/**
 * AutoTuner.h
 *
 *  Created on: 6/11/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include <array>
#include <pairwiseFunctors/Functor.h>
#include "ContainerSelector.h"
namespace autopas {

template<class Particle, class ParticleCell>
class AutoTuner {
 public:
  AutoTuner(std::array<double, 3> boxMin,
            std::array<double, 3> boxMax,
            double cutoff,
            unsigned int retuneInterval,
            std::vector<ContainerOptions> allowedContainerOptions,
            std::vector<TraversalOptions> allowedTraversalOptions
  ) : _containerSelector(boxMin,
                         boxMax,
                         cutoff,
                         retuneInterval,
                         allowedContainerOptions,
                         allowedTraversalOptions) {
  }

//  void iteratePairwise(Functor<Particle, ParticleCell> *functor);

  std::unique_ptr<autopas::ParticleContainer<Particle, ParticleCell>> getContainer() {
    _containerSelector.getOptimalContainer();
  };

 private:
  ContainerSelector<Particle, ParticleCell> _containerSelector;

//  std::vector<ContainerOptions> _allowedContainerOptions;
//  std::vector<TraversalOptions> _allowedTraversalOptions;

};
}