/**
 * ContainerSelector.h
 *
 *  Created on: 6/11/18
 *     Aauthor: F. Gratl
 */

#pragma once

#include <array>
#include <vector>
#include <containers/ParticleContainer.h>
#include <containers/DirectSum.h>
#include <containers/LinkedCells.h>
namespace autopas {

/**
 * Possible choices for the particle container type.
 */
enum ContainerOptions {
  directSum = 0,
  linkedCells = 1,
  verletLists = 2,
};

/**
 * Provides a way to iterate over the possible choices of ContainerOption.
 */
static std::vector<ContainerOptions> allContainerOptions = {ContainerOptions::directSum,
                                                            ContainerOptions::linkedCells,
                                                            ContainerOptions::verletLists};

template<class Particle, class ParticleCell>
class ContainerSelector {
 public:
  ContainerSelector(std::array<double, 3> &boxMin,
                    std::array<double, 3> &boxMax,
                    double cutoff,
                    unsigned int retuneInterval,
                    std::vector<ContainerOptions> &allowedContainerOptions,
                    std::vector<TraversalOptions> &allowedTraversalOptions
  ) : _boxMin(boxMin),
      _boxMax(boxMax),
      _cutoff(cutoff),
      _retuneInterval(retuneInterval),
      _retuneCounter(0),
      _allowedContainerOptions(allowedContainerOptions),
      _allowedTraversalOptions(allowedTraversalOptions) {
  }

  ParticleContainer<Particle, ParticleCell> *getOptimalContainer();

 private:

  std::vector<ParticleContainer<Particle, ParticleCell> *> generateContainers();
  void chooseOptimalContainer(std::vector<ParticleContainer<Particle, ParticleCell> *> containers);

  std::array<double, 3> _boxMin, _boxMax;
  double _cutoff;
  unsigned int _retuneInterval, _retuneCounter;
  std::vector<ContainerOptions> _allowedContainerOptions;
  std::vector<TraversalOptions> _allowedTraversalOptions;
  ParticleContainer<Particle, ParticleCell> *_optimalContainer;
};

template<class Particle, class ParticleCell>
std::vector<ParticleContainer<Particle, ParticleCell> *> ContainerSelector<Particle,
                                                                           ParticleCell>::generateContainers() {

  std::vector<ParticleContainer<Particle, ParticleCell> *> containers;

  for (auto &option: _allowedContainerOptions) {
    switch (option) {
      case directSum : {
        containers.push_back(new DirectSum<Particle, ParticleCell>(_boxMin,
                                                                   _boxMax,
                                                                   _cutoff));
        break;
      }
      case linkedCells : {
        containers.push_back(new LinkedCells<Particle, ParticleCell>(_boxMin,
                                                                     _boxMax,
                                                                     _cutoff,
                                                                     0,
                                                                     _allowedTraversalOptions));
        break;
      }
      case verletLists : {
        containers.push_back(new DirectSum<Particle,
                                           ParticleCell>(_boxMin,
                                                         _boxMax,
                                                         _cutoff));
        break;
      }
      default: {
        AutoPasLogger->warn("Container type {} is not a known type!", option);
      }
    }
  }

  assert(containers.size() > 0);
  return containers;
}

template<class Particle, class ParticleCell>
void ContainerSelector<Particle, ParticleCell>::chooseOptimalContainer(std::vector<ParticleContainer<Particle,
                                                                                                     ParticleCell> *> containers) {
  // TODO: Autotuning goes here
  _optimalContainer = containers.front();
}

template<class Particle, class ParticleCell>
ParticleContainer<Particle, ParticleCell> *ContainerSelector<Particle,
                                                             ParticleCell>::getOptimalContainer() {
  if (_retuneCounter == 0) {
    chooseOptimalContainer(generateContainers());
    _retuneCounter = _retuneInterval;
  }

  --_retuneCounter;

  return _optimalContainer;
}
}