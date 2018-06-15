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
  /**
   * Constructor for the ContainerSelecor class.
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff  Cutoff radius to be used in this container.
   * @param allowedContainers Vector of container types the selector can choose from.
   * @param allowedTraversals Vector of traversals the selector can choose from.
   */
  ContainerSelector(std::array<double, 3> &boxMin,
                    std::array<double, 3> &boxMax,
                    double cutoff,
                    std::vector<ContainerOptions> &allowedContainerOptions,
                    std::vector<TraversalOptions> &allowedTraversalOptions
  ) : _boxMin(boxMin),
      _boxMax(boxMax),
      _cutoff(cutoff),
      _allowedContainerOptions(allowedContainerOptions),
      _allowedTraversalOptions(allowedTraversalOptions),
      _optimalContainer(nullptr) {
  }

  /**
   * Getter for the optimal container.
   * @return Smartpointer to the optimal container.
   */
  std::shared_ptr<ParticleContainer<Particle, ParticleCell> > getOptimalContainer();

  /**
   * Evaluates the optimal container option.
   */
  void tune();

 private:

  std::vector<std::unique_ptr<ParticleContainer<Particle, ParticleCell>>> generateContainers();
  void chooseOptimalContainer(std::vector<std::unique_ptr<ParticleContainer<Particle, ParticleCell> >> containers);

  std::array<double, 3> _boxMin, _boxMax;
  double _cutoff;
  std::vector<ContainerOptions> _allowedContainerOptions;
  std::vector<TraversalOptions> _allowedTraversalOptions;
  std::shared_ptr<ParticleContainer<Particle, ParticleCell> > _optimalContainer;
};

template<class Particle, class ParticleCell>
std::vector<std::unique_ptr<ParticleContainer<Particle, ParticleCell>>> ContainerSelector<Particle,
                                                                                          ParticleCell>::generateContainers() {

  std::vector<std::unique_ptr<ParticleContainer<Particle, ParticleCell> >> containers;

  for (auto &option: _allowedContainerOptions) {
    switch (option) {
      case directSum : {
        containers.push_back(std::make_unique<DirectSum<Particle, ParticleCell>>(_boxMin,
                                                                                 _boxMax,
                                                                                 _cutoff));
        break;
      }
      case linkedCells : {
        containers.push_back(std::make_unique<LinkedCells<Particle, ParticleCell>>(_boxMin,
                                                                                   _boxMax,
                                                                                   _cutoff,
                                                                                   _allowedTraversalOptions));
        break;
      }
      case verletLists : {
        containers.push_back(std::make_unique<DirectSum<Particle,
                                                        ParticleCell>>(_boxMin,
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
void ContainerSelector<Particle, ParticleCell>::chooseOptimalContainer(std::vector<std::unique_ptr<ParticleContainer<
    Particle,
    ParticleCell> >> containers) {
  // TODO: Autotuning goes here
  _optimalContainer = std::move(containers.front());
}

template<class Particle, class ParticleCell>
std::shared_ptr<ParticleContainer<Particle, ParticleCell>> ContainerSelector<Particle,
                                                                             ParticleCell>::getOptimalContainer() {
  if (_optimalContainer == nullptr)
    tune();
  return _optimalContainer;
}
template<class Particle, class ParticleCell>
void ContainerSelector<Particle, ParticleCell>::tune() {
  chooseOptimalContainer(generateContainers());
}
}