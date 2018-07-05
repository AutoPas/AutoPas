/**
 * @file ContainerSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <containers/DirectSum.h>
#include <containers/LinkedCells.h>
#include <containers/ParticleContainer.h>
#include <containers/VerletLists.h>
#include <array>
#include <vector>
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
static std::vector<ContainerOptions> allContainerOptions = {ContainerOptions::directSum, ContainerOptions::linkedCells,
                                                            ContainerOptions::verletLists};

/**
 * Selector for a particle container.
 *
 * The class is given a list of allowed container and traversal options to choose from.
 * This class selects the optimal container and delegates the choice of the optimal traversal down to this container.
 *
 * @tparam Particle
 * @tparam ParticleCell
 */
template <class Particle, class ParticleCell>
class ContainerSelector {
 public:
  /**
   * Constructor for the ContainerSelecor class.
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff  Cutoff radius to be used in this container.
   * @param verletSkin Length added to the cutoff for the verlet lists' skin.
   * @param verletRebuildFrequency Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   * @param allowedContainerOptions Vector of container types the selector can choose from.
   * @param allowedTraversalOptions Vector of traversals the selector can choose from.
   */
  ContainerSelector(std::array<double, 3> &boxMin, std::array<double, 3> &boxMax, double cutoff, double verletSkin,
                    unsigned int verletRebuildFrequency, std::vector<ContainerOptions> allowedContainerOptions,
                    std::vector<TraversalOptions> allowedTraversalOptions)
      : _boxMin(boxMin),
        _boxMax(boxMax),
        _cutoff(cutoff),
        _verletSkin(verletSkin),
        _verletRebuildFrequency(verletRebuildFrequency),
        _allowedContainerOptions(std::move(allowedContainerOptions)),
        _allowedTraversalOptions(std::move(allowedTraversalOptions)),
        _optimalContainer(nullptr) {}

  /**
   * Getter for the optimal container.
   * @return Smartpointer to the optimal container.
   */
  std::shared_ptr<ParticleContainer<Particle, ParticleCell>> getOptimalContainer();

  /**
   * Evaluates the optimal container option.
   */
  void tune();

 private:
  std::vector<std::unique_ptr<ParticleContainer<Particle, ParticleCell>>> generateContainers();
  void chooseOptimalContainer(std::vector<std::unique_ptr<ParticleContainer<Particle, ParticleCell>>> containers);

  std::array<double, 3> _boxMin, _boxMax;
  double _cutoff;
  double _verletSkin;
  unsigned int _verletRebuildFrequency;
  std::vector<ContainerOptions> _allowedContainerOptions;
  std::vector<TraversalOptions> _allowedTraversalOptions;
  std::shared_ptr<ParticleContainer<Particle, ParticleCell>> _optimalContainer;
};

template <class Particle, class ParticleCell>
std::vector<std::unique_ptr<ParticleContainer<Particle, ParticleCell>>>
ContainerSelector<Particle, ParticleCell>::generateContainers() {
  std::vector<std::unique_ptr<ParticleContainer<Particle, ParticleCell>>> containers;

  for (auto &option : _allowedContainerOptions) {
    switch (option) {
      case directSum: {
        containers.push_back(std::make_unique<DirectSum<Particle, ParticleCell>>(_boxMin, _boxMax, _cutoff));
        break;
      }
      case linkedCells: {
        containers.push_back(
            std::make_unique<LinkedCells<Particle, ParticleCell>>(_boxMin, _boxMax, _cutoff, _allowedTraversalOptions));
        break;
      }
      case verletLists: {
        // TODO determine verletSkin and verletRebuildFrequency via tuning
        containers.push_back(
            std::make_unique<VerletLists<Particle>>(_boxMin, _boxMax, _cutoff, _verletSkin, _verletRebuildFrequency));
        break;
      }
      default: { AutoPasLogger->warn("Container type {} is not a known type!", option); }
    }

    // copy particles so they do not get lost when container is switched
    // TODO: optimize this such that we do not save the whole domain x times
    if (_optimalContainer != nullptr) {
      for (auto particleIter = _optimalContainer->begin(IteratorBehavior::ownedOnly); particleIter.isValid();
           ++particleIter) {
        containers[containers.size() - 1]->addParticle(*particleIter);
      }
      for (auto particleIter = _optimalContainer->begin(IteratorBehavior::haloOnly); particleIter.isValid();
           ++particleIter) {
        containers[containers.size() - 1]->addHaloParticle(*particleIter);
      }
    }
  }

  assert(containers.size() > 0);
  return containers;
}

template <class Particle, class ParticleCell>
void ContainerSelector<Particle, ParticleCell>::chooseOptimalContainer(
    std::vector<std::unique_ptr<ParticleContainer<Particle, ParticleCell>>> containers) {
  // TODO: Autotuning goes here
  _optimalContainer = std::move(containers.front());
}

template <class Particle, class ParticleCell>
std::shared_ptr<ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::getOptimalContainer() {
  if (_optimalContainer == nullptr) tune();
  return _optimalContainer;
}
template <class Particle, class ParticleCell>
void ContainerSelector<Particle, ParticleCell>::tune() {
  chooseOptimalContainer(generateContainers());
}
}  // namespace autopas