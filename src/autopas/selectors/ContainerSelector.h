/**
 * @file ContainerSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>
#include <vector>
#include "autopas/containers/DirectSum.h"
#include "autopas/containers/LinkedCells.h"
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/VerletLists.h"

namespace autopas {

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
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> getOptimalContainer();

  /**
   * Evaluates the optimal container option.
   * @return true if still in tuning phase
   */
  bool tune();

  /**
   * Save the runtime of a given container.
   * If a better (from a different traversal) already exists nothing happens.
   * If a worse (from a different traversal) already exists it is overwritten.
   * @param container
   * @param time
   */
  void addTimeMeasurement(ContainerOptions container, long time);

 private:
  /**
   * Container factory that also copies all particles to the new container
   * @param containerChoice container to generate
   * @return smartpointer to new container
   */
  std::unique_ptr<autopas::ParticleContainer<Particle, ParticleCell>> generateContainer(
      ContainerOptions containerChoice);

  /**
   * Chooses the optimal container from a given list.
   * @return true if still in tuning phase
   */
  bool chooseOptimalContainer();

  std::array<double, 3> _boxMin, _boxMax;
  double _cutoff;
  double _verletSkin;
  unsigned int _verletRebuildFrequency;
  std::vector<ContainerOptions> _allowedContainerOptions;
  std::vector<TraversalOptions> _allowedTraversalOptions;
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> _optimalContainer;
  std::vector<std::pair<ContainerOptions, long>> _containerTimes;
  bool _currentlyTuning = false;
};

template <class Particle, class ParticleCell>
std::unique_ptr<autopas::ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::generateContainer(ContainerOptions containerChoice) {
  std::unique_ptr<autopas::ParticleContainer<Particle, ParticleCell>> container;

  switch (containerChoice) {
    case directSum: {
      container = std::make_unique<DirectSum<Particle, ParticleCell>>(_boxMin, _boxMax, _cutoff);
      break;
    }
    case linkedCells: {
      container = std::make_unique<LinkedCells<Particle, ParticleCell>>(_boxMin, _boxMax, _cutoff);
      break;
    }
    case verletLists: {
      // TODO determine verletSkin and verletRebuildFrequency via tuning
      container =
          std::make_unique<VerletLists<Particle>>(_boxMin, _boxMax, _cutoff, _verletSkin, _verletRebuildFrequency);
      break;
    }
    default: { AutoPasLogger->warn("Container type {} is not a known type!", containerChoice); }
  }

  // copy particles so they do not get lost when container is switched
  if (_optimalContainer != nullptr) {
    for (auto particleIter = _optimalContainer->begin(IteratorBehavior::ownedOnly); particleIter.isValid();
         ++particleIter) {
      container->addParticle(*particleIter);
    }
    for (auto particleIter = _optimalContainer->begin(IteratorBehavior::haloOnly); particleIter.isValid();
         ++particleIter) {
      container->addHaloParticle(*particleIter);
    }
  }

  return container;
}

template <class Particle, class ParticleCell>
bool ContainerSelector<Particle, ParticleCell>::chooseOptimalContainer() {
  size_t bestContainerID = 0;

  // Test all options to find the fastest
  // If there is no container chosen yet or no measurements were made by now choose the first
  if (_optimalContainer == nullptr || _containerTimes.size() == 0) {
    _optimalContainer = std::move(generateContainer(_allowedContainerOptions.front()));
  } else if (_currentlyTuning) {
    // if we are in tuning state just select next container
    bestContainerID = std::find(_allowedContainerOptions.begin(), _allowedContainerOptions.end(),
                                _optimalContainer->getContainerType()) -
                      _allowedContainerOptions.begin();
    ++bestContainerID;
    // if the last possible traversal has already been tested choose fastest one and reset timings
    if (bestContainerID >= _allowedContainerOptions.size()) {
      _currentlyTuning = false;
      ContainerOptions fastestContainer;
      long fastestTime = std::numeric_limits<long>::max();
      AutoPasLogger->debug("ContainerSelector.tune(): ContainerSelector: Collected containers:");
      for (auto &&c : _containerTimes) {
        AutoPasLogger->debug("ContainerSelector.tune(): Container {} took {} nanoseconds:", c.first, c.second);
        if (c.second < fastestTime) {
          fastestContainer = c.first;
          fastestTime = c.second;
        }
      }
      // sanity check
      if (fastestTime == std::numeric_limits<long>::max()) {
        utils::ExceptionHandler::exception("ContainerSelector: nothing was faster than max long oO");
      }

      // find id of fastest container in passed container list
      bestContainerID = std::find(_allowedContainerOptions.begin(), _allowedContainerOptions.end(), fastestContainer) -
                        _allowedContainerOptions.begin();
      _containerTimes.clear();
    }

    _optimalContainer = std::move(generateContainer(_allowedContainerOptions[bestContainerID]));
  }
  AutoPasLogger->debug("ContainerSelector.tune(): Selected container {}", _optimalContainer->getContainerType());
  return _currentlyTuning;
}

template <class Particle, class ParticleCell>
std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::getOptimalContainer() {
  if (_optimalContainer == nullptr) tune();
  return _optimalContainer;
}
template <class Particle, class ParticleCell>
bool ContainerSelector<Particle, ParticleCell>::tune() {
  _currentlyTuning = true;
  return chooseOptimalContainer();
}

template <class Particle, class ParticleCell>
void ContainerSelector<Particle, ParticleCell>::addTimeMeasurement(ContainerOptions container, long time) {
  bool found = false;
  for (auto &&p : _containerTimes) {
    if (p.first == container && p.second > time) {
      found = true;
      p.second = time;
    }
  }
  if (not found) {
    _containerTimes.push_back(std::make_pair(container, time));
  }
}
}  // namespace autopas