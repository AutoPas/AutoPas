/**
 * @file ContainerSelector.h
 * @author F. Gratl
 * @date 11.06.18
 */

#pragma once

#include <array>
#include <vector>
#include "autopas/containers/ParticleContainer.h"
#include "autopas/containers/directSum/DirectSum.h"
#include "autopas/containers/linkedCells/LinkedCells.h"
#include "autopas/containers/verletClusterLists/VerletClusterLists.h"
#include "autopas/containers/verletListsCellBased/verletLists/VerletLists.h"
#include "autopas/containers/verletListsCellBased/verletListsCells/VerletListsCells.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/utils/StringUtils.h"

namespace autopas {

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
                    unsigned int verletRebuildFrequency, std::vector<ContainerOption> allowedContainerOptions,
                    std::vector<TraversalOption> allowedTraversalOptions)
      : _boxMin(boxMin),
        _boxMax(boxMax),
        _cutoff(cutoff),
        _verletSkin(verletSkin),
        _verletRebuildFrequency(verletRebuildFrequency),
        _allowedContainerOptions(std::move(allowedContainerOptions)),
        _allowedTraversalOptions(std::move(allowedTraversalOptions)),
        _currentContainer(nullptr) {
    // @FIXME This is a workaround because this container does not yet use traversals like it should
    _allowedTraversalOptions.push_back(TraversalOption::dummyTraversal);
  }

  /**
   * Selects the next allowed container.
   * @return Smartpointer to the selected container.
   */
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> selectNextContainer();

  /**
   * Selects the optimal container based on saved measurements.
   * @return Smartpointer to the selected container.
   */
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> selectOptimalContainer();

  /**
   * Sets the container to the given option.
   * @param containerOption
   */
  void selectContainer(ContainerOption containerOption);

  /**
   * Getter for the optimal container. If no container is chosen yet the first allowed is selected.
   * @return Smartpointer to the optimal container.
   */
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> getOptimalContainer();

  /**
   * Save the runtime of a given container.
   * If a better (from a different traversal) already exists nothing happens.
   * If a worse (from a different traversal) already exists it is overwritten.
   * @param container
   * @param time
   */
  void addTimeMeasurement(ContainerOption container, long time);

 private:
  /**
   * Container factory that also copies all particles to the new container
   * @param containerChoice container to generate
   * @return smartpointer to new container
   */
  std::unique_ptr<autopas::ParticleContainer<Particle, ParticleCell>> generateContainer(
      ContainerOption containerChoice);

  std::array<double, 3> _boxMin, _boxMax;
  double _cutoff;
  double _verletSkin;
  unsigned int _verletRebuildFrequency;
  std::vector<ContainerOption> _allowedContainerOptions;
  std::vector<TraversalOption> _allowedTraversalOptions;
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> _currentContainer;
  std::vector<std::pair<ContainerOption, long>> _containerTimes;
};

template <class Particle, class ParticleCell>
std::unique_ptr<autopas::ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::generateContainer(ContainerOption containerChoice) {
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
      // @todo determine verletSkin and verletRebuildFrequency via tuning
      container =
          std::make_unique<VerletLists<Particle>>(_boxMin, _boxMax, _cutoff, _verletSkin, _verletRebuildFrequency);
      break;
    }
    case verletListsCells: {
      // @todo determine verletSkin and verletRebuildFrequency via tuning
      container = std::make_unique<VerletListsCells<Particle>>(_boxMin, _boxMax, _cutoff, TraversalOption::c08,
                                                               _verletSkin, _verletRebuildFrequency);
      break;
    }
    case verletClusterLists: {
      // @todo determine verletSkin and verletRebuildFrequency via tuning
      container = std::make_unique<VerletClusterLists<Particle>>(_boxMin, _boxMax, _cutoff, _verletSkin,
                                                                 _verletRebuildFrequency);
      break;
    }
    default: {
      AutoPasLog(warn, "Container type {} is not a known type!", utils::StringUtils::to_string(containerChoice));
    }
  }

  // copy particles so they do not get lost when container is switched
  if (_currentContainer != nullptr) {
    for (auto particleIter = _currentContainer->begin(IteratorBehavior::ownedOnly); particleIter.isValid();
         ++particleIter) {
      container->addParticle(*particleIter);
    }
    for (auto particleIter = _currentContainer->begin(IteratorBehavior::haloOnly); particleIter.isValid();
         ++particleIter) {
      container->addHaloParticle(*particleIter);
    }
  }

  return container;
}

template <class Particle, class ParticleCell>
std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::selectNextContainer() {
  ContainerOption nextContainerType;

  // if none is chosen yet create a new
  if (_currentContainer == nullptr) {
    nextContainerType = _allowedContainerOptions.begin().operator*();
  } else {
    auto containerTypeIter = std::find(_allowedContainerOptions.begin(), _allowedContainerOptions.end(),
                                       _currentContainer->getContainerType());
    ++containerTypeIter;

    if (containerTypeIter >= _allowedContainerOptions.end()) {
      return std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>>(nullptr);
    }

    nextContainerType = *containerTypeIter;
  }
  _currentContainer = std::move(generateContainer(nextContainerType));
  AutoPasLog(debug, "Testing Container {}", utils::StringUtils::to_string(nextContainerType));

  return _currentContainer;
}

template <class Particle, class ParticleCell>
std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::selectOptimalContainer() {
  // Time measure strategy
  if (_containerTimes.empty()) {
    utils::ExceptionHandler::exception("ContainerSelector: Trying to determine fastest container before measuring!");
  }

  // choose the fastest container and reset timings
  // Initialize with something. This will be overridden.
  ContainerOption optimalContainerOption = ContainerOption::directSum;
  long optimalContainerTime = std::numeric_limits<long>::max();
  AutoPasLog(debug, "ContainerSelector: Collected container times:");
  for (auto &&c : _containerTimes) {
    AutoPasLog(debug, "Container {} took {} nanoseconds:", utils::StringUtils::to_string(c.first), c.second);
    if (c.second < optimalContainerTime) {
      optimalContainerOption = c.first;
      optimalContainerTime = c.second;
    }
  }
  // measurements are not needed anymore
  _containerTimes.clear();

  // sanity check
  if (optimalContainerTime == std::numeric_limits<long>::max()) {
    utils::ExceptionHandler::exception("ContainerSelector: Nothing was faster than max long! o_O");
  }

  // only convert when best container is not the current or no container exists yet
  if (_currentContainer == nullptr || _currentContainer->getContainerType() != optimalContainerOption) {
    _currentContainer = std::move(generateContainer(optimalContainerOption));
  }
  AutoPasLog(debug, "Selected container {}", utils::StringUtils::to_string(_currentContainer->getContainerType()));
  return _currentContainer;
}

template <class Particle, class ParticleCell>
std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::getOptimalContainer() {
  if (_currentContainer == nullptr) selectNextContainer();
  return _currentContainer;
}

template <class Particle, class ParticleCell>
void ContainerSelector<Particle, ParticleCell>::addTimeMeasurement(ContainerOption container, long time) {
  _containerTimes.push_back(std::make_pair(container, time));
}
template <class Particle, class ParticleCell>
void ContainerSelector<Particle, ParticleCell>::selectContainer(ContainerOption containerOption) {
  // if we already have this container do nothing.
  if (_currentContainer->getContainerType() != containerOption) {
    _currentContainer = std::move(generateContainer(containerOption));
  }
}
}  // namespace autopas
