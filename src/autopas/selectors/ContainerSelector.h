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
#include "autopas/options/ContainerOptions.h"
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
   * @param blackBoxMode Indicates whether the [blackbox mode](@ref md_BlackBoxMode) shall be used.
   */
  ContainerSelector(std::array<double, 3> &boxMin, std::array<double, 3> &boxMax, double cutoff, double verletSkin,
                    unsigned int verletRebuildFrequency, std::vector<ContainerOptions> allowedContainerOptions,
                    std::vector<TraversalOptions> allowedTraversalOptions, bool blackBoxMode)
      : _boxMin(boxMin),
        _boxMax(boxMax),
        _cutoff(cutoff),
        _verletSkin(verletSkin),
        _verletRebuildFrequency(verletRebuildFrequency),
        _allowedContainerOptions(std::move(allowedContainerOptions)),
        _allowedTraversalOptions(std::move(allowedTraversalOptions)),
        _optimalContainer(nullptr),
        _blackBoxMode(blackBoxMode) {
    // @FIXME This is a workaround because this container does not yet use traversals like it should
    _allowedTraversalOptions.push_back(TraversalOptions::dummyTraversal);
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
  void addTimeMeasurement(ContainerOptions container, long time);

 private:
  /**
   * Container factory that also copies all particles to the new container
   * @param containerChoice container to generate
   * @return smartpointer to new container
   */
  std::unique_ptr<autopas::ParticleContainer<Particle, ParticleCell>> generateContainer(
      ContainerOptions containerChoice);

  std::array<double, 3> _boxMin, _boxMax;
  double _cutoff;
  double _verletSkin;
  unsigned int _verletRebuildFrequency;
  std::vector<ContainerOptions> _allowedContainerOptions;
  std::vector<TraversalOptions> _allowedTraversalOptions;
  std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>> _optimalContainer;
  std::vector<std::pair<ContainerOptions, long>> _containerTimes;
  bool _blackBoxMode;
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
      // @todo determine verletSkin and verletRebuildFrequency via tuning
      container =
          std::make_unique<VerletLists<Particle>>(_boxMin, _boxMax, _cutoff, _verletSkin, _verletRebuildFrequency,
                                                  VerletLists<Particle>::BuildVerletListType::VerletSoA, _blackBoxMode);
      break;
    }
    case verletListsCells: {
      // @todo determine verletSkin and verletRebuildFrequency via tuning
      container = std::make_unique<VerletListsCells<Particle>>(_boxMin, _boxMax, _cutoff, TraversalOptions::c08,
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
std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::selectNextContainer() {
  ContainerOptions nextContainerType;

  // if none is chosen yet create a new
  if (_optimalContainer == nullptr) {
    nextContainerType = _allowedContainerOptions.begin().operator*();
  } else {
    auto containerTypeIter = std::find(_allowedContainerOptions.begin(), _allowedContainerOptions.end(),
                                       _optimalContainer->getContainerType());
    ++containerTypeIter;

    if (containerTypeIter >= _allowedContainerOptions.end()) {
      return std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>>(nullptr);
    }

    nextContainerType = *containerTypeIter;
  }
  _optimalContainer = std::move(generateContainer(nextContainerType));
  AutoPasLog(debug, "Testing Container {}", utils::StringUtils::to_string(nextContainerType));

  return _optimalContainer;
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
  ContainerOptions optimalContainerOption = ContainerOptions::directSum;
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
  if (_optimalContainer == nullptr || _optimalContainer->getContainerType() != optimalContainerOption) {
    _optimalContainer = std::move(generateContainer(optimalContainerOption));
  }
  AutoPasLog(debug, "Selected container {}", utils::StringUtils::to_string(_optimalContainer->getContainerType()));
  return _optimalContainer;
}

template <class Particle, class ParticleCell>
std::shared_ptr<autopas::ParticleContainer<Particle, ParticleCell>>
ContainerSelector<Particle, ParticleCell>::getOptimalContainer() {
  if (_optimalContainer == nullptr) selectNextContainer();
  //    utils::ExceptionHandler::exception("ContainerSelector::getOptimalContainer(): No Container selected yet!");
  return _optimalContainer;
}

template <class Particle, class ParticleCell>
void ContainerSelector<Particle, ParticleCell>::addTimeMeasurement(ContainerOptions container, long time) {
  _containerTimes.push_back(std::make_pair(container, time));
}
}  // namespace autopas
