/**
 * @file AutoPas.h
 * Main include file for the AutoPas library.
 *
 */

#pragma once

#include <iostream>
#include <memory>
#include "autopas/autopasIncludes.h"
#include "autopas/selectors/AutoTuner.h"

namespace autopas {

/**
 * instance counter to help track the number of autopas instances. Needed for correct management of the logger.
 */
static unsigned int _instanceCounter = 0;

/**
 * The AutoPas class is intended to be the main point of Interaction for the
 * user. It puts a layer of abstraction over the container and handles the
 * autotuning.
 * @todo autotuning
 * @tparam Particle Class for particles
 * @tparam ParticleCell Class for the particle cells
 */
template <class Particle, class ParticleCell>
class AutoPas {
 public:
  AutoPas()
      : boxMin{0, 0, 0},
        boxMax{0, 0, 0},
        cutoff(1),
        verletSkin(0.2),
        verletRebuildFrequency(20),
        tuningInterval(5000),
        numSamples(3),
        selectorStrategy(SelectorStrategy::fastestAbs),
        allowedContainers(allContainerOptions),
        allowedTraversals(allTraversalOptions),
        allowedDataLayouts(allDataLayoutOptions),
        allowNewton3(false) {
    // count the number of autopas instances. This is needed to ensure that the autopas
    // logger is not unregistered while other instances are still using it.
    _instanceCounter++;
    if (_instanceCounter == 1) {
      // initialize the Logger
      autopas::Logger::create();
      // The logger is normally only flushed on successful program termination.
      // This line ensures flushing when log messages of level warning or more severe are created.
      autopas::Logger::get()->flush_on(spdlog::level::warn);
    }
  }

  ~AutoPas() {
    _instanceCounter--;
    if (_instanceCounter == 0) {
      // remove the Logger from the registry. Do this only if we have no other autopas instances running.
      autopas::Logger::unregister();
    }
  }

  /**
   * Move assignment operator
   * @param other
   * @return
   */
  AutoPas &operator=(AutoPas &&other) noexcept {
    _autoTuner = std::move(other._autoTuner);
    return *this;
  }

  /**
   * Initialize the auto tuner. This will completely reset the container and remove all containing particles!
   *
   * This function needs to be called before any other functions on the AutoPas object.
   *
   * Changing any of the member options only takes effect when init is called.
   *
   */
  void init() {
    _autoTuner = std::make_unique<autopas::AutoTuner<Particle, ParticleCell>>(
        boxMin, boxMax, cutoff, verletSkin, verletRebuildFrequency, allowedContainers, allowedTraversals,
        allowedDataLayouts, allowNewton3, selectorStrategy, tuningInterval, numSamples);
  }

  /**
   * Updates the internal container.
   * This is needed e.g. for linked-cells if particles move from one cell to another.
   * It resorts particles into appropriate cells and moves them to the halo, if necessary.
   */
  void updateContainer() { _autoTuner->getContainer()->updateContainer(); }

  /**
   * Returns a pointer to the actual container.
   * @todo do we need the whole container functionality available to the outside
   * @return container
   */
  // @todo: remove this once we are convinced all necessary container functions are wrapped
  autopas::ParticleContainer<Particle, ParticleCell> *getContainer() const { return _autoTuner->getContainer().get(); }

  /**
   * Adds a particle to the container.
   * @param p Reference to the particle to be added
   */
  void addParticle(Particle &p) { _autoTuner->getContainer()->addParticle(p); }

  /**
   * adds a particle to the container that lies in the halo region of the
   * container
   * @param haloParticle particle to be added
   */
  void addHaloParticle(Particle &haloParticle) { _autoTuner->getContainer()->addHaloParticle(haloParticle); }

  /**
   * deletes all halo particles
   */
  void deleteHaloParticles() { _autoTuner->getContainer()->deleteHaloParticles(); }

  /**
   * deletes all particles
   */
  void deleteAllParticles() { _autoTuner->getContainer()->deleteAllParticles(); }

  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @param f Functor that describes the pair-potential
   */
  template <class Functor>
  void iteratePairwise(Functor *f) {
    _autoTuner->iteratePairwise(f);
  }

  /**
   * Iterate over all particles by using
   * for(auto iter = autoPas.begin(); iter.isValid(); ++iter)
   * @param behavior the behavior of the iterator. You can specify whether to iterate over owned particles, halo
   * particles, or both.
   * @return iterator to the first particle
   */
  autopas::ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) {
    return _autoTuner->getContainer()->begin(behavior);
  }

  /**
   * iterate over all particles in a specified region
   * for(auto iter = container.getRegionIterator(lowCorner,
   * highCorner);iter.isValid();++iter)
   * @param lowerCorner lower corner of the region
   * @param higherCorner higher corner of the region
   * @param behavior the behavior of the iterator. You can specify whether to iterate over owned particles, halo
   * particles, or both.
   * @return iterator to iterate over all particles in a specific region
   */
  autopas::ParticleIteratorWrapper<Particle> getRegionIterator(
      std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) {
    return _autoTuner->getContainer()->getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * Returns the number of particles in this container.
   * @return the number of particles in this container.
   */
  unsigned long getNumberOfParticles() { return _autoTuner->getContainer()->getNumParticles(); }

  /**
   * Get the lower corner of the container.
   * @return lower corner of the container.
   */
  std::array<double, 3> getBoxMin() { return _autoTuner->getContainer()->getBoxMin(); }

  /**
   * Get the upper corner of the container.
   * @return upper corner of the container.
   */
  std::array<double, 3> getBoxMax() { return _autoTuner->getContainer()->getBoxMax(); }

  /**
   * Checks if the container needs to be updated.
   * Will return false if no lists are used.
   * This function can indicate whether you should send only halo particles or whether you should send leaving particles
   * as well.
   * @return True if the lists are valid, false if a rebuild is needed.
   */
  bool needsContainerUpdate() {
    if (_autoTuner->willRebuild()) {
      return true;
    }
    if (auto container = dynamic_cast<VerletLists<Particle> *>(_autoTuner->getContainer().get())) {
      return container->needsRebuild();
    } else {
      return true;
    }
  }

  /**
   * Lower corner of the container.
   */
  std::array<double, 3> boxMin;
  /**
   * Upper corner of the container.
   */
  std::array<double, 3> boxMax;
  /**
   * Cutoff radius to be used in this container.
   */
  double cutoff;
  /**
   * Length added to the cutoff for the verlet lists' skin.
   */
  double verletSkin;
  /**
   * Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   */
  unsigned int verletRebuildFrequency;
  /**
   * Number of timesteps after which the auto-tuner shall reevaluate all selections.
   */
  unsigned int tuningInterval;
  /**
   * Number of samples the tuner should collect for each combination.
   */
  unsigned int numSamples;
  /**
   * Strategy for the container selector.
   * For possible container choices see AutoPas::SelectorStrategy.
   */
  SelectorStrategy selectorStrategy;
  /**
   * List of container types AutoPas can choose from.
   * For possible container choices see AutoPas::ContainerOption.
   */
  std::vector<ContainerOption> allowedContainers;
  /**
   * List of traversals AutoPas can choose from.
   * For possible container choices see AutoPas::TraversalOption.
   */
  std::vector<TraversalOption> allowedTraversals;
  /**
   * List of data layouts AutoPas can choose from.
   * For possible container choices see AutoPas::DataLayoutOption.
   */
  std::vector<DataLayoutOption> allowedDataLayouts;
  /**
   * Whether AutoPas is allowed to exploit Newton's third law of motion.
   */
  bool allowNewton3;

 private:
  std::unique_ptr<autopas::AutoTuner<Particle, ParticleCell>> _autoTuner;
};  // namespace autopas
}  // namespace autopas