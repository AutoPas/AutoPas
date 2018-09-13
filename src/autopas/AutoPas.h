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
  AutoPas() {
    // initialize the Logger
    autopas::Logger::create();
    // The logger is normally only flushed on successful program termination.
    // This line ensures flushing when log messages of level warning or more severe are created.
    AutoPasLogger->flush_on(spdlog::level::warn);
  }

  ~AutoPas() {
    // remove the Logger from the registry
    autopas::Logger::unregister();
  }

  /**
   * Initialize the particle container.
   *
   * For possible container choices see AutoPas::ContainerOption.
   *
   * @param boxMin Lower corner of the container.
   * @param boxMax Upper corner of the container.
   * @param cutoff  Cutoff radius to be used in this container.
   * @param verletSkin Length added to the cutoff for the verlet lists' skin.
   * @param verletRebuildFrequency Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   * @param allowedContainers List of container types AutoPas can choose from.
   * @param allowedTraversals List of traversals AutoPas can choose from.
   * @param tuningInterval Number of timesteps after which the auto-tuner shall reevaluate all selections.
   */
  void init(std::array<double, 3> boxMin, std::array<double, 3> boxMax, double cutoff, double verletSkin,
            unsigned int verletRebuildFrequency,
            const std::vector<autopas::ContainerOptions> &allowedContainers = autopas::allContainerOptions,
            const std::vector<autopas::TraversalOptions> &allowedTraversals = autopas::allTraversalOptions,
            unsigned int tuningInterval = 100) {
    _autoTuner = std::make_unique<autopas::AutoTuner<Particle, ParticleCell>>(boxMin, boxMax, cutoff, verletSkin,
                                                                              verletRebuildFrequency, allowedContainers,
                                                                              allowedTraversals, tuningInterval);
  }

  /**
   * @overload
   *
   */
  void init(std::array<double, 3> boxSize, double cutoff, double verletSkin = 0,
            unsigned int verletRebuildFrequency = 1,
            const std::vector<autopas::ContainerOptions> &allowedContainers = autopas::allContainerOptions,
            const std::vector<autopas::TraversalOptions> &allowedTraversals = autopas::allTraversalOptions,
            unsigned int tuningInterval = 100) {
    init({0, 0, 0}, boxSize, cutoff, verletSkin, verletRebuildFrequency, allowedContainers, allowedTraversals,
         tuningInterval);
  }

  /**
   * Returns a pointer to the actual container.
   * @todo do we need the whole container functionality available to the outside
   * @return container
   */
  // TODO: remove this once we are convinced all necessary container functions are wrapped
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
  void addHaloParticle(Particle &haloParticle) { _autoTuner->getContainer()->addHaloParticle(haloParticle); };

  /**
   * deletes all halo particles
   */
  void deleteHaloParticles() { _autoTuner->getContainer()->deleteHaloParticles(); };

  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @param f Functor that describes the pair-potential
   * @param dataLayoutOption useSoA Bool to decide if SoA or AoS should be used.
   */
  template <class Functor>
  void iteratePairwise(Functor *f, autopas::DataLayoutOption dataLayoutOption) {
    // @todo remove this and let is be handled via a selector
    _autoTuner->iteratePairwise(f, dataLayoutOption);
  }

  /**
   * iterate over all particles by using
   * for(auto iter = autoPas.begin(); iter.isValid(); ++iter)
   * @return iterator to the first particle
   */
  autopas::ParticleIteratorWrapper<Particle> begin() { return _autoTuner->getContainer()->begin(); }

  /**
   * iterate over all particles in a specified region
   * for(auto iter = container.getRegionIterator(lowCorner,
   * highCorner);iter.isValid();++iter)
   * @param lowerCorner lower corner of the region
   * @param higherCorner higher corner of the region
   * @return iterator to iterate over all particles in a specific region
   */
  autopas::ParticleIteratorWrapper<Particle> getRegionIterator(std::array<double, 3> lowerCorner,
                                                               std::array<double, 3> higherCorner) {
    return _autoTuner->getContainer()->getRegionIterator(lowerCorner, higherCorner);
  }
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

 private:
  std::unique_ptr<autopas::AutoTuner<Particle, ParticleCell>> _autoTuner;
};

}  // namespace autopas