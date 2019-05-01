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
  /**
   * Constructor for the autopas class.
   * @param logOutputStream Stream where log output should go to. Default is std::out.
   */
  AutoPas(std::ostream &logOutputStream = std::cout)
      : _boxMin{0, 0, 0},
        _boxMax{0, 0, 0},
        _cutoff(1),
        _cellSizeFactor(1),
        _verletSkin(0.2),
        _verletRebuildFrequency(20),
        _tuningInterval(5000),
        _numSamples(3),
        _selectorStrategy(SelectorStrategy::fastestAbs),
        _allowedContainers(allContainerOptions),
        _allowedTraversals(allTraversalOptions),
        _allowedDataLayouts(allDataLayoutOptions),
        _allowedNewton3Options(allNewton3Options) {
    // count the number of autopas instances. This is needed to ensure that the autopas
    // logger is not unregistered while other instances are still using it.
    _instanceCounter++;
    // remove potentially existing logger
    autopas::Logger::unregister();
    // initialize the Logger
    autopas::Logger::create(logOutputStream);
    // The logger is normally only flushed on successful program termination.
    // This line ensures flushing when log messages of level warning or more severe are created.
    autopas::Logger::get()->flush_on(spdlog::level::warn);
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
   * Initialize AutoPas. This will completely reset the container and remove all containing particles!
   *
   * This function needs to be called before any other function (except setters) on the AutoPas object.
   *
   * Changing any of the member options only takes effect when init is called.
   *
   */
  void init() {
    _autoTuner = std::make_unique<autopas::AutoTuner<Particle, ParticleCell>>(
        _boxMin, _boxMax, _cutoff, _cellSizeFactor, _verletSkin, _verletRebuildFrequency, _allowedContainers,
        _allowedTraversals, _allowedDataLayouts, _allowedNewton3Options, _selectorStrategy, _tuningInterval,
        _numSamples);
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
   * @param f Functor that describes the pair-potential.
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
   * Set coordinates of the lower corner of the domain.
   * @param boxMin
   */
  void setBoxMin(const std::array<double, 3> &boxMin) { AutoPas::_boxMin = boxMin; }

  /**
   * Set coordinates of the upper corner of the domain.
   * @param boxMax
   */
  void setBoxMax(const std::array<double, 3> &boxMax) { AutoPas::_boxMax = boxMax; }

  /**
   * Get cutoff radius.
   * @return
   */
  double getCutoff() const { return _cutoff; }

  /**
   * Set cutoff radius.
   * @param cutoff
   */
  void setCutoff(double cutoff) {
    if (cutoff <= 0.0) {
      AutoPasLog(error, "Cutoff <= 0.0: {}", cutoff);
      utils::ExceptionHandler::exception("Error: Cutoff <= 0.0!");
    }
    AutoPas::_cutoff = cutoff;
  }

  /**
   * Get cell size factor.
   * @return
   */
  double getCellSizeFactor() const { return _cellSizeFactor; }

  /**
   * Set cell size factor.
   * @param cellSizeFactor
   */
  void setCellSizeFactor(double cellSizeFactor) {
    if (cellSizeFactor <= 0.0) {
      AutoPasLog(error, "cell size <= 0.0: {}", cellSizeFactor);
      utils::ExceptionHandler::exception("Error: cell size <= 0.0!");
    }
    AutoPas::_cellSizeFactor = cellSizeFactor;
  }

  /**
   * Get length added to the cutoff for the Verlet lists' skin.
   * @return
   */
  double getVerletSkin() const { return _verletSkin; }

  /**
   * Set length added to the cutoff for the Verlet lists' skin.
   * @param verletSkin
   */
  void setVerletSkin(double verletSkin) { AutoPas::_verletSkin = verletSkin; }

  /**
   * Get Verlet rebuild frequency.
   * @return
   */
  unsigned int getVerletRebuildFrequency() const { return _verletRebuildFrequency; }

  /**
   * Set Verlet rebuild frequency.
   * @param verletRebuildFrequency
   */
  void setVerletRebuildFrequency(unsigned int verletRebuildFrequency) {
    AutoPas::_verletRebuildFrequency = verletRebuildFrequency;
  }

  /**
   * Get tuning interval.
   * @return
   */
  unsigned int getTuningInterval() const { return _tuningInterval; }

  /**
   * Set tuning interval.
   * @param tuningInterval
   */
  void setTuningInterval(unsigned int tuningInterval) { AutoPas::_tuningInterval = tuningInterval; }

  /**
   * Get number of samples taken per configuration during the tuning.
   * @return
   */
  unsigned int getNumSamples() const { return _numSamples; }

  /**
   * Set number of samples taken per configuration during the tuning.
   * @param numSamples
   */
  void setNumSamples(unsigned int numSamples) { AutoPas::_numSamples = numSamples; }

  /**
   * Get the selector configuration strategy.
   * @return
   */
  SelectorStrategy getSelectorStrategy() const { return _selectorStrategy; }

  /**
   * Set the selector configuration strategy.
   * For possible selector strategy choices see AutoPas::SelectorStrategy.
   * @param selectorStrategy
   */
  void setSelectorStrategy(SelectorStrategy selectorStrategy) { AutoPas::_selectorStrategy = selectorStrategy; }

  /**
   * Get the list of allowed containers.
   * @return
   */
  const std::vector<ContainerOption> &getAllowedContainers() const { return _allowedContainers; }

  /**
   * Set the list of allowed containers.
   * For possible container choices see AutoPas::ContainerOption.
   * @param allowedContainers
   */
  void setAllowedContainers(const std::vector<ContainerOption> &allowedContainers) {
    AutoPas::_allowedContainers = allowedContainers;
  }

  /**
   * Get the list of allowed traversals.
   * @return
   */
  const std::vector<TraversalOption> &getAllowedTraversals() const { return _allowedTraversals; }

  /**
   * Set the list of allowed traversals.
   * For possible traversals choices see AutoPas::TraversalOption.
   * @param allowedTraversals
   */
  void setAllowedTraversals(const std::vector<TraversalOption> &allowedTraversals) {
    AutoPas::_allowedTraversals = allowedTraversals;
  }

  /**
   * Get the list of allowed data layouts.
   * @return
   */
  const std::vector<DataLayoutOption> &getAllowedDataLayouts() const { return _allowedDataLayouts; }

  /**
   * Set the list of allowed data layouts.
   * For possible data layout choices see AutoPas::DataLayoutOption.
   * @param allowedDataLayouts
   */
  void setAllowedDataLayouts(const std::vector<DataLayoutOption> &allowedDataLayouts) {
    AutoPas::_allowedDataLayouts = allowedDataLayouts;
  }

  /**
   * Get the list of allowed newton 3 options.
   * @return
   */
  const std::vector<Newton3Option> &getAllowedNewton3Options() const { return _allowedNewton3Options; }

  /**
   * Set the list of allowed newton 3 options.
   * For possible newton 3 choices see AutoPas::Newton3Option
   * @param allowedNewton3Options
   */
  void setAllowedNewton3Options(const std::vector<Newton3Option> &allowedNewton3Options) {
    AutoPas::_allowedNewton3Options = allowedNewton3Options;
  }

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
   * Getter for the currently selected configuration.
   * @return Configuration object currently used.
   */
  const Configuration getCurrentConfig() const { return _autoTuner->getCurrentConfig(); }

 private:
  /**
   * Lower corner of the container.
   */
  std::array<double, 3> _boxMin;
  /**
   * Upper corner of the container.
   */
  std::array<double, 3> _boxMax;
  /**
   * Cutoff radius to be used in this container.
   */
  double _cutoff;
  /**
   * Cell size factor to be used in this container.
   */
  double _cellSizeFactor;
  /**
   * Length added to the cutoff for the verlet lists' skin.
   */
  double _verletSkin;
  /**
   * Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   */
  unsigned int _verletRebuildFrequency;
  /**
   * Number of timesteps after which the auto-tuner shall reevaluate all selections.
   */
  unsigned int _tuningInterval;
  /**
   * Number of samples the tuner should collect for each combination.
   */
  unsigned int _numSamples;
  /**
   * Strategy for the configuration selector.
   * For possible container choices see AutoPas::SelectorStrategy.
   */
  SelectorStrategy _selectorStrategy;
  /**
   * List of container types AutoPas can choose from.
   * For possible container choices see AutoPas::ContainerOption.
   */
  std::vector<ContainerOption> _allowedContainers;
  /**
   * List of traversals AutoPas can choose from.
   * For possible container choices see AutoPas::TraversalOption.
   */
  std::vector<TraversalOption> _allowedTraversals;
  /**
   * List of data layouts AutoPas can choose from.
   * For possible container choices see AutoPas::DataLayoutOption.
   */
  std::vector<DataLayoutOption> _allowedDataLayouts;
  /**
   * Whether AutoPas is allowed to exploit Newton's third law of motion.
   */
  std::vector<Newton3Option> _allowedNewton3Options;

  std::unique_ptr<autopas::AutoTuner<Particle, ParticleCell>> _autoTuner;
};  // namespace autopas
}  // namespace autopas
