/**
 * @file AutoPas.h
 * Main include file for the AutoPas library.
 *
 */

#pragma once

#include <iostream>
#include <memory>
#include <set>
#include <type_traits>
#include "autopas/LogicHandler.h"
#include "autopas/autopasIncludes.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/selectors/AutoTuner.h"
#include "autopas/selectors/tuningStrategy/BayesianSearch.h"
#include "autopas/selectors/tuningStrategy/FullSearch.h"
#include "autopas/selectors/tuningStrategy/RandomSearch.h"
#include "autopas/selectors/tuningStrategy/ActiveHarmony.h"
#include "autopas/utils/NumberSet.h"

namespace autopas {

/**
 * Instance counter to help track the number of autopas instances. Needed for correct management of the logger.
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
   * Particle type to be accessible after initialization.
   */
  using Particle_t = Particle;

  /**
   * Particle Cell type to be accessible after initialization.
   */
  using ParticleCell_t = ParticleCell;

  /**
   * Define the iterator_t for simple use, also from the outside.
   * Helps to, e.g., wrap the AutoPas iterators
   */
  using iterator_t = typename autopas::IteratorTraits<Particle>::iterator_t;

  /**
   * Define the const_iterator_t for simple use, also from the outside.
   * Helps to, e.g., wrap the AutoPas iterators
   */
  using const_iterator_t = typename autopas::IteratorTraits<Particle>::const_iterator_t;

  /**
   * Constructor for the autopas class.
   * @param logOutputStream Stream where log output should go to. Default is std::out.
   */
  explicit AutoPas(std::ostream &logOutputStream = std::cout)
      : _boxMin{0, 0, 0},
        _boxMax{0, 0, 0},
        _cutoff(1.),
        _verletSkin(0.2),
        _verletRebuildFrequency(20),
        _verletClusterSize(64),
        _tuningInterval(5000),
        _numSamples(3),
        _maxEvidence(10),
        _acquisitionFunctionOption(AcquisitionFunctionOption::lowerConfidenceBound),
        _tuningStrategyOption(TuningStrategyOption::fullSearch),
        _selectorStrategy(SelectorStrategyOption::fastestAbs),
        _allowedContainers(ContainerOption::getAllOptions()),
        _allowedTraversals(TraversalOption::getAllOptions()),
        _allowedDataLayouts(DataLayoutOption::getAllOptions()),
        _allowedNewton3Options(Newton3Option::getAllOptions()),
        _allowedCellSizeFactors(std::make_unique<NumberSetFinite<double>>(std::set<double>({1.}))) {
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
    _logicHandler = std::move(other._logicHandler);
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
        _boxMin, _boxMax, _cutoff, _verletSkin, _verletClusterSize, std::move(generateTuningStrategy()),
        _selectorStrategy, _tuningInterval, _numSamples);
    _logicHandler =
        std::make_unique<autopas::LogicHandler<Particle, ParticleCell>>(*(_autoTuner.get()), _verletRebuildFrequency);
  }

  /**
   * Potentially updates the internal container.
   * On an update, halo particles are deleted, the particles are resorted into appropriate cells and particles that do
   * no longer belong into the container will be returned, the lists will be invalidated. If the internal container is
   * still valid and a rebuild of the container is not forced, this will return an empty list of particles and nothing
   * else will happen.
   * @return A pair of a vector of invalid particles that do no belong in the current container and a bool that
   * specifies whether the container was updated. If the bool is false, the vector will be an empty vector. If the
   * returned bool evaluates to true, the vector can both be empty or non-empty, depending on whether particles have
   * left the container or not.
   */
  AUTOPAS_WARN_UNUSED_RESULT
  std::pair<std::vector<Particle>, bool> updateContainer() { return _logicHandler->updateContainer(false); }

  /**
   * Forces a container update.
   * On an update, the particles are resorted into appropriate cells and will return particles that do no longer belong
   * into the container.
   * @return A vector of invalid particles that do no belong in the current container.
   */
  AUTOPAS_WARN_UNUSED_RESULT
  std::vector<Particle> updateContainerForced() { return std::get<0>(_logicHandler->updateContainer(true)); }

  /**
   * Adds a particle to the container.
   * This is only allowed if the neighbor lists are not valid.
   * @param p Reference to the particle to be added
   */
  void addParticle(Particle &p) { _logicHandler->addParticle(p); }

  /**
   * Adds or updates a particle to/in the container that lies in the halo region of the container.
   * If the neighbor lists inside of AutoPas are NOT valid, the halo particle will be added.
   * If the neighbor lists of AutoPas are valid the particle will be used to update an already existing halo particle.
   * In this case if there is no matching halo particle, the given haloParticle will be ignored.
   * @note Exceptions are thrown in the following cases:
   * 1. If the halo particle is added and it is inside of the owned domain (defined by boxmin and boxmax)of the
   * container.
   * 2. If the halo particle should be updated and the given haloParticle is too far inside of the domain (by more than
   * skin/2)
   * 3. If the halo particle should be updated, but no matching particle is found, even though the given haloParticle is
   * close enough to the domain (at most cutoff + skin/2)
   *
   * @param haloParticle particle to be added or updated
   */
  void addOrUpdateHaloParticle(Particle &haloParticle) { _logicHandler->addOrUpdateHaloParticle(haloParticle); }

  /**
   * Deletes all particles.
   * @note This invalidates the container, a rebuild is forced on the next iteratePairwise() call.
   */
  void deleteAllParticles() { _logicHandler->deleteAllParticles(); }

  /**
   * Deletes the particle behind the current iterator position.
   * @param iter Needs to be a modify-able iterator.
   */
  void deleteParticle(ParticleIteratorWrapper<Particle, true> &iter) { _logicHandler->deleteParticle(iter); }

  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @param f Functor that describes the pair-potential.
   * @return true if this was a tuning iteration.
   */
  template <class Functor>
  bool iteratePairwise(Functor *f) {
    static_assert(not std::is_same<Functor, autopas::Functor<Particle, ParticleCell>>::value,
                  "The static type of Functor in iteratePairwise is not allowed to be autopas::Functor. Please use the "
                  "derived type instead, e.g. by using a dynamic_cast.");
    if (f->getCutoff() > this->getCutoff()) {
      utils::ExceptionHandler::exception("Functor cutoff ({}) must not be larger than container cutoff ({})",
                                         f->getCutoff(), this->getCutoff());
    }
    return _logicHandler->iteratePairwise(f);
  }

  /**
   * Iterate over all particles by using
   * for(auto iter = autoPas.begin(); iter.isValid(); ++iter)
   * @param behavior the behavior of the iterator. You can specify whether to iterate over owned particles, halo
   * particles, or both.
   * @return iterator to the first particle.
   */
  iterator_t begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) {
    return _logicHandler->begin(behavior);
  }

  /**
   * @copydoc begin()
   * @note const version
   */
  const_iterator_t begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const {
    return std::as_const(*_logicHandler).begin(behavior);
  }

  /**
   * @copydoc begin()
   * @note cbegin will guarantee to return a const_iterator.
   */
  const_iterator_t cbegin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const { return begin(behavior); }

  /**
   * End of the iterator.
   * This returns a bool, which is false to allow range-based for loops.
   * @return false
   */
  [[nodiscard]] constexpr bool end() const { return false; }

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
  iterator_t getRegionIterator(std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
                               IteratorBehavior behavior = IteratorBehavior::haloAndOwned) {
    return _logicHandler->getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * @copydoc getRegionIterator()
   * @note const version
   */
  const_iterator_t getRegionIterator(std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
                                     IteratorBehavior behavior = IteratorBehavior::haloAndOwned) const {
    return _logicHandler->getRegionIterator(lowerCorner, higherCorner, behavior);
  }

  /**
   * Returns the number of particles in this container.
   * @param behavior Tells this function to report the number of halo, owned or all particles.
   * @return the number of particles in this container.
   */
  unsigned long getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::ownedOnly) const {
    switch (behavior) {
      case IteratorBehavior::ownedOnly: {
        return _logicHandler->getNumParticlesOwned();
      }
      case IteratorBehavior::haloOnly: {
        return _logicHandler->getNumParticlesHalo();
      }
      case IteratorBehavior::haloAndOwned: {
        return _logicHandler->getNumParticlesOwned() + _logicHandler->getNumParticlesHalo();
      }
    }
    return 0;
  }

  /**
   * Returns the type of the currently used container.
   * @return The type of the used container is returned.
   */
  unsigned long getContainerType() const { return _autoTuner->getContainer()->getContainerType(); }

  /**
   * Get the lower corner of the container.
   * @return lower corner of the container.
   */
  std::array<double, 3> getBoxMin() const { return _autoTuner->getContainer()->getBoxMin(); }

  /**
   * Get the upper corner of the container.
   * @return upper corner of the container.
   */
  std::array<double, 3> getBoxMax() const { return _autoTuner->getContainer()->getBoxMax(); }

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
   * Get allowed cell size factors (only relevant for LinkedCells, VerletLists and VerletListsCells).
   * @return
   */
  const NumberSet<double> &getAllowedCellSizeFactors() const { return *_allowedCellSizeFactors; }

  /**
   * Set allowed cell size factors (only relevant for LinkedCells, VerletLists and VerletListsCells).
   * @param allowedCellSizeFactors
   */
  void setAllowedCellSizeFactors(const NumberSet<double> &allowedCellSizeFactors) {
    if (allowedCellSizeFactors.getMin() <= 0.0) {
      AutoPasLog(error, "cell size <= 0.0");
      utils::ExceptionHandler::exception("Error: cell size <= 0.0!");
    }
    AutoPas::_allowedCellSizeFactors = std::move(allowedCellSizeFactors.clone());
  }

  /**
   * Set allowed cell size factors to one element (only relevant for LinkedCells, VerletLists and VerletListsCells).
   * @param cellSizeFactor
   */
  void setCellSizeFactor(double cellSizeFactor) {
    if (cellSizeFactor <= 0.0) {
      AutoPasLog(error, "cell size <= 0.0: {}", cellSizeFactor);
      utils::ExceptionHandler::exception("Error: cell size <= 0.0!");
    }
    AutoPas::_allowedCellSizeFactors = std::make_unique<NumberSetFinite<double>>(std::set<double>{cellSizeFactor});
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
   * Get Verlet cluster size.
   * @return
   */
  unsigned int getVerletClusterSize() const { return _verletClusterSize; }

  /**
   * Set Verlet cluster size.
   * @param verletClusterSize
   */
  void setVerletClusterSize(unsigned int verletClusterSize) { AutoPas::_verletClusterSize = verletClusterSize; }

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
   * Get maximum number of evidence for tuning
   * @return
   */
  unsigned int getMaxEvidence() const { return _maxEvidence; }

  /**
   * Set maximum number of evidence for tuning
   * @param maxEvidence
   */
  void setMaxEvidence(unsigned int maxEvidence) { AutoPas::_maxEvidence = maxEvidence; }

  /**
   * Get acquisition function used for tuning
   * @return
   */
  AcquisitionFunctionOption getAcquisitionFunction() const { return _acquisitionFunctionOption; }

  /**
   * Set acquisition function for tuning
   * @param acqFun acquisition function
   */
  void setAcquisitionFunction(AcquisitionFunctionOption acqFun) { AutoPas::_acquisitionFunctionOption = acqFun; }

  /**
   * Get the selector configuration strategy.
   * @return
   */
  SelectorStrategyOption getSelectorStrategy() const { return _selectorStrategy; }

  /**
   * Set the selector configuration strategy.
   * For possible selector strategy choices see AutoPas::SelectorStrategy.
   * @param selectorStrategy
   */
  void setSelectorStrategy(SelectorStrategyOption selectorStrategy) { AutoPas::_selectorStrategy = selectorStrategy; }

  /**
   * Get the list of allowed containers.
   * @return
   */
  const std::set<ContainerOption> &getAllowedContainers() const { return _allowedContainers; }

  /**
   * Set the list of allowed containers.
   * For possible container choices see AutoPas::ContainerOption.
   * @param allowedContainers
   */
  void setAllowedContainers(const std::set<ContainerOption> &allowedContainers) {
    AutoPas::_allowedContainers = allowedContainers;
  }

  /**
   * Get the list of allowed traversals.
   * @return
   */
  const std::set<TraversalOption> &getAllowedTraversals() const { return _allowedTraversals; }

  /**
   * Set the list of allowed traversals.
   * For possible traversals choices see AutoPas::TraversalOption.
   * @param allowedTraversals
   */
  void setAllowedTraversals(const std::set<TraversalOption> &allowedTraversals) {
    AutoPas::_allowedTraversals = allowedTraversals;
  }

  /**
   * Get the list of allowed data layouts.
   * @return
   */
  const std::set<DataLayoutOption> &getAllowedDataLayouts() const { return _allowedDataLayouts; }

  /**
   * Set the list of allowed data layouts.
   * For possible data layout choices see AutoPas::DataLayoutOption.
   * @param allowedDataLayouts
   */
  void setAllowedDataLayouts(const std::set<DataLayoutOption> &allowedDataLayouts) {
    AutoPas::_allowedDataLayouts = allowedDataLayouts;
  }

  /**
   * Get the list of allowed newton 3 options.
   * @return
   */
  const std::set<Newton3Option> &getAllowedNewton3Options() const { return _allowedNewton3Options; }

  /**
   * Set the list of allowed newton 3 options.
   * For possible newton 3 choices see AutoPas::Newton3Option
   * @param allowedNewton3Options
   */
  void setAllowedNewton3Options(const std::set<Newton3Option> &allowedNewton3Options) {
    AutoPas::_allowedNewton3Options = allowedNewton3Options;
  }

  /**
   * Getter for the currently selected configuration.
   * @return Configuration object currently used.
   */
  const Configuration getCurrentConfig() const { return _autoTuner->getCurrentConfig(); }

  /**
   * Getter for the tuning strategy option.
   * @return
   */
  TuningStrategyOption getTuningStrategyOption() const { return _tuningStrategyOption; }

  /**
   * Setter for the tuning strategy option.
   * @param tuningStrategyOption
   */
  void setTuningStrategyOption(TuningStrategyOption tuningStrategyOption) {
    _tuningStrategyOption = tuningStrategyOption;
  }

 private:
  /**
   * Generates a new Tuning Strategy object from the member variables of this autopas object.
   * @return Pointer to the tuning strategy object or the nullpointer if an exception was suppressed.
   */
  std::unique_ptr<TuningStrategyInterface> generateTuningStrategy() {
    // clang compiler bug requires static cast
    switch (static_cast<TuningStrategyOption>(_tuningStrategyOption)) {
      case TuningStrategyOption::randomSearch: {
        return std::make_unique<RandomSearch>(_allowedContainers, *_allowedCellSizeFactors, _allowedTraversals,
                                              _allowedDataLayouts, _allowedNewton3Options, _maxEvidence);
      }
      case TuningStrategyOption::fullSearch: {
        if (not _allowedCellSizeFactors->isFinite()) {
          autopas::utils::ExceptionHandler::exception(
              "AutoPas::generateTuningStrategy: fullSearch can not handle infinite cellSizeFactors!");
          return nullptr;
        }

        return std::make_unique<FullSearch>(_allowedContainers, _allowedCellSizeFactors->getAll(), _allowedTraversals,
                                            _allowedDataLayouts, _allowedNewton3Options);
      }

      case TuningStrategyOption::bayesianSearch: {
        return std::make_unique<BayesianSearch>(_allowedContainers, *_allowedCellSizeFactors, _allowedTraversals,
                                                _allowedDataLayouts, _allowedNewton3Options, _maxEvidence,
                                                _acquisitionFunctionOption);
      }

      case TuningStrategyOption::activeHarmony: {
        return std::make_unique<ActiveHarmony>(_allowedContainers, NumberInterval(_allowedCellSizeFactors->getMin(), _allowedCellSizeFactors->getMax()), _allowedTraversals, _allowedDataLayouts, _allowedNewton3Options);
      }
    }

    autopas::utils::ExceptionHandler::exception("AutoPas::generateTuningStrategy: Unknown tuning strategy {}!",
                                                _tuningStrategyOption);
    return nullptr;
  }

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
   * Length added to the cutoff for the verlet lists' skin.
   */
  double _verletSkin;
  /**
   * Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   */
  unsigned int _verletRebuildFrequency;
  /**
   * Specifies the size of clusters for verlet lists.
   */
  unsigned int _verletClusterSize;
  /**
   * Number of timesteps after which the auto-tuner shall reevaluate all selections.
   */
  unsigned int _tuningInterval;
  /**
   * Number of samples the tuner should collect for each combination.
   */
  unsigned int _numSamples;
  /**
   * Tuning Strategies which work on a fixed number of evidence should use this value.
   */
  unsigned int _maxEvidence;
  /**
   * Acquisition function used for tuning.
   * For possible acquisition function choices see AutoPas::AcquisitionFunction.
   */
  AcquisitionFunctionOption _acquisitionFunctionOption;

  /**
   * Strategy option for the auto tuner.
   * For possible tuning strategy choices see AutoPas::TuningStrategy.
   */
  TuningStrategyOption _tuningStrategyOption;

  /**
   * Strategy for the configuration selector.
   * For possible container choices see AutoPas::SelectorStrategy.
   */
  SelectorStrategyOption _selectorStrategy;
  /**
   * List of container types AutoPas can choose from.
   * For possible container choices see AutoPas::ContainerOption.
   */
  std::set<ContainerOption> _allowedContainers;
  /**
   * List of traversals AutoPas can choose from.
   * For possible container choices see AutoPas::TraversalOption.
   */
  std::set<TraversalOption> _allowedTraversals;
  /**
   * List of data layouts AutoPas can choose from.
   * For possible container choices see AutoPas::DataLayoutOption.
   */
  std::set<DataLayoutOption> _allowedDataLayouts;
  /**
   * Whether AutoPas is allowed to exploit Newton's third law of motion.
   */
  std::set<Newton3Option> _allowedNewton3Options;
  /**
   * Cell size factor to be used in this container (only relevant for LinkedCells, VerletLists and VerletListsCells).
   */
  std::unique_ptr<NumberSet<double>> _allowedCellSizeFactors;

  /**
   * LogicHandler of autopas.
   */
  std::unique_ptr<autopas::LogicHandler<Particle, ParticleCell>> _logicHandler;

  /**
   * This is the AutoTuner that owns the container, ...
   */
  std::unique_ptr<autopas::AutoTuner<Particle, ParticleCell>> _autoTuner;
};  // class AutoPas
}  // namespace autopas
