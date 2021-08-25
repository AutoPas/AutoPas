/**
 * @file AutoPasDecl.h
 *
 */

#pragma once

#include <memory>
#include <set>

#include "autopas/iterators/ParticleIteratorWrapper.h"
#include "autopas/options//ExtrapolationMethodOption.h"
#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/MPIStrategyOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/selectors/Configuration.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/WrapMPI.h"

// These next three includes have dependencies to all of AutoPas and thus are moved here from AutoPasDecl.h.
#include "autopas/LogicHandler.h"
#include "autopas/selectors/AutoTuner.h"
#include "autopas/selectors/tuningStrategy/TuningStrategyFactory.h"

namespace autopas {

// Forward declare AutoTuner and LogicHandler so that including this header does not include the whole library with all
// containers and traversals.
template <class Particle>
class AutoTuner;

template <class Particle>
class LogicHandler;

/**
 * The AutoPas class is intended to be the main point of Interaction for the user.
 * It acts as an interface from where all features of the library can be triggered and configured.
 * @tparam Particle Class for particles
 * @tparam ParticleCell Class for the particle cells
 */
template <class Particle>
class AutoPas {
 public:
  /**
   * Particle type to be accessible after initialization.
   */
  using Particle_t = Particle;

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
  explicit AutoPas(std::ostream &logOutputStream = std::cout);

  ~AutoPas();

  /**
   * Move assignment operator
   * @param other
   * @return
   */
  AutoPas &operator=(AutoPas &&other) noexcept;

  /**
   * Initialize AutoPas. This will completely reset the container and remove all containing particles!
   *
   * This function needs to be called before any other function (except setters) on the AutoPas object.
   *
   * Changing any of the member options only takes effect when init is called.
   *
   */
  void init();

  /**
   * Resizes the bounding box of the AutoPas object.
   * @param boxMin
   * @param boxMax
   * @return Vector of particles that are outside the box after the resize.
   */
  std::vector<Particle> resizeBox(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax);

  /**
   * Force the internal tuner to enter a new tuning phase upon the next call to iteratePairwise().
   */
  void forceRetune();

  /**
   * Free the AutoPas MPI communicator.
   * To be called before MPI_Finalize.
   * If no MPI is used just call this at the end of the program.
   */
  void finalize() {
    if (not _externalMPICommunicator) {
      AutoPas_MPI_Comm_free(&_autopasMPICommunicator);
    }
  }

  /**
   * Potentially updates the internal container.
   * On an update, halo particles are deleted, the particles are resorted into appropriate cells and particles that do
   * no longer belong into the container will be returned, the lists will be invalidated. If the internal container is
   * still valid and a rebuild of the container is not forced, this will return an empty list of particles and nothing
   * else will happen.
   * @param forced specifies whether an update of the container is enforced.
   * @return A pair of a vector of invalid particles that do no belong in the current container and a bool that
   * specifies whether the container was updated. If the bool is false, the vector will be an empty vector. If the
   * returned bool evaluates to true, the vector can both be empty or non-empty, depending on whether particles have
   * left the container or not.
   */
  [[nodiscard]] std::pair<std::vector<Particle>, bool> updateContainer(bool forced = false);

  /**
   * Adds a particle to the container.
   * This is only allowed if the neighbor lists are not valid.
   * @param p Reference to the particle to be added
   */
  void addParticle(const Particle &p);

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
  void addOrUpdateHaloParticle(const Particle &haloParticle);

  /**
   * Deletes all particles.
   * @note This invalidates the container, a rebuild is forced on the next iteratePairwise() call.
   */
  void deleteAllParticles();

  /**
   * Deletes the particle behind the current iterator position.
   * @param iter Needs to be a modify-able iterator.
   */
  void deleteParticle(ParticleIteratorWrapper<Particle, true> &iter);

  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @param f Functor that describes the pair-potential.
   * @return true if this was a tuning iteration.
   */
  template <class Functor>
  bool iteratePairwise(Functor *f);

  /**
   * Iterate over all particles by using
   * for(auto iter = autoPas.begin(); iter.isValid(); ++iter)
   * @param behavior The behavior of the iterator. You can specify whether to iterate over owned particles, halo
   * particles, or both.
   * @return iterator to the first particle.
   */
  iterator_t begin(IteratorBehavior behavior = IteratorBehavior::ownedOrHalo);

  /**
   * @copydoc begin()
   * @note const version
   */
  const_iterator_t begin(IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const;

  /**
   * execute code on all particles in parallel as defined by a lambda function
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   * @note not actually parallel until kokkos integration
   */
  template <typename Lambda>
  void forEachParallel(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    _logicHandler->forEachParallel(forEachLambda, behavior);
  }

  /**
   * @copydoc forEachParallel()
   * @note const version
   */
  template <typename Lambda>
  void forEachParallel(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    std::as_const(*_logicHandler).forEachParallel(forEachLambda, behavior);
  }

  /**
   * execute code on all particles as defined by a lambda function
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    _logicHandler->forEach(forEachLambda, behavior);
  }

  /**
   * @copydoc forEach()
   * @note const version
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    std::as_const(*_logicHandler).forEach(forEachLambda, behavior);
  }

  /**
   * @copydoc begin()
   * @note cbegin will guarantee to return a const_iterator.
   */
  const_iterator_t cbegin(IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const { return begin(behavior); }

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
                               IteratorBehavior behavior = IteratorBehavior::ownedOrHalo);
  /**
   * @copydoc getRegionIterator()
   * @note const version
   */
  const_iterator_t getRegionIterator(std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
                                     IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const;

  /**
   * execute code on all particles in a certain region in parallel as defined by a lambda function
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   * @note not actually parallel until kokkos integration
   */
  template <typename Lambda>
  void forEachInRegionParallel(Lambda forEachLambda, std::array<double, 3> lowerCorner,
                               std::array<double, 3> higherCorner,
                               IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    _logicHandler->forEachInRegionParallel(forEachLambda, lowerCorner, higherCorner, behavior);
  }

  /**
   * @copydoc forEachInRegionParallel()
   * @note const version
   */
  template <typename Lambda>
  void forEachInRegionParallel(Lambda forEachLambda, std::array<double, 3> lowerCorner,
                               std::array<double, 3> higherCorner,
                               IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    std::as_const(*_logicHandler).forEachInRegionParallel(forEachLambda, lowerCorner, higherCorner, behavior);
  }

  /**
   * execute code on all particles in a certain region as defined by a lambda function
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
                       IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    _logicHandler->forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
  }

  /**
   * @copydoc forEachInRegion()
   * @note const version
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
                       IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    std::as_const(*_logicHandler).forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
  }

  /**
   * Returns the number of particles in this container.
   * @param behavior Tells this function to report the number of halo, owned or all particles.
   * @return the number of particles in this container.
   */
  [[nodiscard]] size_t getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::owned) const;

  /**
   * Returns the type of the currently used container.
   * @return The type of the used container is returned.
   */
  [[nodiscard]] unsigned long getContainerType() const;

  /**
   * Get the lower corner of the container without the halo.
   * @return lower corner of the container.
   */
  [[nodiscard]] std::array<double, 3> getBoxMin() const;

  /**
   * Get the upper corner of the container without the halo.
   * @return upper corner of the container.
   */
  [[nodiscard]] std::array<double, 3> getBoxMax() const;

  /**
   * Set coordinates of the lower corner of the domain.
   * @param boxMin
   */
  void setBoxMin(const std::array<double, 3> &boxMin) { _boxMin = boxMin; }

  /**
   * Set coordinates of the upper corner of the domain.
   * @param boxMax
   */
  void setBoxMax(const std::array<double, 3> &boxMax) { _boxMax = boxMax; }

  /**
   * Get cutoff radius.
   * @return
   */
  [[nodiscard]] double getCutoff() const { return _cutoff; }

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
  [[nodiscard]] const NumberSet<double> &getAllowedCellSizeFactors() const { return *_allowedCellSizeFactors; }

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
  [[nodiscard]] double getVerletSkin() const { return _verletSkin; }

  /**
   * Set length added to the cutoff for the Verlet lists' skin.
   * @param verletSkin
   */
  void setVerletSkin(double verletSkin) { AutoPas::_verletSkin = verletSkin; }

  /**
   * Get Verlet rebuild frequency.
   * @return
   */
  [[nodiscard]] unsigned int getVerletRebuildFrequency() const { return _verletRebuildFrequency; }

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
  [[nodiscard]] unsigned int getVerletClusterSize() const { return _verletClusterSize; }

  /**
   * Set Verlet cluster size.
   * @param verletClusterSize
   */
  void setVerletClusterSize(unsigned int verletClusterSize) { AutoPas::_verletClusterSize = verletClusterSize; }

  /**
   * Get tuning interval.
   * @return
   */
  [[nodiscard]] unsigned int getTuningInterval() const { return _tuningInterval; }

  /**
   * Set tuning interval.
   * @param tuningInterval
   */
  void setTuningInterval(unsigned int tuningInterval) { AutoPas::_tuningInterval = tuningInterval; }

  /**
   * Get number of samples taken per configuration during the tuning.
   * @return
   */
  [[nodiscard]] unsigned int getNumSamples() const { return _numSamples; }

  /**
   * Set number of samples taken per configuration during the tuning.
   * @param numSamples
   */
  void setNumSamples(unsigned int numSamples) { AutoPas::_numSamples = numSamples; }

  /**
   * Get maximum number of evidence for tuning
   * @return
   */
  [[nodiscard]] unsigned int getMaxEvidence() const { return _maxEvidence; }

  /**
   * Set maximum number of evidence for tuning
   * @param maxEvidence
   */
  void setMaxEvidence(unsigned int maxEvidence) { AutoPas::_maxEvidence = maxEvidence; }

  /**
   * Get the range for the optimum in which has to be to be tested
   * @return
   */
  [[nodiscard]] double getRelativeOptimumRange() const { return _relativeOptimumRange; }

  /**
   * Set the range for the optimum in which has to be to be tested
   * @param relativeOptimumRange
   */
  void setRelativeOptimumRange(double relativeOptimumRange) { AutoPas::_relativeOptimumRange = relativeOptimumRange; }

  /**
   * Get the maximum number of tuning phases a configuration can not be tested.
   * @return
   */
  [[nodiscard]] unsigned int getMaxTuningPhasesWithoutTest() const { return _maxTuningPhasesWithoutTest; }

  /**
   * For Predictive tuning: Set the relative cutoff for configurations to be blacklisted.
   * E.g. 2.5 means all configurations that take 2.5x the time of the optimum are blacklisted.
   * @param maxTuningPhasesWithoutTest
   */
  void setMaxTuningPhasesWithoutTest(unsigned int maxTuningPhasesWithoutTest) {
    AutoPas::_maxTuningPhasesWithoutTest = maxTuningPhasesWithoutTest;
  }

  /**
   * For Predictive tuning: Get the relative cutoff for configurations to be blacklisted.
   * E.g. 2.5 means all configurations that take 2.5x the time of the optimum are blacklisted.
   * @return
   */
  [[nodiscard]] double getRelativeBlacklistRange() const { return _relativeBlacklistRange; }

  /**
   * Set the range of the configurations that are not going to be blacklisted.
   * @param relativeBlacklistRange
   */
  void setRelativeBlacklistRange(double relativeBlacklistRange) {
    AutoPas::_relativeBlacklistRange = relativeBlacklistRange;
  }

  /**
   * Get the number of tests that need to have happened for a configuration until the first predictions are going to be
   * calculated.
   * @return
   */
  [[nodiscard]] unsigned int getEvidenceFirstPrediction() const { return _evidenceFirstPrediction; }

  /**
   * Set the number of tests that need to have happened for a configuration until the first predictions are going to be
   * calculated.
   * @param evidenceFirstPrediction
   */
  void setEvidenceFirstPrediction(unsigned int evidenceFirstPrediction) {
    AutoPas::_evidenceFirstPrediction = evidenceFirstPrediction;
  }

  /**
   * Get acquisition function used for tuning
   * @return
   */
  [[nodiscard]] AcquisitionFunctionOption getAcquisitionFunction() const { return _acquisitionFunctionOption; }

  /**
   * Set acquisition function for tuning.
   * For possible acquisition function choices see options::AcquisitionFunctionOption::Value.
   * @note This function is only relevant for the bayesian based searches.
   * @param acqFun acquisition function
   */
  void setAcquisitionFunction(AcquisitionFunctionOption acqFun) { AutoPas::_acquisitionFunctionOption = acqFun; }

  /**
   * Get extrapolation method for the prediction of the configuration performance.
   * @return
   */
  ExtrapolationMethodOption getExtrapolationMethodOption() const { return _extrapolationMethodOption; }

  /**
   * Set extrapolation method for the prediction of the configuration performance.
   * @param extrapolationMethodOption
   */
  void setExtrapolationMethodOption(ExtrapolationMethodOption extrapolationMethodOption) {
    AutoPas::_extrapolationMethodOption = extrapolationMethodOption;
  }

  /**
   * Get the selector configuration strategy.
   * @return
   */
  [[nodiscard]] SelectorStrategyOption getSelectorStrategy() const { return _selectorStrategy; }

  /**
   * Set the strategy of how to select a performance value for a piece of evidence from multiple time measurements
   * (=samples). For possible selector strategy choices see options::SelectorStrategyOption::Value.
   * @param selectorStrategy
   */
  void setSelectorStrategy(SelectorStrategyOption selectorStrategy) { AutoPas::_selectorStrategy = selectorStrategy; }

  /**
   * Get the list of allowed load estimation algorithms.
   * @return
   */
  const std::set<LoadEstimatorOption> &getAllowedLoadEstimators() const { return _allowedLoadEstimators; }

  /**
   * Set the list of allowed load estimation algorithms.
   * For possible container choices see AutoPas::LoadEstimatorOption.
   * @param allowedLoadEstimators
   */
  void setAllowedLoadEstimators(const std::set<LoadEstimatorOption> &allowedLoadEstimators) {
    AutoPas::_allowedLoadEstimators = allowedLoadEstimators;
  }

  /**
   * Get the list of allowed containers.
   * @return
   */
  [[nodiscard]] const std::set<ContainerOption> &getAllowedContainers() const { return _allowedContainers; }

  /**
   * Set the list of allowed containers.
   * For possible container choices see options::ContainerOption::Value.
   * @param allowedContainers
   */
  void setAllowedContainers(const std::set<ContainerOption> &allowedContainers) {
    AutoPas::_allowedContainers = allowedContainers;
  }

  /**
   * Get the list of allowed traversals.
   * @return
   */
  [[nodiscard]] const std::set<TraversalOption> &getAllowedTraversals() const { return _allowedTraversals; }

  /**
   * Set the list of allowed traversals.
   * For possible traversals choices see options::TraversalOption::Value.
   * @param allowedTraversals
   */
  void setAllowedTraversals(const std::set<TraversalOption> &allowedTraversals) {
    AutoPas::_allowedTraversals = allowedTraversals;
  }

  /**
   * Get the list of allowed data layouts.
   * @return
   */
  [[nodiscard]] const std::set<DataLayoutOption> &getAllowedDataLayouts() const { return _allowedDataLayouts; }

  /**
   * Set the list of allowed data layouts.
   * For possible data layout choices see options::DataLayoutOption::Value.
   * @param allowedDataLayouts
   */
  void setAllowedDataLayouts(const std::set<DataLayoutOption> &allowedDataLayouts) {
    AutoPas::_allowedDataLayouts = allowedDataLayouts;
  }

  /**
   * Get the list of allowed newton 3 options.
   * @return
   */
  [[nodiscard]] const std::set<Newton3Option> &getAllowedNewton3Options() const { return _allowedNewton3Options; }

  /**
   * Set the list of allowed newton 3 options.
   * For possible newton 3 choices see options::Newton3Option::Value.
   * @param allowedNewton3Options
   */
  void setAllowedNewton3Options(const std::set<Newton3Option> &allowedNewton3Options) {
    AutoPas::_allowedNewton3Options = allowedNewton3Options;
  }

  /**
   * Getter for the currently selected configuration.
   * @return Configuration object currently used.
   */
  [[nodiscard]] const Configuration &getCurrentConfig() const { return _autoTuner->getCurrentConfig(); }

  /**
   * Getter for the tuning strategy option.
   * @return
   */
  [[nodiscard]] TuningStrategyOption getTuningStrategyOption() const { return _tuningStrategyOption; }

  /**
   * Setter for the tuning strategy option.
   * For possible tuning strategy choices see options::TuningStrategyOption::Value.
   * @param tuningStrategyOption
   */
  void setTuningStrategyOption(TuningStrategyOption tuningStrategyOption) {
    _tuningStrategyOption = tuningStrategyOption;
  }

  /**
   * Setter for the mpi strategy option
   * @param mpiStrategyOption
   */
  void setMPIStrategy(MPIStrategyOption mpiStrategyOption) { _mpiStrategyOption = mpiStrategyOption; }

// Only define the interface for the MPI communicator if AUTOPAS_INCLUDE_MPI=ON
// The internal implementation will use _autopasMPICommunicator with WrapMPI regardless of AUTOPAS_INCLUDE_MPI
#if defined(AUTOPAS_INCLUDE_MPI)
  /**
   * Setter for the MPI communicator that AutoPas uses for potential MPI calls.
   * If not set, MPI_COMM_WORLD will be used.
   * @param comm: communicator (handle)
   */
  void setMPICommunicator(MPI_Comm comm) { _autopasMPICommunicator = comm; }

  /**
   * Getter for the AutoPas MPI communicator
   * @return communicator
   */
  MPI_Comm getMPICommunicator() { return _autopasMPICommunicator; }
#endif

  /**
   * Suffix for all output files produced by this instance of AutoPas, e.g. from csv loggers.
   * This is useful when multiple instances of AutoPas exist, especially in an MPI context.
   * @param suffix
   */
  void setOutputSuffix(const std::string &suffix) { _outputSuffix = suffix; }

 private:
  /**
   * Lower corner of the container.
   */
  std::array<double, 3> _boxMin{0, 0, 0};
  /**
   * Upper corner of the container.
   */
  std::array<double, 3> _boxMax{0, 0, 0};
  /**
   * Cutoff radius to be used in this container.
   */
  double _cutoff{1.0};
  /**
   * Length added to the cutoff for the Verlet lists' skin.
   */
  double _verletSkin{0.2};
  /**
   * Specifies after how many pair-wise traversals the neighbor lists are to be rebuild.
   */
  unsigned int _verletRebuildFrequency{20};
  /**
   * Specifies the size of clusters for Verlet lists.
   */
  unsigned int _verletClusterSize{4};
  /**
   * Number of timesteps after which the auto-tuner shall reevaluate all selections.
   */
  unsigned int _tuningInterval{5000};
  /**
   * Number of samples the tuner should collect for each combination.
   */
  unsigned int _numSamples{3};
  /**
   * Tuning Strategies which work on a fixed number of evidence should use this value.
   */
  unsigned int _maxEvidence{10};
  /**
   * Specifies the factor of the range of the optimal configurations in PredicitveTuning.
   */
  double _relativeOptimumRange{1.2};
  /**
   * Specifies how many tuning phases a configuration can not be tested in PredicitveTuning.
   */
  unsigned int _maxTuningPhasesWithoutTest{5};
  /**
   * The relative cutoff for configurations to be blacklisted.
   * E.g. 2.5 means all configurations that take 2.5x the time of the optimum are blacklisted.
   */
  double _relativeBlacklistRange{0};
  /**
   * Specifies how many tests that need to have happened for a configuration until the first prediction is calculated in
   * PredictiveTuning.
   */
  unsigned int _evidenceFirstPrediction{3};
  /**
   * Acquisition function used for tuning.
   * For possible acquisition function choices see options::AcquisitionFunction::Value.
   */
  AcquisitionFunctionOption _acquisitionFunctionOption{AcquisitionFunctionOption::upperConfidenceBound};

  /**
   * Extrapolation method used in predictiveTuning.
   * For possible extrapolation method choices see autopas/options/ExtrapolationMethodOption.
   */
  ExtrapolationMethodOption _extrapolationMethodOption{ExtrapolationMethodOption::linearRegression};

  /**
   * Strategy option for the auto tuner.
   * For possible tuning strategy choices see options::TuningStrategyOption::Value.
   */
  TuningStrategyOption _tuningStrategyOption{TuningStrategyOption::fullSearch};

  /**
   * Strategy for the configuration selector.
   * For possible selector strategies see options::SelectorStrategyOption::Value.
   */
  SelectorStrategyOption _selectorStrategy{SelectorStrategyOption::fastestAbs};

  /**
   * List of container types AutoPas can choose from.
   * For possible container choices see options::ContainerOption::Value.
   */
  std::set<ContainerOption> _allowedContainers{ContainerOption::getMostOptions()};

  /**
   * List of traversals AutoPas can choose from.
   * For possible traversal choices see options::TraversalOption::Value.
   */
  std::set<TraversalOption> _allowedTraversals{TraversalOption::getMostOptions()};

  /**
   * List of data layouts AutoPas can choose from.
   * For possible data layout choices see options::DataLayoutOption::Value.
   */
  std::set<DataLayoutOption> _allowedDataLayouts{DataLayoutOption::getMostOptions()};

  /**
   * Whether AutoPas is allowed to exploit Newton's third law of motion.
   */
  std::set<Newton3Option> _allowedNewton3Options{Newton3Option::getMostOptions()};

  /**
   * Whether the chosen tuning strategy will be parallelized by MPI
   */
  MPIStrategyOption _mpiStrategyOption{MPIStrategyOption::noMPI};

  /**
   * Cell size factor to be used in this container (only relevant for LinkedCells, VerletLists and VerletListsCells).
   */
  std::unique_ptr<NumberSet<double>> _allowedCellSizeFactors{
      std::make_unique<NumberSetFinite<double>>(std::set<double>({1.}))};

  /***
   * Load estimation algorithm to be used for efficient parallelisation (only relevant for LCSlicedBalancedTraversal and
   * VLCSlicedBalancedTraversal).
   */
  std::set<LoadEstimatorOption> _allowedLoadEstimators{LoadEstimatorOption::getAllOptions()};

  /**
   * LogicHandler of autopas.
   */
  std::unique_ptr<autopas::LogicHandler<Particle>> _logicHandler;

  /**
   * This is the AutoTuner that owns the container, ...
   */
  std::unique_ptr<autopas::AutoTuner<Particle>> _autoTuner;

  /**
   * Communicator that should be used for MPI calls inside of AutoPas
   */
  AutoPas_MPI_Comm _autopasMPICommunicator{AUTOPAS_MPI_COMM_NULL};

  /**
   * Stores whether the mpi communicator was provided externally or not
   */
  bool _externalMPICommunicator{false};

  /**
   * Suffix for all output files produced by this instance of AutoPas, e.g. from csv loggers.
   * This is useful when multiple instances of AutoPas exist, especially in an MPI context.
   */
  std::string _outputSuffix{""};

};  // class AutoPas
}  // namespace autopas
