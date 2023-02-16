/**
 * @file AutoPasDecl.h
 *
 */

#pragma once

#include <memory>
#include <set>

#include "autopas/containers/ParticleContainerInterface.h"
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
#include "autopas/utils/StaticContainerSelector.h"
#include "autopas/utils/WrapMPI.h"

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
   * Define the iterator type for simple use, also from the outside.
   * Helps to, e.g., wrap the AutoPas iterators
   */
  using IteratorT = autopas::ContainerIterator<Particle, true, false>;

  /**
   * Define the const iterator type for simple use, also from the outside.
   * Helps to, e.g., wrap the AutoPas iterators
   */
  using ConstIteratorT = autopas::ContainerIterator<Particle, false, false>;

  /**
   * Define the region iterator type for simple use, also from the outside.
   * Helps to, e.g., wrap the AutoPas iterators
   */
  using RegionIteratorT = autopas::ContainerIterator<Particle, true, true>;

  /**
   * Define the const region iterator type for simple use, also from the outside.
   * Helps to, e.g., wrap the AutoPas iterators
   */
  using RegionConstIteratorT = autopas::ContainerIterator<Particle, false, true>;

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
   * Updates the container.
   * On an update, halo particles are deleted and particles that do no longer belong into the container will be removed
   * and returned.
   * @return A vector of invalid particles that do no longer belong in the current container.
   */
  [[nodiscard]] std::vector<Particle> updateContainer();

  /**
   * Adds a particle to the container.
   * This is only allowed if the neighbor lists are not valid.
   * @param p Reference to the particle to be added
   * @note An exception is thrown if the particle is added and it is not inside of the owned domain (defined by
   * boxmin and boxmax) of the container.
   * @note This function is NOT thread-safe.
   */
  void addParticle(const Particle &p);

  /**
   * Adds a particle to the container that lies in the halo region of the container.
   * @param haloParticle Particle to be added.
   * @note An exception is thrown if the halo particle is added and it is inside of the owned domain (defined by boxmin
   * and boxmax) of the container.
   * @note This function is NOT thread-safe.
   */
  void addHaloParticle(const Particle &haloParticle);

  /**
   * Deletes all particles.
   * @note This invalidates the container, a rebuild is forced on the next iteratePairwise() call.
   */
  void deleteAllParticles();

  /**
   * Deletes the particle behind the current iterator position and leaves the container in a valid state.
   *
   * Internally, depending on the container, this might just mark the particle as deleted without actually removing it.
   * If this can not be done without compromising e.g. a VerletList reference structure the particle is only marked.
   *
   * @param iter Needs to be a modify-able iterator.
   */
  void deleteParticle(ContainerIterator<Particle, true, false> &iter);

  /**
   * @copydoc deleteParticle(ContainerIterator<Particle, true, false> &iter)
   *
   * Region Iterator version.
   */
  void deleteParticle(ContainerIterator<Particle, true, true> &iter);

  /**
   * Deletes the given particle and leaves the container in a valid state.
   *
   * Internally, depending on the container, this might just mark the particle as deleted without actually removing it.
   * If this can not be done without compromising e.g. a VerletList reference structure the particle is only marked.
   *
   * @note This function might invalidate iterators.
   *
   * @param particle Reference to the particle that should be deleted.
   *
   * @return True iff the reference still points to a valid particle.
   */
  bool deleteParticle(Particle &particle);

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
  IteratorT begin(IteratorBehavior behavior = IteratorBehavior::ownedOrHalo);

  /**
   * @copydoc begin()
   * @note const version
   */
  ConstIteratorT begin(IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const;

  /**
   * execute code on all particles in parallel as defined by a lambda function
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   * @note not actually parallel until kokkos integration
   */
  template <typename Lambda>
  void forEachParallel(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto containerPtr) { containerPtr->forEach(forEachLambda, behavior); });
  }

  /**
   * @copydoc forEachParallel()
   * @note const version
   */
  template <typename Lambda>
  void forEachParallel(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto containerPtr) { containerPtr->forEach(forEachLambda, behavior); });
  }

  /**
   * Execute code on all particles as defined by a lambda function.
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    withStaticContainerType(getContainer(), [&](auto containerPtr) { containerPtr->forEach(forEachLambda, behavior); });
  }

  /**
   * @copydoc forEach()
   * @note const version
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    withStaticContainerType(getContainer(), [&](auto containerPtr) { containerPtr->forEach(forEachLambda, behavior); });
  }

  /**
   * Reduce properties of particles in parallel as defined by a lambda function.
   * @tparam Lambda (Particle p, A &initialValue) -> void
   * @tparam reference to result of type A
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownedOrHalo
   * @note not actually parallel until kokkos integration
   */
  template <typename Lambda, typename A>
  void reduceParallel(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(),
                            [&](auto containerPtr) { containerPtr->reduce(reduceLambda, result, behavior); });
  }

  /**
   * @copydoc reduce()
   * @note const version
   */
  template <typename Lambda, typename A>
  void reduceParallel(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(),
                            [&](auto containerPtr) { containerPtr->reduce(reduceLambda, result, behavior); });
  }

  /**
   * Reduce properties of particles as defined by a lambda function.
   * @tparam Lambda (Particle p, A &initialValue) -> void
   * @tparam reference to result of type A
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownedOrHalo
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    withStaticContainerType(getContainer(),
                            [&](auto containerPtr) { containerPtr->reduce(reduceLambda, result, behavior); });
  }

  /**
   * @copydoc reduce()
   * @note const version
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    withStaticContainerType(getContainer(),
                            [&](auto containerPtr) { containerPtr->reduce(reduceLambda, result, behavior); });
  }

  /**
   * @copydoc begin()
   * @note cbegin will guarantee to return a const_iterator.
   */
  ConstIteratorT cbegin(IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const { return begin(behavior); }

  /**
   * Helper to enable range-based for loops for the AutoPas object.
   * ParticleIterator::operator==() compares its own validity state against this value. Hence, as soon as the iterator
   * is invalid the loop ends.
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
  RegionIteratorT getRegionIterator(const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
                                    IteratorBehavior behavior = IteratorBehavior::ownedOrHalo);
  /**
   * @copydoc getRegionIterator()
   * @note const version
   */
  RegionConstIteratorT getRegionIterator(const std::array<double, 3> &lowerCorner,
                                         const std::array<double, 3> &higherCorner,
                                         IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const;

  /**
   * Execute code on all particles in a certain region in parallel as defined by a lambda function.
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
    // TODO (lgaertner): parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto containerPtr) {
      containerPtr->forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * @copydoc forEachInRegionParallel()
   * @note const version
   */
  template <typename Lambda>
  void forEachInRegionParallel(Lambda forEachLambda, std::array<double, 3> lowerCorner,
                               std::array<double, 3> higherCorner,
                               IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    // TODO (lgaertner): parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto containerPtr) {
      containerPtr->forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * Execute code on all particles in a certain region as defined by a lambda function.
   * @tparam Lambda (Particle &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
                       IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    withStaticContainerType(getContainer(), [&](auto containerPtr) {
      containerPtr->forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * @copydoc forEachInRegion()
   * @note const version
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
                       IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    withStaticContainerType(getContainer(), [&](auto containerPtr) {
      containerPtr->forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * Execute code on all particles in a certain region in parallel as defined by a lambda function.
   * @tparam Lambda (Particle &p, A &result) -> void
   * @tparam A type of reduction value
   * @param reduceLambda code to be executed on all particles
   * @param result reference to starting and final value of reduction
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   * @note not actually parallel until kokkos integration
   */
  template <typename Lambda, typename A>
  void reduceInRegionParallel(Lambda reduceLambda, A &result, std::array<double, 3> lowerCorner,
                              std::array<double, 3> higherCorner,
                              IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto containerPtr) {
      containerPtr->reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * @copydoc reduceInRegion()
   * @note const version
   */
  template <typename Lambda, typename A>
  void reduceInRegionParallel(Lambda reduceLambda, A &result, std::array<double, 3> lowerCorner,
                              std::array<double, 3> higherCorner,
                              IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto containerPtr) {
      containerPtr->reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * Execute code on all particles in a certain region as defined by a lambda function.
   * @tparam Lambda (Particle &p, A &result) -> void
   * @tparam A type of reduction value
   * @param reduceLambda code to be executed on all particles
   * @param result reference to starting and final value of reduction
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, std::array<double, 3> lowerCorner,
                      std::array<double, 3> higherCorner, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    withStaticContainerType(getContainer(), [&](auto containerPtr) {
      containerPtr->reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * @copydoc reduceInRegion()
   * @note const version
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, std::array<double, 3> lowerCorner,
                      std::array<double, 3> higherCorner,
                      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    withStaticContainerType(getContainer(), [&](auto containerPtr) {
      containerPtr->reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
    });
  }
  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @return _verletSkin
   */
  double getVerletSkin() {
    double _verletSkin = AutoPas::_verletSkinPerTimestep * AutoPas::_verletRebuildFrequency;
    return _verletSkin;
  };

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
   * get the bool value indicating if the search space is trivial (not more than one configuration to test).
   * @return bool indicating if search space is trivial.
   */
  [[nodiscard]] bool searchSpaceIsTrivial();

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
      AutoPasLog(ERROR, "Cutoff <= 0.0: {}", cutoff);
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
      AutoPasLog(ERROR, "cell size <= 0.0");
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
      AutoPasLog(ERROR, "cell size <= 0.0: {}", cellSizeFactor);
      utils::ExceptionHandler::exception("Error: cell size <= 0.0!");
    }
    AutoPas::_allowedCellSizeFactors = std::make_unique<NumberSetFinite<double>>(std::set<double>{cellSizeFactor});
  }

  /**
   * Get length added to the cutoff for the Verlet lists' skin per timestep.
   * @return _verletSkinPerTimestep
   */
  [[nodiscard]] double getVerletSkinPerTimestep() const { return _verletSkinPerTimestep; }

  /**
   * Set length added to the cutoff for the Verlet lists' skin per timestep.
   * @param verletSkinPerTimestep
   */
  void setVerletSkinPerTimestep(double verletSkinPerTimestep) {
    AutoPas::_verletSkinPerTimestep = verletSkinPerTimestep;
  }

  /**
   * Get Verlet rebuild frequency.
   * @return _verletRebuildFrequency
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

  /**
   * Setter for the maximal Difference for the bucket distribution
   * @param MPITuningMaxDifferenceForBucket
   */
  void setMPITuningMaxDifferenceForBucket(double MPITuningMaxDifferenceForBucket) {
    _mpiTuningMaxDifferenceForBucket = MPITuningMaxDifferenceForBucket;
  }

  /**
   * Setter for the maxDensity-Weight in calculation for bucket distribution
   * @param MPITuningWeightForMaxDensity
   */
  void setMPITuningWeightForMaxDensity(double MPITuningWeightForMaxDensity) {
    _mpiTuningWeightForMaxDensity = MPITuningWeightForMaxDensity;
  }

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
  std::shared_ptr<autopas::ParticleContainerInterface<Particle>> getContainer();

  std::shared_ptr<const autopas::ParticleContainerInterface<Particle>> getContainer() const;

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
   * Length added to the cutoff for the Verlet lists' skin per Timestep.
   */
  double _verletSkinPerTimestep{0.01};
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
   * For MPI-tuning: Maximum of the relative difference in the comparison metric for two ranks which exchange their
   * tuning information.
   */
  double _mpiTuningMaxDifferenceForBucket{0.3};

  /**
   * For MPI-tuning: Weight for maxDensity in the calculation for bucket distribution.
   */
  double _mpiTuningWeightForMaxDensity{0.0};

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
