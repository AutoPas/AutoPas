/**
 * @file AutoPasDecl.h
 *
 */

#pragma once

#include <memory>
#include <set>

#include "autopas/LogicHandlerInfo.h"
#include "autopas/containers/ParticleContainerInterface.h"
#include "autopas/options//ExtrapolationMethodOption.h"
#include "autopas/options/AcquisitionFunctionOption.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/EnergySensorOption.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/options/LoadEstimatorOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/TuningMetricOption.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/tuning/AutoTuner.h"
#include "autopas/tuning/Configuration.h"
#include "autopas/tuning/tuningStrategy/TuningStrategyFactoryInfo.h"
#include "autopas/utils/NumberSet.h"
#include "autopas/utils/StaticContainerSelector.h"
#include "autopas/utils/WrapMPI.h"

namespace autopas {

// Forward declare Handler so that including this header does not include the whole library with all
// containers and traversals.
template <class Particle_T>
class LogicHandler;

/**
 * The AutoPas class is intended to be the main point of Interaction for the user.
 * It acts as an interface from where all features of the library can be triggered and configured.
 * @tparam Particle_T Class for particles
 * @tparam ParticleCell Class for the particle cells
 */
template <class Particle_T>
class AutoPas {
 public:
  /**
   * Particle type to be accessible after initialization.
   */
  using ParticleType = Particle_T;

  /**
   * Define the iterator type for ease of use. Also for external use.
   * Helps to, e.g., wrap the AutoPas iterators
   */
  using IteratorT = autopas::ContainerIterator<Particle_T, true, false>;

  /**
   * Define the const iterator type for ease of use. Also for external use.
   * Helps to, e.g., wrap the AutoPas iterators
   */
  using ConstIteratorT = autopas::ContainerIterator<Particle_T, false, false>;

  /**
   * Define the region iterator type for ease of use. Also for external use.
   * Helps to, e.g., wrap the AutoPas iterators
   */
  using RegionIteratorT = autopas::ContainerIterator<Particle_T, true, true>;

  /**
   * Define the const region iterator type for ease of use. Also for external use.
   * Helps to, e.g., wrap the AutoPas iterators
   */
  using RegionConstIteratorT = autopas::ContainerIterator<Particle_T, false, true>;

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
  std::vector<Particle_T> resizeBox(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax);

  /**
   * Force the internal tuner to enter a new tuning phase upon the next call to computeInteractions().
   */
  void forceRetune();

  /**
   * Free the AutoPas MPI communicator.
   * To be called before MPI_Finalize.
   * If no MPI is used just call this at the end of the program.
   */
  void finalize() {
    if (not _externalMPICommunicator) {
      AutoPas_MPI_Comm_free(&_tuningStrategyFactoryInfo.autopasMpiCommunicator);
    }
  }

  /**
   * Updates the container.
   * On an update, halo particles are deleted and particles that do no longer belong into the container will be removed
   * and returned.
   * @return A vector of invalid particles that do no longer belong in the current container.
   */
  [[nodiscard]] std::vector<Particle_T> updateContainer();

  /**
   * Reserve memory for a given number of particles in the container and logic layers.
   * This function assumes a uniform distribution of particles throughout the domain.
   * For example, this means that in a LinkedCells Container in each cell vector.reserve(numParticles/numCells) is
   * called.
   * @note This functions will create an estimate for the number of halo particles.
   * @param numParticles No buffer factor is applied. It is probably wise to slightly over-reserve to account for
   * imbalance or particle movement.
   */
  void reserve(size_t numParticles);

  /**
   * Reserve memory for a given number of particles in the container and logic layers
   * (e.g. LogicHandler::_particleBuffer).
   * This function assumes a uniform distribution of particles throughout the domain.
   * For example, this means that in a LinkedCells Container in each cell vector.reserve(numParticles/numCells) is
   * called.
   * @param numParticles
   * @param numHaloParticles
   */
  void reserve(size_t numParticles, size_t numHaloParticles);

  /**
   * Adds a particle to the container.
   * This is only allowed if the neighbor lists are not valid.
   * @param p Reference to the particle to be added
   * @note An exception is thrown if the particle is added and it is not inside of the owned domain (defined by
   * boxMin and boxMax) of the container.
   * @note This function is NOT thread-safe if the container is Octree.
   */
  void addParticle(const Particle_T &p);

  /**
   * Adds all particles from the collection to the container.
   * @note This function uses reserve().
   * @note This function uses addParticle().
   * @tparam Collection Collection type that contains the particles (e.g. std::vector). Needs to support `.size()`.
   * @param particles
   */
  template <class Collection>
  void addParticles(Collection &&particles);

  /**
   * Adds all particles for which predicate(particle) == true to the container.
   * @note This function uses reserve().
   * @note This function uses addParticle().
   * @tparam Collection Collection type that contains the particles (e.g. std::vector). Needs to support `.size()`.
   * @tparam F Function type of predicate. Should be of the form: (const Particle_T &) -> bool.
   * @param particles Particles that are potentially added.
   * @param predicate Condition that determines if an individual particle should be added.
   */
  template <class Collection, class F>
  void addParticlesIf(Collection &&particles, F predicate);

  /**
   * Adds a particle to the container that lies in the halo region of the container.
   * @param haloParticle Particle to be added.
   * @note An exception is thrown if the halo particle is added and it is inside of the owned domain (defined by boxMin
   * and boxMax) of the container.
   * @note This function is NOT thread-safe if the container is Octree.
   */
  void addHaloParticle(const Particle_T &haloParticle);

  /**
   * Adds all halo particles from the collection to the container.
   * @note This function uses reserve().
   * @note This function uses addHaloParticle().
   * @tparam Collection Collection type that contains the particles (e.g. std::vector). Needs to support `.size()`.
   * @param particles
   */
  template <class Collection>
  void addHaloParticles(Collection &&particles);

  /**
   * Adds all halo particles for which predicate(particle) == true to the container.
   * @note This function uses reserve().
   * @note This function uses addHaloParticle().
   * @tparam Collection Collection type that contains the particles (e.g. std::vector). Needs to support `.size()`.
   * @tparam F Function type of predicate. Should be of the form: (const Particle_T &) -> bool.
   * @param particles Particles that are potentially added.
   * @param predicate Condition that determines if an individual particle should be added.
   */
  template <class Collection, class F>
  void addHaloParticlesIf(Collection &&particles, F predicate);

  /**
   * Deletes all particles.
   * @note This invalidates the container, a rebuild is forced on the next computeInteractions() call.
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
  void deleteParticle(IteratorT &iter);

  /**
   * @copydoc deleteParticle(IteratorT &iter)
   *
   * Region Iterator version.
   */
  void deleteParticle(RegionIteratorT &iter);

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
  bool deleteParticle(Particle_T &particle);

  /**
   * Function to iterate over all inter-particle interactions in the container
   * This function only handles short-range interactions.
   * @param f Functor that describes the interaction (e.g. force).
   * @return true if this was a tuning iteraction.
   */
  template <class Functor>
  bool computeInteractions(Functor *f);

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
   * @tparam Lambda (Particle_T &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   * @note not actually parallel until kokkos integration
   */
  template <typename Lambda>
  void forEachParallel(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto &container) { container.forEach(forEachLambda, behavior); });
  }

  /**
   * @copydoc forEachParallel()
   * @note const version
   */
  template <typename Lambda>
  void forEachParallel(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto &container) { container.forEach(forEachLambda, behavior); });
  }

  /**
   * Execute code on all particles as defined by a lambda function.
   * @tparam Lambda (Particle_T &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    withStaticContainerType(getContainer(), [&](auto &container) { container.forEach(forEachLambda, behavior); });
  }

  /**
   * @copydoc forEach()
   * @note const version
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    withStaticContainerType(getContainer(), [&](auto &container) { container.forEach(forEachLambda, behavior); });
  }

  /**
   * Reduce properties of particles in parallel as defined by a lambda function.
   * @tparam Lambda (Particle_T p, A &initialValue) -> void
   * @tparam reference to result of type A
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownedOrHalo
   * @note not actually parallel until kokkos integration
   */
  template <typename Lambda, typename A>
  void reduceParallel(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto &container) { container.reduce(reduceLambda, result, behavior); });
  }

  /**
   * @copydoc reduce()
   * @note const version
   */
  template <typename Lambda, typename A>
  void reduceParallel(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto &container) { container.reduce(reduceLambda, result, behavior); });
  }

  /**
   * Reduce properties of particles as defined by a lambda function.
   * @tparam Lambda (Particle_T p, A &initialValue) -> void
   * @tparam reference to result of type A
   * @param reduceLambda code to reduce properties of particles
   * @param result reference to result of type A
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownedOrHalo
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    withStaticContainerType(getContainer(), [&](auto &container) { container.reduce(reduceLambda, result, behavior); });
  }

  /**
   * @copydoc reduce()
   * @note const version
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    withStaticContainerType(getContainer(), [&](auto &container) { container.reduce(reduceLambda, result, behavior); });
  }

  /**
   * @copydoc begin()
   * @note cbegin will guarantee to return a const_iterator.
   */
  ConstIteratorT cbegin(IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const { return begin(behavior); }

  /**
   * Dummy to make range-based for loops work.
   *
   * Range-Based for loops use the incremented begin() expression and compare it against the end() expression.
   * ContainerIterator implements ContainerIterator::operator==() that accepts a bool as right hand side argument,
   * which is triggered by this end() function.
   * This operator then proceeds to check the validity of the iterator itself.
   *
   * @return false
   */
  [[nodiscard]] constexpr bool end() const { return false; }

  /**
   * Iterate over all particles in a specified region
   * ```c++
   * for (auto iter = container.getRegionIterator(lowCorner, highCorner); iter.isValid(); ++iter) { }
   * ```
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
   * @tparam Lambda (Particle_T &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   * @note not actually parallel until kokkos integration
   */
  template <typename Lambda>
  void forEachInRegionParallel(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                               const std::array<double, 3> &higherCorner,
                               IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    // TODO (lgaertner): parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto &container) {
      container.forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * @copydoc forEachInRegionParallel()
   * @note const version
   */
  template <typename Lambda>
  void forEachInRegionParallel(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                               const std::array<double, 3> &higherCorner,
                               IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    // TODO (lgaertner): parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto &container) {
      container.forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * Execute code on all particles in a certain region as defined by a lambda function.
   * @tparam Lambda (Particle_T &p) -> void
   * @param forEachLambda code to be executed on all particles
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner,
                       IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    withStaticContainerType(getContainer(), [&](auto &container) {
      container.forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * @copydoc forEachInRegion()
   * @note const version
   */
  template <typename Lambda>
  void forEachInRegion(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
                       const std::array<double, 3> &higherCorner,
                       IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    withStaticContainerType(getContainer(), [&](auto &container) {
      container.forEachInRegion(forEachLambda, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * Execute code on all particles in a certain region in parallel as defined by a lambda function.
   * @tparam Lambda (Particle_T &p, A &result) -> void
   * @tparam A type of reduction value
   * @param reduceLambda code to be executed on all particles
   * @param result reference to starting and final value of reduction
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   * @note not actually parallel until kokkos integration
   */
  template <typename Lambda, typename A>
  void reduceInRegionParallel(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                              const std::array<double, 3> &higherCorner,
                              IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto &container) {
      container.reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * @copydoc reduceInRegion()
   * @note const version
   */
  template <typename Lambda, typename A>
  void reduceInRegionParallel(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                              const std::array<double, 3> &higherCorner,
                              IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    // TODO lgaertner: parallelize with kokkos integration
    withStaticContainerType(getContainer(), [&](auto &container) {
      container.reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * Execute code on all particles in a certain region as defined by a lambda function.
   * @tparam Lambda (Particle_T &p, A &result) -> void
   * @tparam A type of reduction value
   * @param reduceLambda code to be executed on all particles
   * @param result reference to starting and final value of reduction
   * @param lowerCorner lower corner of bounding box
   * @param higherCorner higher corner of bounding box
   * @param behavior @see IteratorBehavior default: @see IteratorBehavior::ownerOrHalo
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner,
                      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) {
    withStaticContainerType(getContainer(), [&](auto &container) {
      container.reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * @copydoc reduceInRegion()
   * @note const version
   */
  template <typename Lambda, typename A>
  void reduceInRegion(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
                      const std::array<double, 3> &higherCorner,
                      IteratorBehavior behavior = IteratorBehavior::ownedOrHalo) const {
    withStaticContainerType(getContainer(), [&](auto &container) {
      container.reduceInRegion(reduceLambda, result, lowerCorner, higherCorner, behavior);
    });
  }

  /**
   * Function to iterate over all pairs of particles in the container.
   * This function only handles short-range interactions.
   * @return _verletSkin
   */
  double getVerletSkin() { return _logicHandlerInfo.verletSkin; };

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
  [[nodiscard]] const std::array<double, 3> &getBoxMin() const;

  /**
   * Get the upper corner of the container without the halo.
   * @return upper corner of the container.
   */
  [[nodiscard]] const std::array<double, 3> &getBoxMax() const;

  /**
   * get the bool value indicating if the search space is trivial (not more than one configuration to test).
   * @return bool indicating if search space is trivial.
   */
  [[nodiscard]] bool searchSpaceIsTrivial();

  /**
   * Set coordinates of the lower corner of the domain.
   * @param boxMin
   */
  void setBoxMin(const std::array<double, 3> &boxMin) { _logicHandlerInfo.boxMin = boxMin; }

  /**
   * Set coordinates of the upper corner of the domain.
   * @param boxMax
   */
  void setBoxMax(const std::array<double, 3> &boxMax) { _logicHandlerInfo.boxMax = boxMax; }

  /**
   * Get cutoff radius.
   * @return
   */
  [[nodiscard]] double getCutoff() const { return _logicHandlerInfo.cutoff; }

  /**
   * Set cutoff radius.
   * @param cutoff
   */
  void setCutoff(double cutoff) {
    if (cutoff <= 0.0) {
      utils::ExceptionHandler::exception("Error: Cutoff has to be positive {} <= 0.0!", cutoff);
    }
    _logicHandlerInfo.cutoff = cutoff;
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
      utils::ExceptionHandler::exception("Error: minimum cell size factor has to be positive {} <= 0.0!",
                                         allowedCellSizeFactors.getMin());
    }
    _allowedCellSizeFactors = std::move(allowedCellSizeFactors.clone());
  }

  /**
   * Set allowed cell size factors to one element (only relevant for LinkedCells, VerletLists and VerletListsCells).
   * @param cellSizeFactor
   */
  void setCellSizeFactor(double cellSizeFactor) {
    if (cellSizeFactor <= 0.0) {
      utils::ExceptionHandler::exception("Error: cell size factor has to be positive! {}<= 0.0!", cellSizeFactor);
    }
    _allowedCellSizeFactors = std::make_unique<NumberSetFinite<double>>(std::set<double>{cellSizeFactor});
  }

  /**
   * Set length added to the cutoff for the Verlet lists' skin per timestep.
   * @param verletSkin
   */
  void setVerletSkin(double verletSkin) { _logicHandlerInfo.verletSkin = verletSkin; }

  /**
   * Set time step of the simulation.
   * This is currently used for estimating the rebuild frequency.
   * @param deltaT
   */
  void setDeltaT(double deltaT) { _logicHandlerInfo.deltaT = deltaT; }

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
    _verletRebuildFrequency = verletRebuildFrequency;
  }
  /**
   * Get Verlet cluster size.
   * @return
   */
  [[nodiscard]] unsigned int getVerletClusterSize() const { return _logicHandlerInfo.verletClusterSize; }

  /**
   * Set Verlet cluster size.
   * @param verletClusterSize
   */
  void setVerletClusterSize(unsigned int verletClusterSize) { _logicHandlerInfo.verletClusterSize = verletClusterSize; }

  /**
   * Get tuning interval.
   * @return
   */
  [[nodiscard]] unsigned int getTuningInterval() const { return _autoTunerInfo.tuningInterval; }

  /**
   * Set tuning interval.
   * @param tuningInterval
   */
  void setTuningInterval(unsigned int tuningInterval) { _autoTunerInfo.tuningInterval = tuningInterval; }

  /**
   * Get number of samples taken per configuration during the tuning.
   * @return
   */
  [[nodiscard]] unsigned int getNumSamples() const { return _autoTunerInfo.maxSamples; }

  /**
   * Set number of samples taken per configuration during the tuning.
   * @param numSamples
   */
  void setNumSamples(unsigned int numSamples) { _autoTunerInfo.maxSamples = numSamples; }

  /**
   * Set the earlyStoppingFactor for the auto tuner. If a configuration seems to be slower than the optimum
   * configuration found so far by more than this factor, it will not be sampled again during that tuning phase.
   * @param earlyStoppingFactor
   */
  void setEarlyStoppingFactor(double earlyStoppingFactor) { _autoTunerInfo.earlyStoppingFactor = earlyStoppingFactor; }

  /**
   * Get flag for whether a LOESS-based smoothening is used.
   * @return
   */
  [[nodiscard]] bool getUseLOESSSmoothening() const { return _autoTunerInfo.useLOESSSmoothening; }

  /**
   * Set flag for whether a LOESS-based smoothening is used.
   * @param useLOESSSmoothening
   */
  void setUseLOESSSmoothening(bool useLOESSSmoothening) { _autoTunerInfo.useLOESSSmoothening = useLOESSSmoothening; }

  /**
   * Get maximum number of evidence for tuning
   * @return
   */
  [[nodiscard]] unsigned int getMaxEvidence() const { return _tuningStrategyFactoryInfo.maxEvidence; }

  /**
   * Set maximum number of evidence for tuning
   * @param maxEvidence
   */
  void setMaxEvidence(unsigned int maxEvidence) { _tuningStrategyFactoryInfo.maxEvidence = maxEvidence; }

  /**
   * Get the range for the optimum in which has to be to be tested
   * @return
   */
  [[nodiscard]] double getRelativeOptimumRange() const { return _tuningStrategyFactoryInfo.relativeOptimumRange; }

  /**
   * Set the range for the optimum in which has to be to be tested
   * @param relativeOptimumRange
   */
  void setRelativeOptimumRange(double relativeOptimumRange) {
    _tuningStrategyFactoryInfo.relativeOptimumRange = relativeOptimumRange;
  }

  /**
   * Get the maximum number of tuning phases before a configuration is certainly tested again.
   * @return
   */
  [[nodiscard]] unsigned int getMaxTuningPhasesWithoutTest() const {
    return _tuningStrategyFactoryInfo.maxTuningPhasesWithoutTest;
  }

  /**
   * Set the maximum number of tuning phases before a configuration is certainly tested again.
   * @param maxTuningPhasesWithoutTest
   */
  void setMaxTuningPhasesWithoutTest(unsigned int maxTuningPhasesWithoutTest) {
    _tuningStrategyFactoryInfo.maxTuningPhasesWithoutTest = maxTuningPhasesWithoutTest;
  }

  /**
   * For Predictive tuning: Get the relative cutoff for configurations to be blacklisted.
   * E.g. 2.5 means all configurations that take 2.5x the time of the optimum are blacklisted.
   * @return
   */
  [[nodiscard]] double getRelativeBlacklistRange() const { return _tuningStrategyFactoryInfo.relativeBlacklistRange; }

  /**
   * Set the range of the configurations that are not going to be blacklisted.
   * @param relativeBlacklistRange
   */
  void setRelativeBlacklistRange(double relativeBlacklistRange) {
    _tuningStrategyFactoryInfo.relativeBlacklistRange = relativeBlacklistRange;
  }

  /**
   * Get the number of tests that need to have happened for a configuration until the first predictions are going to be
   * calculated.
   * @return
   */
  [[nodiscard]] unsigned int getEvidenceFirstPrediction() const {
    return _tuningStrategyFactoryInfo.minNumberOfEvidence;
  }

  /**
   * Set the number of tests that need to have happened for a configuration until the first predictions are going to be
   * calculated.
   * @param evidenceFirstPrediction
   */
  void setEvidenceFirstPrediction(unsigned int evidenceFirstPrediction) {
    _tuningStrategyFactoryInfo.minNumberOfEvidence = evidenceFirstPrediction;
  }

  /**
   * Get acquisition function used for tuning
   * @return
   */
  [[nodiscard]] AcquisitionFunctionOption getAcquisitionFunction() const {
    return _tuningStrategyFactoryInfo.acquisitionFunctionOption;
  }

  /**
   * Set acquisition function for tuning.
   * For possible acquisition function choices see options::AcquisitionFunctionOption::Value.
   * @note This function is only relevant for the bayesian based searches.
   * @param acqFun acquisition function
   */
  void setAcquisitionFunction(AcquisitionFunctionOption acqFun) {
    _tuningStrategyFactoryInfo.acquisitionFunctionOption = acqFun;
  }

  /**
   * Get extrapolation method for the prediction of the configuration performance.
   * @return
   */
  ExtrapolationMethodOption getExtrapolationMethodOption() const {
    return _tuningStrategyFactoryInfo.extrapolationMethodOption;
  }

  /**
   * Set extrapolation method for the prediction of the configuration performance.
   * @param extrapolationMethodOption
   */
  void setExtrapolationMethodOption(ExtrapolationMethodOption extrapolationMethodOption) {
    _tuningStrategyFactoryInfo.extrapolationMethodOption = extrapolationMethodOption;
  }

  /**
   * Get the selector configuration strategy.
   * @return
   */
  [[nodiscard]] SelectorStrategyOption getSelectorStrategy() const { return _autoTunerInfo.selectorStrategy; }

  /**
   * Set the strategy of how to select a performance value for a piece of evidence from multiple time measurements
   * (=samples). For possible selector strategy choices see options::SelectorStrategyOption::Value.
   * @param selectorStrategy
   */
  void setSelectorStrategy(SelectorStrategyOption selectorStrategy) {
    _autoTunerInfo.selectorStrategy = selectorStrategy;
  }

  /**
   * Get the list of allowed load estimation algorithms.
   * @return
   */
  const std::set<LoadEstimatorOption> &getAllowedLoadEstimators() const { return _allowedLoadEstimators; }

  /**
   * Set the list of allowed load estimation algorithms.
   * For possible container choices see LoadEstimatorOption.
   * @param allowedLoadEstimators
   */
  void setAllowedLoadEstimators(const std::set<LoadEstimatorOption> &allowedLoadEstimators) {
    _allowedLoadEstimators = allowedLoadEstimators;
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
    _allowedContainers = allowedContainers;
  }

  /**
   * Get the list of allowed traversals.
   * @param interactionType Get allowed traversals for this interaction type. Defaults to
   * InteractionTypeOption::pairwise.
   * @return
   */
  [[nodiscard]] const std::set<TraversalOption> &getAllowedTraversals(
      const InteractionTypeOption interactionType = InteractionTypeOption::pairwise) const {
    return _allowedTraversals.at(interactionType);
  }

  /**
   * Set the list of allowed traversals.
   * For possible traversals choices see options::TraversalOption::Value.
   * @param allowedTraversals
   * @param interactionType Set allowed traversals for this interaction type. Defaults to
   * InteractionTypeOption::pairwise.
   */
  void setAllowedTraversals(const std::set<TraversalOption> &allowedTraversals,
                            const InteractionTypeOption interactionType = InteractionTypeOption::pairwise) {
    if (interactionType == InteractionTypeOption::all) {
      for (auto iType : InteractionTypeOption::getMostOptions()) {
        _allowedTraversals[iType] = allowedTraversals;
      }
    } else {
      _allowedTraversals[interactionType] = allowedTraversals;
    }
  }

  /**
   * Get the list of allowed data layouts.
   * @return
   * @param interactionType Get allowed data layouts for this interaction type. Defaults to
   * InteractionTypeOption::pairwise.
   */
  [[nodiscard]] const std::set<DataLayoutOption> &getAllowedDataLayouts(
      const InteractionTypeOption interactionType = InteractionTypeOption::pairwise) const {
    return _allowedDataLayouts.at(interactionType);
  }

  /**
   * Set the list of allowed data layouts.
   * For possible data layouts choices see options::DataLayoutOption::Value.
   * @param allowedDataLayouts
   * @param interactionType Set allowed data layouts for this interaction type. Defaults to
   * InteractionTypeOption::pairwise.
   */
  void setAllowedDataLayouts(const std::set<DataLayoutOption> &allowedDataLayouts,
                             const InteractionTypeOption interactionType = InteractionTypeOption::pairwise) {
    if (interactionType == InteractionTypeOption::all) {
      for (auto iType : InteractionTypeOption::getMostOptions()) {
        _allowedDataLayouts[iType] = allowedDataLayouts;
      }
    } else {
      _allowedDataLayouts[interactionType] = allowedDataLayouts;
    }
  }

  /**
   * Get the list of allowed newton 3 options.
   * @param interactionType Get allowed newton 3 options for this interaction type. Defaults to
   * InteractionTypeOption::pairwise.
   * @return
   */
  [[nodiscard]] const std::set<Newton3Option> &getAllowedNewton3Options(
      const InteractionTypeOption interactionType = InteractionTypeOption::pairwise) const {
    return _allowedNewton3Options.at(interactionType);
  }

  /**
   * Set the list of allowed newton 3 options.
   * For possible newton 3 choices see options::Newton3Option::Value.
   * @param allowedNewton3Options
   * @param interactionType Set allowed newton 3 options for this interaction type. Defaults to
   * InteractionTypeOption::pairwise.
   */
  void setAllowedNewton3Options(const std::set<Newton3Option> &allowedNewton3Options,
                                const InteractionTypeOption interactionType = InteractionTypeOption::pairwise) {
    if (interactionType == InteractionTypeOption::all) {
      for (auto iType : InteractionTypeOption::getMostOptions()) {
        _allowedNewton3Options[iType] = allowedNewton3Options;
      }
    } else {
      _allowedNewton3Options[interactionType] = allowedNewton3Options;
    }
  }

  /**
   * Set the list of allowed interaction types.
   * AutoPas will initialize AutoTuners for the allowed interaction types.
   * For possible newton 3 choices see options::interactionTypeOption::Value.
   * @param allowedInteractionTypeOptions
   */
  void setAllowedInteractionTypeOptions(const std::set<InteractionTypeOption> &allowedInteractionTypeOptions) {
    _allowedInteractionTypeOptions = allowedInteractionTypeOptions;
  }

  /**
   * Getter for the currently selected configuration.
   * @return Configuration objects currently used for respective interaction types.
   */
  [[nodiscard]] std::unordered_map<InteractionTypeOption::Value, std::reference_wrapper<const Configuration>>
  getCurrentConfigs() const {
    std::unordered_map<InteractionTypeOption::Value, std::reference_wrapper<const Configuration>> currentConfigs;
    currentConfigs.reserve(_autoTuners.size());

    for (const auto &[type, tuner] : _autoTuners) {
      currentConfigs.emplace(type, std::cref(tuner->getCurrentConfig()));
    }
    return currentConfigs;
  }

  /**
   * Getter for the tuning strategy option.
   * @return
   */
  [[nodiscard]] const std::vector<TuningStrategyOption> &getTuningStrategyOptions() const {
    return _tuningStrategyOptions;
  }

  /**
   * Setter for the tuning strategy option.
   * For possible tuning strategy choices see options::TuningStrategyOption::Value.
   * @param tuningStrategyOptions
   */
  void setTuningStrategyOption(const std::vector<TuningStrategyOption> &tuningStrategyOptions) {
    _tuningStrategyOptions = tuningStrategyOptions;
  }

  /**
   * Getter for the tuning metric option.
   * @return
   */
  [[nodiscard]] const TuningMetricOption &getTuningMetricOption() const { return _autoTunerInfo.tuningMetric; }

  /**
   * Setter for the tuning metric option.
   * For possible tuning metric choices see options::TuningMetricOption::Value.
   * @param tuningMetricOption
   */
  void setTuningMetricOption(TuningMetricOption tuningMetricOption) {
    _autoTunerInfo.tuningMetric = tuningMetricOption;
  }

  /**
   * Getter for the energy sensor
   * @return
   */
  [[nodiscard]] const EnergySensorOption &getEnergySensorOption() const { return _autoTunerInfo.energySensor; }

  /**
   * Setter for the energy sensor
   * @param energySensorOption
   */
  void setEnergySensorOption(EnergySensorOption energySensorOption) {
    _autoTunerInfo.energySensor = energySensorOption;
  }

  /**
   * Setter for the maximal Difference for the bucket distribution.
   * @param MPITuningMaxDifferenceForBucket
   */
  void setMPITuningMaxDifferenceForBucket(double MPITuningMaxDifferenceForBucket) {
    _tuningStrategyFactoryInfo.mpiTuningMaxDifferenceForBucket = MPITuningMaxDifferenceForBucket;
  }

  /**
   * Setter for the maxDensity-Weight in calculation for bucket distribution.
   * @param MPITuningWeightForMaxDensity
   */
  void setMPITuningWeightForMaxDensity(double MPITuningWeightForMaxDensity) {
    _tuningStrategyFactoryInfo.mpiTuningWeightForMaxDensity = MPITuningWeightForMaxDensity;
  }

// Only define the interface for the MPI communicator if AUTOPAS_INCLUDE_MPI=ON
// The internal implementation will use _autopasMPICommunicator with WrapMPI regardless of AUTOPAS_INCLUDE_MPI
#if defined(AUTOPAS_INCLUDE_MPI)
  /**
   * Setter for the MPI communicator that AutoPas uses for potential MPI calls.
   * If not set, MPI_COMM_WORLD will be used.
   * @param comm: communicator (handle)
   */
  void setMPICommunicator(MPI_Comm comm) { _tuningStrategyFactoryInfo.autopasMpiCommunicator = comm; }

  /**
   * Getter for the AutoPas MPI communicator
   * @return communicator
   */
  MPI_Comm getMPICommunicator() { return _tuningStrategyFactoryInfo.autopasMpiCommunicator; }
#endif

  /**
   * Suffix for all output files produced by this instance of AutoPas, e.g. from csv loggers.
   * This is useful when multiple instances of AutoPas exist, especially in an MPI context.
   * @param suffix
   */
  void setOutputSuffix(const std::string &suffix) { _outputSuffix = suffix; }

  /**
   * Getter for the mean rebuild frequency.
   * Helpful for determining the frequency for the dynamic containers as well as for determining fast particles by
   * computing skinPerStep for static container
   * @return Value of the mean rebuild frequency as double
   */
  double getMeanRebuildFrequency() { return _logicHandler->getMeanRebuildFrequency(); }

  /**
   * Set if the tuning information should be logged to a file. It can then be replayed to test other tuning strategies.
   * @param useTuningLogger
   */
  void setUseTuningLogger(bool useTuningLogger) { _useTuningStrategyLoggerProxy = useTuningLogger; }

  /**
   * Set rule file name for the RuleBasedTuning.
   * @param ruleFileName The name of the rule file to use during rule based tuning.
   */
  void setRuleFileName(const std::string &ruleFileName) { _tuningStrategyFactoryInfo.ruleFileName = ruleFileName; }

  /**
   * Set fuzzy rule file name for the RuleBasedTuning.
   * @param fuzzyRuleFileName The name of the fuzzy rule file to use during rule based tuning.
   */
  void setFuzzyRuleFileName(const std::string &fuzzyRuleFileName) {
    _tuningStrategyFactoryInfo.fuzzyRuleFileName = fuzzyRuleFileName;
  }

  /**
   * Get the name / path of the rule file for the RuleBasedTuning.
   * @return
   */
  const std::string &getRuleFileName() const { return _tuningStrategyFactoryInfo.ruleFileName; }

  /**
   * Set the name / path of the model file for the PythonBasedDecisionTreeTuning.
   * @param modelFileName The name of the model file to use during decision tree tuning.
   */
  void setModelFileName(const std::string &modelFileName) { _tuningStrategyFactoryInfo.modelFileName = modelFileName; }

  /**
   * Get the name / path of the model file for the PythonBasedDecisionTreeTuning.
   * @return
   */
  const std::string &getModelFileName() const { return _tuningStrategyFactoryInfo.modelFileName; }

  /**
   * Set the name / path of the pairwise model file for the TreeliteBasedDecisionTreeTuning.
   * @param modelPairwiseFileName The name of the model file to use during decision tree tuning.
   */
  void setModelPairwiseFileName(const std::string &modelPairwiseFileName) {
    _tuningStrategyFactoryInfo.modelPairwiseFileName = modelPairwiseFileName;
  }

  /**
   * Get the name / path of the pairwise model file for the TreeliteBasedDecisionTreeTuning.
   * @return
   */
  const std::string &getModelPairwiseFileName() const { return _tuningStrategyFactoryInfo.modelPairwiseFileName; }

  /**
   * Set the name / path of the triwise model file for the TreeliteBasedDecisionTreeTuning.
   * @param modelTriwiseFileName The name of the model file to use during decision tree tuning.
   */
  void setModelTriwiseFileName(const std::string &modelTriwiseFileName) {
    _tuningStrategyFactoryInfo.modelTriwiseFileName = modelTriwiseFileName;
  }

  /**
   * Get the name / path of the triwise model file for the TreeliteBasedDecisionTreeTuning.
   * @return
   */
  const std::string &getModelTriwiseFileName() const { return _tuningStrategyFactoryInfo.modelTriwiseFileName; }

  /**
   * Set the confidence threshold for the PythonBasedDecisionTreeTuning and TreeliteBasedDecisionTreeTuning.
   * @param confidenceThreshold The confidence threshold to use during decision tree tuning.
   */
  void setConfidenceThreshold(double confidenceThreshold) {
    _tuningStrategyFactoryInfo.confidenceThreshold = confidenceThreshold;
  }
  /**
   * Get the confidence threshold for the PythonBasedDecisionTreeTuning and TreeliteBasedDecisionTreeTuning.
   * @return
   */
  double getConfidenceThreshold() const { return _tuningStrategyFactoryInfo.confidenceThreshold; }
  /**
   * Set the sorting-threshold for traversals that use the CellFunctor
   * If the sum of the number of particles in two cells is greater or equal to that value, the CellFunctor creates a
   * sorted view of the particles to avoid unnecessary distance checks.
   * @param sortingThreshold Sum of the number of particles in two cells from which sorting should be enabled.
   */
  void setSortingThreshold(size_t sortingThreshold) { _sortingThreshold = sortingThreshold; }

  /**
   * Get the sorting-threshold for traversals that use the CellFunctor.
   * @return sorting-threshold
   */
  size_t getSortingThreshold() const { return _sortingThreshold; }

 private:
  autopas::ParticleContainerInterface<Particle_T> &getContainer();

  const autopas::ParticleContainerInterface<Particle_T> &getContainer() const;
  /**
   * Information needed for TuningStrategyFactory::generateTuningStrategy().
   */
  TuningStrategyFactoryInfo _tuningStrategyFactoryInfo{};
  /**
   * Information needed for the AutoTuner.
   */
  AutoTunerInfo _autoTunerInfo{};
  /**
   * Information needed for the LogicHandler.
   */
  LogicHandlerInfo _logicHandlerInfo{};
  /**
   * Whether to insert the tuning strategy logger into the list of strategies.
   */
  bool _useTuningStrategyLoggerProxy{false};
  /**
   * Specifies after how many pair-wise traversals the neighbor lists are to be rebuild, if a rebuild is not triggered
   * earlier by the dynamic rebuild mechanic.
   */
  unsigned int _verletRebuildFrequency{100};
  /**
   * Strategy option for the auto tuner.
   * For possible tuning strategy choices see options::TuningStrategyOption::Value.
   */
  std::vector<TuningStrategyOption> _tuningStrategyOptions{};
  /**
   * List of container types AutoPas can choose from.
   * For possible container choices see options::ContainerOption::Value.
   */
  std::set<ContainerOption> _allowedContainers{ContainerOption::getMostOptions()};
  /**
   * List of pairwise traversals AutoPas can choose from.
   * For possible traversal choices see options::TraversalOption::Value.
   */
  std::unordered_map<InteractionTypeOption::Value, std::set<TraversalOption>> _allowedTraversals{
      {InteractionTypeOption::pairwise, TraversalOption::getMostPairwiseOptions()},
      {InteractionTypeOption::triwise, TraversalOption::getMostTriwiseOptions()}};
  /**
   * List of data layouts AutoPas can choose from for pairwise interactions.
   * For possible data layout choices see options::DataLayoutOption::Value.
   */
  std::unordered_map<InteractionTypeOption::Value, std::set<DataLayoutOption>> _allowedDataLayouts{
      {InteractionTypeOption::pairwise, DataLayoutOption::getMostOptions()},
      {InteractionTypeOption::triwise, DataLayoutOption::getMostOptions()}};
  /**
   * Whether AutoPas is allowed to exploit Newton's third law of motion for pairwise traversals.
   */
  std::unordered_map<InteractionTypeOption::Value, std::set<Newton3Option>> _allowedNewton3Options{
      {InteractionTypeOption::pairwise, Newton3Option::getMostOptions()},
      {InteractionTypeOption::triwise, Newton3Option::getMostOptions()}};
  /**
   * What kind of interactions AutoPas should expect.
   * By default AutoPas is configured to only use pairwise interactions.
   */
  std::set<InteractionTypeOption> _allowedInteractionTypeOptions{InteractionTypeOption::pairwise};
  /**
   * Cell size factor to be used in this container (only relevant for LinkedCells, VerletLists and VerletListsCells).
   */
  std::unique_ptr<NumberSet<double>> _allowedCellSizeFactors{
      std::make_unique<NumberSetFinite<double>>(std::set<double>({1.}))};
  /**
   * Load estimation algorithm to be used for efficient parallelization (only relevant for LCSlicedBalancedTraversal and
   * VLCSlicedBalancedTraversal).
   */
  std::set<LoadEstimatorOption> _allowedLoadEstimators{LoadEstimatorOption::getAllOptions()};
  /**
   * LogicHandler of autopas.
   */
  std::unique_ptr<autopas::LogicHandler<Particle_T>> _logicHandler;

  /**
   * All AutoTuners used in this instance of AutoPas.
   * There can be up to one per interaction type.
   */
  std::unordered_map<InteractionTypeOption::Value, std::unique_ptr<autopas::AutoTuner>> _autoTuners;

  /**
   * Stores whether the mpi communicator was provided externally or not
   */
  bool _externalMPICommunicator{false};
  /**
   * Suffix for all output files produced by this instance of AutoPas, e.g. from csv loggers.
   * This is useful when multiple instances of AutoPas exist, especially in an MPI context.
   */
  std::string _outputSuffix{""};
  /**
   * Number of particles in two cells from which sorting should be performed for traversal that use the CellFunctor
   */
  size_t _sortingThreshold{8};
  /**
   * Helper function to reduce code duplication for all forms of addParticle while minimizing overhead through loops.
   * Triggers reserve() and provides a parallel loop with deliberate scheduling.
   * @tparam F Function type of loopBody: (int) -> void.
   * @param numParticlesToAdd For how many new owned particles should space be allocated.
   * @param numHalosToAdd For how many new halo particles should space be allocated.
   * @param collectionSize Size of the collection from which particles are added.
   * @param loopBody Function to be called in the parallel loop over collectionSize.
   * Typically `[&](auto i) {addParticle(collection[i]);}`.
   */
  template <class F>
  void addParticlesAux(size_t numParticlesToAdd, size_t numHalosToAdd, size_t collectionSize, F loopBody);
};  // class AutoPas
}  // namespace autopas
