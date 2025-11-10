/**
 * @file ParticleContainerInterface.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <vector>

#include "autopas/cells/ParticleCell.h"
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/containers/TraversalInterface.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/tuning/selectors/TraversalSelectorInfo.h"
#include "autopas/utils/AutoPasMacros.h"
#include "autopas/utils/inBox.h"

namespace autopas {

// forward declaration
template <class Particle_T, bool modifiable, bool regionIter>
class ContainerIterator;

/**
 * The ParticleContainerInterface class provides a basic interface for all Containers within AutoPas.
 * It defines method interfaces for addition and deletion of particles, accessing general container
 * properties and creating iterators.
 *
 * @tparam Particle_T Class for particle.
 */
template <class Particle_T>
class ParticleContainerInterface {
 public:
  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle_T;

  /**
   * Get the ParticleCell type as an Enum
   * @return The Cell type as an Enum
   */
  virtual CellType getParticleCellTypeEnum() const = 0;

  /**
   * Default constructor
   */
  ParticleContainerInterface() = default;

  /**
   * Destructor of ParticleContainerInterface.
   */
  virtual ~ParticleContainerInterface() = default;

  /**
   * Delete the copy constructor to prevent unwanted copies.
   * No particle container should ever be copied.
   * @param obj
   */
  ParticleContainerInterface(const ParticleContainerInterface &obj) = delete;

  /**
   * Delete the copy assignment operator to prevent unwanted copies.
   * No particle container should ever be copied.
   * @param other
   * @return
   */
  ParticleContainerInterface &operator=(const ParticleContainerInterface &other) = delete;

  /**
   * Get the ContainerType.
   * @return ContainerOption of the type of this container.
   */
  [[nodiscard]] virtual ContainerOption getContainerType() const = 0;

  /**
   * @copydoc AutoPas::reserve()
   * @param numParticlesHaloEstimate Estimate for the number of halo particles.
   * Reserves space in the container data structure.
   */
  virtual void reserve(size_t numParticles, size_t numParticlesHaloEstimate) = 0;

  /**
   * Adds a particle to the container.
   * @tparam checkInBox Specifies whether a boundary check should be performed. Only disable this if the check has
   * already been performed.
   * @param p The particle to be added.
   */
  template <bool checkInBox = true>
  void addParticle(const Particle_T &p) {
    if constexpr (checkInBox) {
      if (utils::notInBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
        utils::ExceptionHandler::exception(
            "ParticleContainerInterface: Trying to add a particle that is not in the bounding box.\n"
            "Box Min {}\n"
            "Box Max {}\n"
            "{}",
            this->getBoxMin(), this->getBoxMax(), p.toString());
      }
    }
    addParticleImpl(p);
  };

 protected:
  /**
   * Adds a particle to the container.
   * This is an unsafe version of addParticle() and does not perform a boundary check.
   * @param p The particle to be added. This particle is already checked to be inside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be inside of the bounding box!
   */
  virtual void addParticleImpl(const Particle_T &p) = 0;

 public:
  /**
   * Adds a particle to the container that lies in the halo region of the container.
   * @param haloParticle Particle to be added.
   * @tparam checkInBox Specifies whether a boundary check should be performed. Only disable this if the check has
   * already been performed.
   */
  template <bool checkInBox = true>
  void addHaloParticle(const Particle_T &haloParticle) {
    if constexpr (checkInBox) {
      /// @todo do we want a check of the particle not being too far away in here as well?
      if (utils::inBox(haloParticle.getR(), this->getBoxMin(), this->getBoxMax())) {
        utils::ExceptionHandler::exception(
            "ParticleContainerInterface: Trying to add a halo particle that is not outside of in the bounding box.\n"
            "Box Min {}\n"
            "Box Max {}\n"
            "{}",
            this->getBoxMin(), this->getBoxMax(), haloParticle.toString());
      }
    }
    addHaloParticleImpl(haloParticle);
  }

 protected:
  /**
   * Adds a particle to the container that lies in the halo region of the container.
   * This is an unsafe version of addParticle() and does not perform a boundary check.
   * @param haloParticle Particle to be added. This particle is already checked to be outside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be outside of the bounding box!
   */
  virtual void addHaloParticleImpl(const Particle_T &haloParticle) = 0;

 public:
  /**
   * Update a halo particle of the container with the given haloParticle.
   * @param haloParticle Particle to be updated.
   * @return Returns true if the particle was updated, false if no particle could be found.
   */
  virtual bool updateHaloParticle(const Particle_T &haloParticle) = 0;

  /**
   * Rebuilds the neighbor lists for the next traversals.
   * @param traversal The used traversal.
   */
  virtual void rebuildNeighborLists(TraversalInterface *traversal) = 0;

  /**
   * Deletes all halo particles.
   */
  virtual void deleteHaloParticles() = 0;

  /**
   * Deletes all particles.
   */
  virtual void deleteAllParticles() = 0;

  /**
   * Get the number of particles with respect to the specified IteratorBehavior.
   * @warning: Since this function counts the number of the respective particles in the internal particle storage, this
   * is in O(n) + lock is required. Only use it when it is absolutely necessary to have the exact number of different
   * particle types like owned or halo. If it is enough to have the whole number of particles (owned + halo + dummy),
   * the function size() can be used.
   * @param behavior Behavior of the iterator, see IteratorBehavior.
   * @return The number of particles with respect to the specified IteratorBehavior.
   */
  [[nodiscard]] virtual size_t getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::owned) const = 0;

  /**
   * Get the total number of particles saved in the container (owned + halo + dummy).
   * @return Number of particles saved in the container (owned + halo + dummy).
   */
  [[nodiscard]] virtual size_t size() const = 0;

  /**
   * Iterates over all particle multiples (e.g. pairs, triplets) in the container using the given traversal.
   * @param traversal The traversal to use for the iteration.
   */
  virtual void computeInteractions(TraversalInterface *traversal) = 0;

  /**
   * Get the upper corner of the container without halo.
   * @return Upper corner of the container.
   */
  [[nodiscard]] virtual const std::array<double, 3> &getBoxMax() const = 0;

  /**
   * Get the lower corner of the container without halo.
   * @return Lower corner of the container.
   */
  [[nodiscard]] virtual const std::array<double, 3> &getBoxMin() const = 0;

  /**
   * Return the cutoff of the container.
   * @return Cutoff radius.
   */
  [[nodiscard]] virtual double getCutoff() const = 0;

  /**
   * Set the cutoff of the container.
   * @param cutoff
   */
  virtual void setCutoff(double cutoff) = 0;

  /**
   * Return the verletSkin of the container verletSkin
   * @return verletSkin
   */
  [[nodiscard]] virtual double getVerletSkin() const = 0;

  /**
   * Return the number of time-steps since last neighbor list rebuild
   * @note: The value has to be set by setStepsSinceLastRebuild() from outside the container. Otherwise this will always
   * return 0
   * @return steps since last rebuild
   */
  [[nodiscard]] virtual size_t getStepsSinceLastRebuild() const { return _stepsSinceLastRebuild; }

  /**
   * Set the number of time-steps since last neighbor list rebuild
   * @param stepsSinceLastRebuild steps since last neighbor list rebuild
   */
  virtual void setStepsSinceLastRebuild(size_t stepsSinceLastRebuild) {
    _stepsSinceLastRebuild = stepsSinceLastRebuild;
  }

  /**
   * Return the interaction length (cutoff+skin) of the container.
   * @return interaction length
   */
  [[nodiscard]] virtual double getInteractionLength() const = 0;

  /**
   * Updates the container.
   * This deletes halo particles, resorts particles into appropriate cells and might remove particles from the
   * container, if necessary.
   * @param keepNeighborListsValid Defines whether the neighbor lists have to be kept valid.
   * @return A vector of invalid particles that do not belong into the container.
   */
  [[nodiscard]] virtual std::vector<ParticleType> updateContainer(bool keepNeighborListsValid) = 0;

  /**
   * Generates a traversal selector info for this container.
   * @return Traversal selector info for this container.
   */
  [[nodiscard]] virtual TraversalSelectorInfo getTraversalSelectorInfo() const = 0;

  /**
   * Generates a list of all traversals that are theoretically applicable to this container.
   *
   * Traversals might still be not applicable for other reasons so call traversal.isApplicable to be safe!
   *
   * @param interactionType interaction type for which to get all traversals
   * @return Vector of traversal options.
   */
  [[nodiscard]] std::set<TraversalOption> getAllTraversals(const InteractionTypeOption interactionType) const {
    return compatibleTraversals::allCompatibleTraversals(this->getContainerType(), interactionType);
  }

  /**
   * Deletes the given particle as long as this does not compromise the validity of the container.
   * If this is not possible the particle is just marked as deleted.
   * @note This function might be implemented via swap-delete and is thus not completely thread safe.
   * @param particle Reference to the particle that is to be deleted.
   * @return True if the given pointer still points to a new particle.
   */
  virtual bool deleteParticle(Particle_T &particle) = 0;

 protected:
  /**
   * Stores the number of time-steps since last neighbor list rebuild
   * @note: The value has to be set by setStepsSinceLastRebuild() from outside the container. Otherwise this will always
   * be 0
   */
  size_t _stepsSinceLastRebuild{0};
};

}  // namespace autopas
