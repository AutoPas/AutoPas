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
#include "autopas/iterators/ParticleIteratorWrapper.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/IteratorBehavior.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/selectors/TraversalSelectorInfo.h"
#include "autopas/utils/AutoPasMacros.h"
#include "autopas/utils/inBox.h"

namespace autopas {
/**
 * The ParticleContainerInterface class provides a basic interface for all Containers within AutoPas.
 * It defines method interfaces for addition and deletion of particles, accessing general container
 * properties and creating iterators.
 *
 * @tparam Particle Class for particle.
 */
template <class Particle>
class ParticleContainerInterface {
 public:
  /**
   *  Type of the Particle.
   */
  using ParticleType = Particle;

  /**
   * Get the ParticleCell type as an Enum
   * @return The Cell type as an Enum
   */
  virtual CellType getParticleCellTypeEnum() = 0;

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
   * Adds a particle to the container.
   * @tparam checkInBox Specifies whether a boundary check should be performed. Only disable this if the check has
   * already been performed.
   * @param p The particle to be added.
   */
  template <bool checkInBox = true>
  void addParticle(const Particle &p) {
    if constexpr (not checkInBox) {
      addParticleImpl(p);
    } else {
      if (utils::inBox(p.getR(), this->getBoxMin(), this->getBoxMax())) {
        addParticleImpl(p);
      } else {
        utils::ExceptionHandler::exception(
            "ParticleContainerInterface: Trying to add a particle that is not in the bounding box.\n"
            "Box Min {}\n"
            "Box Max {}\n"
            "{}",
            this->getBoxMin(), this->getBoxMax(), p.toString());
      }
    }
  };

 protected:
  /**
   * Adds a particle to the container.
   * This is an unsafe version of addParticle() and does not perform a boundary check.
   * @param p The particle to be added. This particle is already checked to be inside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be inside of the bounding box!
   */
  virtual void addParticleImpl(const Particle &p) = 0;

 public:
  /**
   * Adds a particle to the container that lies in the halo region of the container.
   * @param haloParticle Particle to be added.
   * @tparam checkInBox Specifies whether a boundary check should be performed. Only disable this if the check has
   * already been performed.
   */
  template <bool checkInBox = true>
  void addHaloParticle(const Particle &haloParticle) {
    if constexpr (not checkInBox) {
      addHaloParticleImpl(haloParticle);
    } else {
      /// @todo do we want a check of the particle not being too far away in here as well?
      if (utils::inBox(haloParticle.getR(), this->getBoxMin(), this->getBoxMax())) {
        std::stringstream error;
        utils::ExceptionHandler::exception(
            "Trying to add a halo particle that is not outside of in the bounding box.\n"
            "Box Min {}\n"
            "Box Max {}\n",
            utils::ArrayUtils::to_string(this->getBoxMin()), utils::ArrayUtils::to_string(this->getBoxMax()),
            haloParticle.toString());
      } else {
        addHaloParticleImpl(haloParticle);
      }
    }
  }

 protected:
  /**
   * Adds a particle to the container that lies in the halo region of the container.
   * This is an unsafe version of addParticle() and does not perform a boundary check.
   * @param haloParticle Particle to be added. This particle is already checked to be outside of the bounding box.
   * @note Only call this function if the position of the particle is guaranteed to be outside of the bounding box!
   */
  virtual void addHaloParticleImpl(const Particle &haloParticle) = 0;

 public:
  /**
   * Update a halo particle of the container with the given haloParticle.
   * @param haloParticle Particle to be updated.
   * @return Returns true if the particle was updated, false if no particle could be found.
   */
  virtual bool updateHaloParticle(const Particle &haloParticle) = 0;

  /**
   * Rebuilds the neighbor lists.
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
   * Get the total number of particles saved in the container (owned and halo but not dummies).
   * @return Number of particles in the container.
   */
  [[nodiscard]] virtual unsigned long getNumberOfParticles() const = 0;

  /**
   * Iterate over all particles using
   * for(auto iter = container.begin(); iter.isValid(); ++iter) .
   * @param behavior Behavior of the iterator, see IteratorBehavior.
   * @note Default argument necessary to enable range based for loops.
   * @return Iterator to the first particle.
   */
  [[nodiscard]] virtual ParticleIteratorWrapper<ParticleType, true> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) = 0;

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   * @note const version
   */
  [[nodiscard]] virtual ParticleIteratorWrapper<ParticleType, false> begin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) const = 0;

  /**
   * @copydoc autopas::ParticleContainerInterface::begin()
   * @note cbegin will guarantee to return a const_iterator.
   */
  [[nodiscard]] virtual ParticleIteratorWrapper<ParticleType, false> cbegin(
      IteratorBehavior behavior = autopas::IteratorBehavior::ownedOrHalo) const final {
    return begin(behavior);
  };

  /**
   * Iterate over all particles in a specified region
   * for(auto iter = container.getRegionIterator(lowCorner, highCorner);iter.isValid();++iter) .
   * @param lowerCorner Lower corner of the region
   * @param higherCorner Higher corner of the region
   * @param behavior The behavior of the iterator (shall it iterate over halo particles as well?).
   * @return Iterator to iterate over all particles in a specific region.
   */
  [[nodiscard]] virtual ParticleIteratorWrapper<ParticleType, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior) = 0;

  /**
   * @copydoc autopas::ParticleContainerInterface::getRegionIterator()
   * @note const version
   */
  [[nodiscard]] virtual ParticleIteratorWrapper<ParticleType, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior) const = 0;

  /**
   * End expression for all containers, this simply returns false.
   * Allows range-based for loops.
   * @return false
   */
  [[nodiscard]] constexpr bool end() const { return false; }

  /**
   * Iterates over all particle pairs in the container.
   * @param traversal The traversal to use for the iteration.
   */
  virtual void iteratePairwise(TraversalInterface *traversal) = 0;

  /**
   * Get the upper corner of the container.
   * @return Upper corner of the container.
   */
  [[nodiscard]] virtual const std::array<double, 3> &getBoxMax() const = 0;

  /**
   * Set the upper corner of the container.
   * @param boxMax Upper corner to be set.
   */
  virtual void setBoxMax(const std::array<double, 3> &boxMax) = 0;

  /**
   * Get the lower corner of the container.
   * @return Lower corner of the container.
   */
  [[nodiscard]] virtual const std::array<double, 3> &getBoxMin() const = 0;

  /**
   * Set the lower corner of the container.
   * @param boxMin Lower corner to be set.
   */
  virtual void setBoxMin(const std::array<double, 3> &boxMin) = 0;

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
   * Return the verletSkin of the container verletSkinPerTimestep*rebuildFrequency
   * @return verletSkin
   */
  [[nodiscard]] virtual double getVerletSkin() const = 0;

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
   * @return Vector of traversal options.
   */
  [[nodiscard]] std::set<TraversalOption> getAllTraversals() const {
    return compatibleTraversals::allCompatibleTraversals(this->getContainerType());
  }
};

}  // namespace autopas
