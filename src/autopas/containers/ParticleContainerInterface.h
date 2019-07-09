/**
 * @file ParticleContainerInterface.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <vector>
#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/iterators/ParticleIteratorWrapper.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/selectors/TraversalSelectorInfo.h"

namespace autopas {

/**
 * The ParticleContainerInterface class provides a basic interface for all Containers within AutoPas.
 * It defines method interfaces for addition and deletion of particles, accessing general container
 * properties and creating iterators.
 *
 * @tparam Particle Class for particles
 */
template <class Particle, class ParticleCell>
class ParticleContainerInterface {
 public:
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
   * Return a enum representing the name of the container class.
   * @return Enum representing the container.
   */
  virtual ContainerOption getContainerType() = 0;

  /**
   * Adds a particle to the container.
   * @param p The particle to be added.
   */
  virtual void addParticle(Particle &p) = 0;

  /**
   * Adds a particle to the container that lies in the halo region of the container.
   * @param haloParticle Particle to be added.
   */
  virtual void addHaloParticle(Particle &haloParticle) = 0;

  /**
   * Deletes all halo particles.
   */
  virtual void deleteHaloParticles() = 0;

  /**
   * Deletes all particles.
   */
  virtual void deleteAllParticles() = 0;

  /**
   * Get the number of particles saved in the container.
   * @return Number of particles in the container.
   */
  virtual unsigned long getNumParticles() = 0;

  /**
   * Iterate over all particles using
   * for(auto iter = container.begin(); iter.isValid(); ++iter) .
   * @param behavior Behavior of the iterator, see IteratorBehavior.
   * @return Iterator to the first particle.
   * @todo implement IteratorBehavior.
   */
  virtual ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) = 0;

  /**
   * Iterate over all particles in a specified region
   * for(auto iter = container.getRegionIterator(lowCorner, highCorner);iter.isValid();++iter) .
   * @param lowerCorner Lower corner of the region
   * @param higherCorner Higher corner of the region
   * @param behavior The behavior of the iterator (shall it iterate over halo particles as well?).
   * @param incSearchRegion Whether to increase the search space (e.g. include more cells)
   * @return Iterator to iterate over all particles in a specific region.
   */
  virtual ParticleIteratorWrapper<Particle> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned, bool incSearchRegion = false) = 0;

  /**
   * Get the upper corner of the container.
   * @return Upper corner of the container.
   */
  virtual const std::array<double, 3> &getBoxMax() const = 0;

  /**
   * Set the upper corner of the container.
   * @param boxMax Upper corner to be set.
   */
  virtual void setBoxMax(const std::array<double, 3> &boxMax) = 0;

  /**
   * Get the lower corner of the container.
   * @return Lower corner of the container.
   */
  virtual const std::array<double, 3> &getBoxMin() const = 0;

  /**
   * Set the lower corner of the container.
   * @param boxMin Lower corner to be set.
   */
  virtual void setBoxMin(const std::array<double, 3> &boxMin) = 0;

  /**
   * Return the cutoff of the container.
   * @return Cutoff radius.
   */
  virtual double getCutoff() const = 0;

  /**
   * Set the cutoff of the container.
   * @param cutoff
   */
  virtual void setCutoff(double cutoff) = 0;

  /**
   * Updates the container.
   * This resorts particles into appropriate cells and moves them to the halo, if necessary.
   */
  virtual void updateContainer() = 0;

  /**
   * Check whether a container is valid, i.e. whether it is safe to use
   * pair-wise interactions or the RegionParticleIteraor right now.
   * @return true if an update is needed, false otherwise
   */
  virtual bool isContainerUpdateNeeded() = 0;

  /**
   * Generates a traversal selector info for this container.
   * @return Traversal selector info for this container.
   */
  virtual std::unique_ptr<TraversalSelectorInfo<ParticleCell>> getTraversalSelectorInfo() = 0;

  /**
   * Generates a list of all traversals that are theoretically applicable to this container.
   *
   * Traversals might still be not applicable for other reasons so call traversal.isApplicable to be safe!
   *
   * @return Vector of traversal options.
   */
  std::set<TraversalOption> getAllTraversals() {
    return compatibleTraversals::allCompatibleTraversals(this->getContainerType());
  }
};

}  // namespace autopas
