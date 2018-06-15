/**
 * @file ParticleContainerInterface.h
 *
 * @date 17 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>

#include "iterators/ParticleIteratorWrapper.h"

namespace autopas {

// consider multiple inheritance or delegation to avoid virtual call to Functor
// FIXME: WRONG DOCUMENTATION
/**
 * The ParticleContainer class stores particles in some object and provides
 * methods to iterate over its particles.
 * @tparam Particle Class for particles
 * @tparam ParticleCell Class for the particle cells
 */
template <class Particle>
class ParticleContainerInterface {
 public:
  ParticleContainerInterface() {}

  /**
   * destructor of ParticleContainerInterface
   */
  virtual ~ParticleContainerInterface() = default;

  /**
   * delete the copy constructor to prevent unwanted copies.
   * No particle container should ever be copied.
   * @param obj
   */
  ParticleContainerInterface(const ParticleContainerInterface &obj) = delete;

  /**
   * delete the copy assignment operator to prevent unwanted copies
   * No particle container should ever be copied.
   * @param other
   * @return
   */
  ParticleContainerInterface &operator=(const ParticleContainerInterface &other) = delete;

  /**
   * adds a particle to the container
   * @param p the particle to be added
   */
  virtual void addParticle(Particle &p) = 0;

  /**
   * adds a particle to the container that lies in the halo region of the
   * container
   * @param haloParticle particle to be added
   */
  virtual void addHaloParticle(Particle &haloParticle) = 0;

  /**
   * deletes all halo particles
   */
  virtual void deleteHaloParticles() = 0;

  /**
   * iterate over all particles using
   * for(auto iter = container.begin(); iter.isValid(); ++iter)
   * @param behavior behavior of the iterator, see IteratorBehavior
   * @return iterator to the first particle
   * @todo implement IteratorBehavior
   */
  virtual ParticleIteratorWrapper<Particle> begin(IteratorBehavior behavior = IteratorBehavior::haloAndOwned) = 0;

  /**
   * iterate over all particles in a specified region
   * for(auto iter = container.getRegionIterator(lowCorner,
   * highCorner);iter.isValid();++iter)
   * @param lowerCorner lower corner of the region
   * @param higherCorner higher corner of the region
   * @param behavior the behavior of the iterator (shall it iterate over halo particles as well?
   * @return iterator to iterate over all particles in a specific region
   */
  virtual ParticleIteratorWrapper<Particle> getRegionIterator(
      std::array<double, 3> lowerCorner, std::array<double, 3> higherCorner,
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) = 0;

  /**
   * get the upper corner of the container
   * @return upper corner of the container
   */
  virtual const std::array<double, 3> &getBoxMax() const = 0;

  /**
   * set the upper corner of the container
   * @param boxMax upper corner to be set
   */
  virtual void setBoxMax(const std::array<double, 3> &boxMax) = 0;

  /**
   * get the lower corner of the container
   * @return lower corner of the container
   */
  virtual const std::array<double, 3> &getBoxMin() const = 0;

  /**
   * set the lower corner of the container
   * @param boxMin lower corner to be set
   */
  virtual void setBoxMin(const std::array<double, 3> &boxMin) = 0;

  /**
   * return the cutoff of the container
   * @return
   */
  virtual double getCutoff() const = 0;

  /**
   * set the cutoff of the container
   * @param cutoff
   */
  virtual void setCutoff(double cutoff) = 0;

  /**
   * updates the container.
   * this resorts particles into appropriate cells, if necessary
   */
  virtual void updateContainer() = 0;

  /**
   * check whether a container is valid, i.e. whether it is safe to use
   * pair-wise interactions or the RegionParticleIteraor right now.
   * @return true if an update is needed, false otherwise
   */
  virtual bool isContainerUpdateNeeded() = 0;
};

}  // namespace autopas
