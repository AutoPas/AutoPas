/*
 * ParticleContainer.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#pragma once

#include <array>
#include "ParticleContainerInterface.h"
#include "iterators/ParticleIterator.h"
#include "iterators/RegionParticleIterator.h"
#include "pairwiseFunctors/Functor.h"

namespace autopas {

// consider multiple inheritance or delegation vor avoidane of virtual call to
// Functor
/**
 * The ParticleContainer class stores particles in some object and provides
 * methods to iterate over its particles.
 * @tparam Particle Class for particles
 * @tparam ParticleCell Class for the particle cells
 */
template <class Particle, class ParticleCell>
class ParticleContainer : public ParticleContainerInterface<Particle> {
 public:
  /**
   * Constructor of ParticleContainer
   * @param boxMin
   * @param boxMax
   * @param cutoff
   */
  ParticleContainer(const std::array<double, 3> boxMin,
                    const std::array<double, 3> boxMax, double cutoff)
      : _data(), _boxMin(boxMin), _boxMax(boxMax), _cutoff(cutoff) {}

  /**
   * destructor of ParticleContainer
   */
  virtual ~ParticleContainer() = default;

  /**
   * delete the copy constructor to prevent unwanted copies.
   * No particle container should ever be copied.
   * @param obj
   */
  ParticleContainer(const ParticleContainer &obj) = delete;

  /**
   * delete the copy assignment operator to prevent unwanted copies
   * No particle container should ever be copied.
   * @param other
   * @return
   */
  ParticleContainer &operator=(const ParticleContainer &other) = delete;

  // virtual void init() {}

  /**
   * adds a particle to the container
   * @param p the particle to be added
   */
  virtual void addParticle(Particle &p) override = 0;

  /**
   * adds a particle to the container that lies in the halo region of the
   * container
   * @param haloParticle particle to be added
   */
  virtual void addHaloParticle(Particle &haloParticle) override = 0;

  /**
   * deletes all halo particles
   */
  virtual void deleteHaloParticles() override = 0;

  /**
   * function to iterate over all pairs of particles in an array of structures
   * setting. This function only handles short-range interactions.
   * @param f functor that describes the pair-potential
   * @param useNewton3 defines whether newton3 should be used
   */
  virtual void iteratePairwiseAoS(Functor<Particle, ParticleCell> *f,
                                  bool useNewton3 = true) = 0;

  /**
   * function to iterate over all pairs of particles in a structure of array
   * setting. This function is often better vectorizable.
   * @param f functor that describes the pair-potential
   * @param useNewton3 defines whether newton3 should be used
   */
  virtual void iteratePairwiseSoA(Functor<Particle, ParticleCell> *f,
                                  bool useNewton3 = true) = 0;

  /**
   * iterate over all particles using
   * for(auto iter = container.begin(); iter.isValid(); ++iter)
   * @return iterator to the first particle
   * @todo implement IteratorBehavior
   */
  ParticleIteratorWrapper<Particle> begin(
      IteratorBehavior behavior = IteratorBehavior::haloAndOwned) override {
    return ParticleIteratorWrapper<Particle>(
        new ParticleIterator<Particle, ParticleCell>(&_data));
  }

  /**
   * iterate over all particles in a specified region
   * for(auto iter = container.getRegionIterator(lowCorner,
   * highCorner);iter.isValid();++iter)
   * @param lowerCorner lower corner of the region
   * @param higherCorner higher corner of the region
   * @return iterator to iterate over all particles in a specific region
   */
  ParticleIteratorWrapper<Particle> getRegionIterator(
      std::array<double, 3> lowerCorner,
      std::array<double, 3> higherCorner) override {
    return ParticleIteratorWrapper<Particle>(
        new RegionParticleIterator<Particle, ParticleCell>(&_data, lowerCorner,
                                                           higherCorner));
  }

  /**
   * get the upper corner of the container
   * @return upper corner of the container
   */
  const std::array<double, 3> &getBoxMax() const override final {
    return _boxMax;
  }

  /**
   * set the upper corner of the container
   * @param boxMax upper corner to be set
   */
  void setBoxMax(const std::array<double, 3> &boxMax) override final {
    _boxMax = boxMax;
  }

  /**
   * get the lower corner of the container
   * @return lower corner of the container
   */
  const std::array<double, 3> &getBoxMin() const override final {
    return _boxMin;
  }

  /**
   * set the lower corner of the container
   * @param boxMin lower corner to be set
   */
  void setBoxMin(const std::array<double, 3> &boxMin) override final {
    _boxMin = boxMin;
  }

  /**
   * return the cutoff of the container
   * @return
   */
  double getCutoff() const override final { return _cutoff; }

  /**
   * set the cutoff of the container
   * @param cutoff
   */
  void setCutoff(double cutoff) override final { _cutoff = cutoff; }

  /**
   * updates the container.
   * this resorts particles into appropriate cells, if necessary
   */
  virtual void updateContainer() override = 0;

  /**
   * check whether a container is valid, i.e. whether it is safe to use
   * pair-wise interactions or the RegionParticleIteraor right now.
   */
  virtual bool isContainerUpdateNeeded() override = 0;

 protected:
  /**
   * vector of particle cells.
   * All particle containers store their particles in ParticleCells. This is the
   * common vector for this purpose.
   */
  std::vector<ParticleCell> _data;

 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _cutoff;
};

} /* namespace autopas */
