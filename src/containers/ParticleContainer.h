/*
 * ParticleContainer.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECONTAINER_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECONTAINER_H_

#include <array>
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
class ParticleContainer {
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

  //virtual void init() {}

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
   * function to iterate over all pairs of particles in an array of structures
   * setting. This function only handles short-range interactions.
   * @param f functor that describes the pair-potential
   */
  virtual void iteratePairwiseAoS(Functor<Particle, ParticleCell> *f) = 0;

  /**
   * function to iterate over all pairs of particles in a structure of array
   * setting. This function is often better vectorizable.
   * @param f
   */
  virtual void iteratePairwiseSoA(Functor<Particle, ParticleCell> *f) = 0;

  /**
   * iterate over all particles using
   * for(auto iter = container.begin(); iter.isValid(); ++iter)
   * @return iterator to the first particle
   */
  ParticleIterator<Particle, ParticleCell> begin() { return ParticleIterator<Particle, ParticleCell>(&_data); }

  /**
   * iterate over all particles in a specified region
   * for(auto iter = container.getRegionIterator(lowCorner, highCorner);iter.isValid();++iter)
   * @param lowerCorner lower corner of the region
   * @param higherCorner higher corner of the region
   * @return iterator to iterate over all particles in a specific region
   */
  RegionParticleIterator<Particle, ParticleCell> getRegionIterator(std::array<double, 3> lowerCorner,
                                   std::array<double, 3> higherCorner) {
    return RegionParticleIterator<Particle, ParticleCell>(&_data, lowerCorner, higherCorner);
  }

  /**
   * get the upper corner of the container
   * @return upper corner of the container
   */
  const std::array<double, 3> &getBoxMax() const { return _boxMax; }

  /**
   * set the upper corner of the container
   * @param boxMax upper corner to be set
   */
  void setBoxMax(const std::array<double, 3> &boxMax) { _boxMax = boxMax; }

  /**
   * get the lower corner of the container
   * @return lower corner of the container
   */
  const std::array<double, 3> &getBoxMin() const { return _boxMin; }

  /**
   * set the lower corner of the container
   * @param boxMin lower corner to be set
   */
  void setBoxMin(const std::array<double, 3> &boxMin) { _boxMin = boxMin; }

  /**
   * return the cutoff of the container
   * @return
   */
  double getCutoff() const { return _cutoff; }

  /**
   * set the cutoff of the container
   * @param cutoff
   */
  void setCutoff(double cutoff) { _cutoff = cutoff; }

  /**
   * updates the container.
   * this resorts particles into appropriate cells, if necessary
   */
  virtual void updateContainer() = 0;

 protected:
  /**
   * vector of particle cells.
   * All particle containers store their particles in ParticleCells. This is the common vector for this purpose.
   */
  std::vector<ParticleCell> _data;

 private:
  std::array<double, 3> _boxMin;
  std::array<double, 3> _boxMax;
  double _cutoff;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECONTAINER_H_ */
