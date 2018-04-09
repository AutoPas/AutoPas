/*
 * ParticleCell.h
 *
 *  Created on: 17 Jan 2018
 *      Author: tchipevn
 */

#ifndef DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECELL_H_
#define DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECELL_H_

namespace autopas {

/**
 * Class for Cells of Particles.
 * The class handles to storage particles and provides an interface to add the
 * particles
 * @tparam Particle the type of particles to be stored in the cells
 */
template <class Particle, class Iterator, class Derived>
class ParticleCell {
 public:
  /**
   * destructor of ParticleCell
   */
  virtual ~ParticleCell() = default;

  /**
   * adds a Particle to the cell
   * @param p the particle to be added
   */
  virtual void addParticle(Particle &p) = 0;

  /**
   * Get an iterator to the start of a ParticleCell
   * normal use:
   * for(auto iter = cell.begin(); iter.isValid; ++iter){...}
   * @return the iterator
   */
  virtual Iterator begin() {
    return Iterator(reinterpret_cast<Derived *>(this));
  }

  /**
   * Get the number of particles stored in this cell
   * @return number of particles stored in this cell
   */
  virtual unsigned long numParticles() const = 0;

  /**
   * Check if the cell is not empty.
   * @return true if at least one particle is stored in this cell
   */
  virtual bool isNotEmpty() const = 0;

  /**
   * Deletes all particles in this cell
   */
  virtual void clear() = 0;

  /**
   * Deletes the index-th particle
   * @param index the index of the particle that shall be deleted
   */
  virtual void deleteByIndex(int index) = 0;
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECELL_H_ */
