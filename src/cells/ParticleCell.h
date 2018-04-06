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
 * The class handles to storage particles and provides an interface to add the particles
 * @tparam Particle the type of particles to be stored in the cells
 */
template <class Particle>
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
   * Get a particle at a specific index
   * @todo enhance documentation: how does the input parameter look like?
   * @param i index of the particle
   * @param rmm_or_not_pointer returns the particle.
   */
  virtual void moleculesAt(
      int i,
      Particle *&rmm_or_not_pointer) = 0;  // TODO: consider
                                           // making return
                                           // type Particle*

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
};

} /* namespace autopas */

#endif /* DEPENDENCIES_EXTERNAL_AUTOPAS_SRC_PARTICLECELL_H_ */
