/**
 * @file ParticleCell.h
 * @date 17.01.2018
 * @author tchipevn
 */

#pragma once

#include "autopas/iterators/SingleCellIteratorWrapper.h"

/**
 * The ParticleCell Type as an Enum.
 */
enum ParticleCellTypeEnum {
  FullParticleCellEnum,
  ReferenceParticleCellEnum,
  ClusterTowerEnum,
  SortedCellViewEnum,
  IsNoCellEnum
};

namespace autopas {

/**
 * Class for Cells of Particles.
 * The class handles to storage particles and provides an interface to add the
 * particles
 * @tparam Particle the type of particles to be stored in the cells
 * @tparam Iterator the type of the iterator iterate through the particles in
 * this cell
 */
template <class Particle>
class ParticleCell {
 public:
  /**
   * The particle type for this cell.
   */
  using ParticleType = Particle;

  /**
   * Destructor of ParticleCell.
   */
  virtual ~ParticleCell() = default;

  /**
   * Adds a Particle to the cell.
   * @param p the particle to be added
   */
  virtual void addParticle(const Particle &p) = 0;

  /**
   * Get an iterator to the start of a ParticleCell.
   * normal use:
   * for(auto iter = cell.begin(); iter.isValid; ++iter){...}
   * @return the iterator
   */
  virtual SingleCellIteratorWrapper<Particle, true> begin() = 0;

  /**
   * @copydoc begin()
   * @note const version
   */
  virtual SingleCellIteratorWrapper<Particle, false> begin() const = 0;

  /**
   * End expression for all cells, this simply returns false.
   * Allows range-based for loops.
   * @return false
   */
  constexpr bool end() const { return false; }

  /**
   * Get the number of particles stored in this cell.
   * @return number of particles stored in this cell
   */
  virtual unsigned long numParticles() const = 0;

  /**
   * Check if the cell is not empty.
   * @return true if at least one particle is stored in this cell
   */
  virtual bool isNotEmpty() const = 0;

  /**
   * Deletes all particles in this cell.
   */
  virtual void clear() = 0;

  /**
   * Deletes all dummy particles in this cell.
   */
  virtual void deleteDummyParticles() = 0;

  /**
   * Get the ParticleCell type as an Enum
   * @return The Cell type as an Enum
   */
  virtual ParticleCellTypeEnum getParticleCellTypeAsEnum() = 0;

  /**
   * Deletes the index-th particle.
   * @param index the index of the particle that shall be deleted
   */
  virtual void deleteByIndex(size_t index) = 0;

  /**
   * Set the side lengths of this cell.
   * @param cellLength cell side length
   */
  virtual void setCellLength(std::array<double, 3> &cellLength) = 0;

  /**
   * Get the side lengths of this cell.
   * @return cell side length
   */
  virtual std::array<double, 3> getCellLength() const = 0;
};

}  // namespace autopas
