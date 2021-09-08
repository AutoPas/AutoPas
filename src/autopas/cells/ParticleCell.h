/**
 * @file ParticleCell.h
 * @date 17.01.2018
 * @author tchipevn
 */

#pragma once

#include <autopas/utils/inBox.h>

#include "autopas/iterators/SingleCellIteratorWrapper.h"

namespace autopas {

/**
 * The ParticleCell Type as an Enum.
 */
enum class CellType {
  /**
   * FullParticleCell : Default cell type for almost everything.
   */
  FullParticleCell,
  /**
   * ReferenceParticleCell : Cell holding only references instead of actual particle objects.
   */
  ReferenceParticleCell,
  /**
   * ClusterTower : Tower for the 2D tower structure of VerletClusterLists.
   */
  ClusterTower,
  /**
   * SortedCellView : Holds pointers to particles sorted by their position projected along a vector.
   */
  SortedCellView,
  /**
   * Currently unused.
   */
  IsNoCell
};

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
  [[nodiscard]] constexpr bool end() const { return false; }

  /**
   * Get the number of particles stored in this cell.
   * @return number of particles stored in this cell
   */
  [[nodiscard]] virtual unsigned long numParticles() const = 0;

  /**
   * Check if the cell is not empty.
   * @return true if at least one particle is stored in this cell
   */
  [[nodiscard]] virtual bool isNotEmpty() const = 0;

  /**
   * Deletes all particles in this cell.
   */
  virtual void clear() = 0;

  /**
   * Deletes all dummy particles in this cell.
   */
  virtual void deleteDummyParticles() = 0;

  /**
   * Get the ParticleCell type as an ParticleCellTypeEnum
   * @return The Cell type as an Enum
   */
  virtual CellType getParticleCellTypeAsEnum() = 0;

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
  [[nodiscard]] virtual std::array<double, 3> getCellLength() const = 0;

  /**
   * Check whether particle is positioned within box
   * @param p particle to check
   * @param lowerCorner lower corner of box
   * @param higherCorner higher corner of box
   * @returns true if the given particle is positioned within the box
   */
  inline bool isParticleInRegion(Particle &p, const std::array<double, 3> &lowerCorner,
                                 const std::array<double, 3> &higherCorner) {
    return utils::inBox(p.getR(), lowerCorner, higherCorner);
  }
};

}  // namespace autopas
