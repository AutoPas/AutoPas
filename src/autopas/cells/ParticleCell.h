/**
 * @file ParticleCell.h
 * @date 17.01.2018
 * @author tchipevn
 */

#pragma once

#include "autopas/particles/OwnershipState.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/inBox.h"

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
   * Default destructor.
   */
  virtual ~ParticleCell() = default;

  /**
   * Default default constructor.
   */
  explicit ParticleCell() = default;

  /**
   * Default move constructor.
   * @param other
   */
  ParticleCell(ParticleCell &&other) noexcept = default;

  /**
   * Copy constructor that creates a new default constructed lock for the new cell.
   * @param other
   */
  ParticleCell(const ParticleCell &other) : _cellLock(){};

  /**
   * Adds a Particle to the cell.
   * @param p the particle to be added
   */
  virtual void addParticle(const Particle &p) = 0;

  /**
   * Get the number of particles stored in this cell.
   * @return number of particles stored in this cell
   */
  [[nodiscard]] virtual unsigned long numParticles() const = 0;

  /**
   * Check if the cell is empty.
   * @return true if no particles are stored in this cell.
   */
  [[nodiscard]] virtual bool isEmpty() const = 0;

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
   * Get the type of particles contained in this cell. Possible values:
   * dummy: this cell is empty
   * owned: this cell can ONLY contain owned particles
   * halo: this cell can ONLY contain halo particles
   * ownedOrHalo: this cell can contain owned or halo particles
   * @return type of particles inside this cell
   */
  const OwnershipState getPossibleParticleOwnerships() {
    if (_numOwnedParticles > 0 and _numHaloParticles == 0) {
      return OwnershipState::owned;
    } else if (_numHaloParticles > 0 and _numOwnedParticles == 0) {
      return OwnershipState::halo;
    } else {
      return (OwnershipState::owned | OwnershipState::halo);
    }
  }

  /**
   * Lock object for exclusive access to this cell.
   */
  AutoPasLock _cellLock{};

 protected:
  int _numOwnedParticles{0};
  int _numHaloParticles{0};
};

}  // namespace autopas
