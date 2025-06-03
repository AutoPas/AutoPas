/**
 * @file ParticleCell.h
 * @date 17.01.2018
 * @author tchipevn
 */

#pragma once

#include "autopas/options/IteratorBehavior.h"
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
 * @tparam Particle_T the type of particles to be stored in the cells
 */
template <class Particle_T>
class ParticleCell {
 public:
  /**
   * The particle type for this cell. Used to refer to the Particle_T typename for an instantiation of ParticleCell.
   */
  using ParticleType = Particle_T;

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
  virtual void addParticle(const Particle_T &p) = 0;

  /**
   * Get the number of all particles stored in this cell (owned, halo and dummy).
   * @return number of particles stored in this cell (owned, halo and dummy).
   */
  [[nodiscard]] virtual size_t size() const = 0;

  /**
   * Get the number of particles with respect to the specified IteratorBehavior.
   * @warning: Since this function counts the number of the respective particles in the internal particle storage, this
   * is in O(n) + lock is required. Only use it when it is absolutely necessary to have the exact number of different
   * particle types like owned or halo. If it is enough to have the whole number of particles (owned + halo + dummy),
   * the function size() can be used.
   * @param behavior Behavior of the iterator, see IteratorBehavior.
   * @return The number of particles with respect to the specified IteratorBehavior.
   */
  [[nodiscard]] virtual size_t getNumberOfParticles(IteratorBehavior behavior = IteratorBehavior::owned) const = 0;

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
  OwnershipState getPossibleParticleOwnerships() const { return _ownershipState; }

  /**
   * Set the type of particles contained in this cell. Possible values:
   * dummy: this cell is empty
   * owned: this cell can ONLY contain owned particles
   * halo: this cell can ONLY contain halo particles
   * ownedOrHalo: this cell can contain owned or halo particles
   * @note: At the moment an ownership of a cell can only be set once.
   * @param state type of particles inside this cell
   */
  void setPossibleParticleOwnerships(OwnershipState state) {
    std::lock_guard<AutoPasLock> guard(this->_cellLock);
    if (_ownershipStateDefined) {
      autopas::utils::ExceptionHandler::exception(
          "ParticleCell::setPossibleParticleOwnerships() can not set OwnershipState of a cell after "
          "it has been set");
    }
    _ownershipState = state;
    _ownershipStateDefined = true;
  }

  /**
   * Get a reference to the lock object for exclusive access to this cell
   * @return reference to the cell lock
   */
  AutoPasLock &getCellLock() const { return _cellLock; }

 protected:
  /**
   * Lock object for exclusive access to this cell. This lock has to be mutable since it is set in
   * getNumberOfParticles() const.
   */
  mutable AutoPasLock _cellLock{};

  /**
   * The particles which can be contained in this cell are determined by the OwnershipState
   */
  OwnershipState _ownershipState{autopas::OwnershipState::owned | autopas::OwnershipState::halo};
  /**
   * Flag that is set to true once OwnershipState has been set to avoid resetting the OwnershipState
   */
  bool _ownershipStateDefined{false};
};

}  // namespace autopas
