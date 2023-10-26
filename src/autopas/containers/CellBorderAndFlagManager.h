/**
 * @file CellBorderAndFlagManager.h
 * @author seckler
 * @date 14.05.18
 */

#pragma once
#include <cstdlib>

namespace autopas::internal {
/**
 * Interface class to handle cell borders and cell types of cells.
 * @todo: add cell border handling
 */
class CellBorderAndFlagManager {
 public:
  /**
   * The index type to access the particle cells.
   */
  using index_t = std::size_t;

  /**
   * Cestructor
   */
  virtual ~CellBorderAndFlagManager() = default;

  /**
   * Checks if the cell with the one-dimensional index index1d can contain halo particles.
   * @param index1d the one-dimensional index of the cell that should be checked
   * @return true if the cell can contain halo particles
   */
  [[nodiscard]] virtual bool cellCanContainHaloParticles(index_t index1d) const = 0;

  /**
   * Checks if the cell with the one-dimensional index index1d can contain owned particles.
   * @param index1d the one-dimensional index of the cell that should be checked
   * @return true if the cell can contain owned particles
   */
  [[nodiscard]] virtual bool cellCanContainOwnedParticles(index_t index1d) const = 0;

  /**
   * Checks if cell with index1d can be ignored for iteration with currently selected behavior.
   * @param index1d 1d index of checked cell
   * @param behavior @see IteratorBehavior
   * @return false if this cell can contain particles that would be affected by current behavior
   */
  [[nodiscard]] bool ignoreCellForIteration(index_t index1d, IteratorBehavior behavior) const {
    if ((behavior & IteratorBehavior::halo) and cellCanContainHaloParticles(index1d)) {
      return false;
    }
    if ((behavior & IteratorBehavior::owned) and cellCanContainOwnedParticles(index1d)) {
      return false;
    }
    if (behavior & IteratorBehavior::dummy) {
      return false;
    }
    return true;
  }
};
}  // namespace autopas::internal
