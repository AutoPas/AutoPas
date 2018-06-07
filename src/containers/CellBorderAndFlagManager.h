/**
 * @file CellBorderAndFlagManager.h
 * @author seckler
 * @date 14.05.18
 */

#pragma once

namespace autopas {
/**
 * interface class to handle cell borders and cell types of cells.
 * @todo: add cell border handling
 */
class CellBorderAndFlagManager {
  /**
   * the index type to access the particle cells
   */
  typedef std::size_t index_t;

 public:

  /**
   * destructor
   */
  virtual ~CellBorderAndFlagManager() = default;

  /**
   * checks if the cell with the one-dimensional index index1d is a halocell.
   * @param index1d the one-dimensional index of the cell that should be checked
   * @return true if the cell is a halo cell
   */
  virtual bool isHaloCell(index_t index1d) const = 0;

  /**
   * checks if the cell with the one-dimensional index index1d is owned by this
   * process, i.e. it is NOT a halo cell.
   * @param index1d the one-dimensional index of the cell that should be checked
   * @return true if the cell is owned by the current process.
   */
  virtual bool isOwningCell(index_t index1d) const = 0;
};
}  // namespace autopas
