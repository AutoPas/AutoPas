/**
 * @file LinkedCellTraversalInterface.h
 * @author seckler
 * @date 09.01.19
 */

#pragma once

#include <vector>

namespace autopas {

/**
 * Interface for traversals used by the LinkedCell class.
 * @tparam ParticleCell
 */
template <class ParticleCell>
class LinkedCellTraversalInterface {
 public:
  /**
   * Traverse all pairs of cells.
   * This function needs to be implemented by derived classes and handles to
   * order in which the cells are traversed.
   * @param cells Vector of cells to traverse
   */
  virtual void traverseCellPairs(std::vector<ParticleCell> &cells) = 0;
};

}  // namespace autopas