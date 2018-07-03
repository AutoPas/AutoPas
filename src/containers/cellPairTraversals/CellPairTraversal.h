/**
 * @file CellPairTraversal.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <vector>
#include <selectors/TraversalSelector.h>
#include "CellPairTraversalInterface.h"

namespace autopas {

/**
 * A cell pair traversal.
 * This class handles traversals through the cell structures.
 * Derived classes handle the order through which the cells are traversed.
 * @tparam ParticleCell type of cells.
 */
template <class ParticleCell>
class CellPairTraversal : public CellPairTraversalInterface {
 public:
  /**
   * Constructor of CellPairTraversal.
   * @param dims the dimensions of the cellblock.
   */
  CellPairTraversal(const std::array<unsigned long, 3> &dims)
      : _cellsPerDimension(dims) {}

  /**
   * Destructor of CellPairTraversal.
   */
  virtual ~CellPairTraversal() = default;


  /**
   * Resets the cell structure of the traversal.
   * @param dims
   */
  virtual void rebuild(const std::array<unsigned long, 3> &dims) { _cellsPerDimension = dims; };

  /**
   * Traverse all pairs of cells.
   * This function needs to be implemented by derived classes and handles to
   * order in which the cells are traversed.
   * @param cells Vector of cells to traverse
   */
  virtual void traverseCellPairs(std::vector<ParticleCell> &cells) = 0;

 protected:
  /**
   * The dimensions of the cellblock.
   * The dimensions are the number of cells in x, y and z direction.
   */
  std::array<unsigned long, 3> _cellsPerDimension;

};

}  // namespace autopas
