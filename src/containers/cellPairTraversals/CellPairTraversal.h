/**
 * @file CellPairTraversal.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <vector>

namespace autopas {

/**
 * A cell pair traversal.
 * This class handles traversals through the cell structures.
 * Derived classes handle the order through which the cells are traversed.
 * @tparam ParticleCell type of cells.
 * @tparam CellFunctor a cell functor.
 */
template <class ParticleCell, class CellFunctor>
class CellPairTraversal {
 public:
  /**
   * Constructor of CellPairTraversal.
   * @param dims the dimensions of the cellblock.
   * @param cellFunctor the cell functor.
   */
  CellPairTraversal(const std::array<unsigned long, 3> &dims, CellFunctor *cellFunctor)
      : _cellsPerDimension(dims), _cellFunctor(cellFunctor) {}

  /**
   * Destructor of CellPairTraversal.
   */
  virtual ~CellPairTraversal() = default;

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  virtual bool isApplicable() = 0;

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

  /**
   * The cell functor which defines the interaction of the particles between two
   * specific cells.
   */
  CellFunctor *_cellFunctor;
};

}  // namespace autopas
