/*
 * CellPairTraversal.h
 *
 *  Created on: 22 Jan 2018
 *      Author: tchipevn
 */

#ifndef SRC_CONTAINERS_CELLPAIRTRAVERSALS_CELLPAIRTRAVERSAL_H_
#define SRC_CONTAINERS_CELLPAIRTRAVERSALS_CELLPAIRTRAVERSAL_H_

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
class CellPairTraversals {
 public:
  /**
   * Constructor of CellPairTraversals.
   * @param cells the vector of cells.
   * @param dims the dimensions of the cellblock.
   * @param cellFunctor the cell functor.
   */
  CellPairTraversals(std::vector<ParticleCell> &cells,
                     const std::array<unsigned long, 3> &dims,
                     CellFunctor *cellFunctor)
      : _cells(&cells), _cellsPerDimension(dims), _cellFunctor(cellFunctor) {}

  /**
   * Destructor of CellPairTraversals.
   */
  virtual ~CellPairTraversals() = default;

  /**
   * Checks if the traversal is applicable to the current state of the domain.
   * @return true iff the traversal can be applied.
   */
  virtual bool isApplicable() = 0;

  /**
   * Resets the cell structure of the traversal.
   * @param cells
   * @param dims
   */
  virtual void rebuild(std::vector<ParticleCell> &cells,
                       const std::array<unsigned long, 3> &dims) {
    _cells = &cells;
    _cellsPerDimension = dims;
  };

  /**
   * Traverse all pairs of cells.
   * This function needs to be implemented by derived classes and handles to
   * order in which the cells are traversed.
   */
  virtual void traverseCellPairs() = 0;

 protected:
  /**
   * The vector of cells.
   */
  std::vector<ParticleCell> *_cells;

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

} /* namespace autopas */

#endif /* SRC_CONTAINERS_CELLPAIRTRAVERSALS_CELLPAIRTRAVERSAL_H_ */
