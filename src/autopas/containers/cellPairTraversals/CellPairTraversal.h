/**
 * @file CellPairTraversal.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <vector>

#include "autopas/containers/TraversalInterface.h"

namespace autopas {

/**
 * A cell pair traversal.
 * This class handles traversals through the cell structures.
 * Derived classes handle the order through which the cells are traversed.
 * @tparam ParticleCell type of cells.
 */
template <class ParticleCell>
class CellPairTraversal : public TraversalInterface {
 public:
  /**
   * Constructor of CellPairTraversal.
   * @param dims the dimensions of the cellblock.
   * @param dataLayout The data layout with which this traversal should be initialised.
   * @param useNewton3 Parameter to specify whether the traversal makes use of newton3 or not.
   */
  explicit CellPairTraversal(const std::array<unsigned long, 3> &dims, DataLayoutOption dataLayout, bool useNewton3)
      : TraversalInterface(dataLayout, useNewton3), _cellsPerDimension(dims), _cells(nullptr) {}

  /**
   * Destructor of CellPairTraversal.
   */
  ~CellPairTraversal() override = default;

  /**
   * Sets the cells to iterate over. Should always be called before initTraversal().
   * @param cells The cells to iterate over.
   */
  virtual void setCellsToTraverse(std::vector<ParticleCell> &cells) { _cells = &cells; }

  /**
   * Set the sorting-threshold for traversals that use the CellFunctor
   * If the sum of the number of particles in two cells is greater or equal to that value, the CellFunctor creates a
   * sorted view of the particles to avoid unnecessary distance checks.
   * @param sortingThreshold Sum of the number of particles in two cells from which sorting should be enabled
   */
  virtual void setSortingThreshold(size_t sortingThreshold) = 0;

 protected:
  /**
   * The dimensions of the cellblock.
   * The dimensions are the number of cells in x, y and z direction.
   */
  std::array<unsigned long, 3> _cellsPerDimension;

  /**
   * The cells to traverse.
   */
  std::vector<ParticleCell> *_cells;
};

}  // namespace autopas