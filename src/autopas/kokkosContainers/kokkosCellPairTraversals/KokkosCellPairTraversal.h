/**
 * @file KokkosCellPairTraversal.h
 *
 * @date 13 Dec 2021
 * @author lgaertner
 */

#pragma once

#include <array>
#include <vector>

#include "autopas/cells/KokkosParticleCell.h"
#include "autopas/containers/TraversalInterface.h"

namespace autopas {

/**
 * A cell pair traversal.
 * This class handles traversals through the cell structures.
 * Derived classes handle the order through which the cells are traversed.
 * @tparam ParticleCell type of cells.
 */
template <class ParticleCell>
class KokkosCellPairTraversal : public TraversalInterface {
 public:
  /**
   * Constructor of CellPairTraversal.
   * @param dims the dimensions of the cellblock.
   */
  explicit KokkosCellPairTraversal(const std::array<unsigned long, 3> &dims) : _cellsPerDimension(dims) {}

  /**
   * Destructor of CellPairTraversal.
   */
  ~KokkosCellPairTraversal() override = default;

  /**
   * Sets the cells to iterate over. Should always be called before initTraversal().
   * @param cells The cells to iterate over.
   */
  virtual void setCellsToTraverse(Kokkos::View<ParticleCell *> &c) { this->_cells = c; }

  ParticleCell getCell(size_t index) { return _cells(index); }

 protected:
  /**
   * The dimensions of the cellblock.
   * The dimensions are the number of cells in x, y and z direction.
   */
  std::array<unsigned long, 3> _cellsPerDimension;

  /**
   * The cells to traverse.
   */
  Kokkos::View<ParticleCell *> _cells;
};

}  // namespace autopas
