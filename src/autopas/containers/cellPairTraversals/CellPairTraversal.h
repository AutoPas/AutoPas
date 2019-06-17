/**
 * @file CellPairTraversal.h
 *
 * @date 22 Jan 2018
 * @author tchipevn
 */

#pragma once

#include <array>
#include <vector>
#include "autopas/containers/cellPairTraversals/TraversalInterface.h"

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
   */
  explicit CellPairTraversal(const std::array<unsigned long, 3> &dims) : _cellsPerDimension(dims) {}

  /**
   * Destructor of CellPairTraversal.
   */
  ~CellPairTraversal() override = default;

  /**
   * Resets the cell structure of the traversal.
   * @param dims
   */
  virtual void rebuild(const std::array<unsigned long, 3> &dims) { _cellsPerDimension = dims; };

  /**
   * load Data Layouts required for this Traversal.
   * @param cells where the data should be loaded
   */
  virtual void initTraversal(std::vector<ParticleCell> &cells) = 0;

  /**
   * write Data to AoS.
   * @param cells for which the data should be written back
   */
  virtual void endTraversal(std::vector<ParticleCell> &cells) = 0;

 protected:
  /**
   * The dimensions of the cellblock.
   * The dimensions are the number of cells in x, y and z direction.
   */
  std::array<unsigned long, 3> _cellsPerDimension;
};

}  // namespace autopas
