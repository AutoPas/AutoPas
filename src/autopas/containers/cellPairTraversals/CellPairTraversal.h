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
 * @tparam dataLayout
 * @tparam useNewton3
 */
template <class ParticleCell, DataLayoutOption dataLayout, bool useNewton3>
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
   * Load Data Layouts required for this Traversal.
   * @param Cells where the data should be loaded.
   */
  virtual void initTraversal(std::vector<ParticleCell> &cells) = 0;

  /**
   * Write Data to AoS.
   * @param Cells for which the data should be written back.
   */
  virtual void endTraversal(std::vector<ParticleCell> &cells) = 0;

  bool getUseNewton3() const override { return useNewton3; };

  DataLayoutOption getDataLayout() const override { return dataLayout; };

 protected:
  /**
   * The dimensions of the cellblock.
   * The dimensions are the number of cells in x, y and z direction.
   */
  std::array<unsigned long, 3> _cellsPerDimension;
};

}  // namespace autopas
