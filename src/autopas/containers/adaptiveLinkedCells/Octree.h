/**
 * @file Octree.h
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include <array>
#include <memory>

#include "autopas/containers/adaptiveLinkedCells/OctreeExternalNode.h"
#include "autopas/containers/adaptiveLinkedCells/OctreeInternalNode.h"
#include "autopas/utils/ExceptionHandler.h"

namespace autopas {

/**
 * Class representing an octree.
 */
template <class Particle, class ParticleCell>
class Octree {
 public:
  Octree() = default;
  /**
   * Constructor for Octree.
   * @param cells Underlaying cells.
   * @param boxMin
   * @param boxMax
   */
  Octree(std::vector<ParticleCell> &cells, const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax)
      : _boxMin(boxMin), _boxMax(boxMax), _cells(&cells) {}

  /**
   * Destructor of Octree.
   */
  virtual ~Octree() { delete root; }

  /**
   * Get the containing cell of a specified position.
   * @param pos the position for which the cell is needed
   * @return cell at the given position
   */
  ParticleCell &getContainingCell(const std::array<double, 3> &pos) const { return root->getContainingCell(pos); }

  /**
   * Return the number of stored elements.
   * @return Number of elements.
   */
  size_t getSize() const { return root->getSize(); }

  /**
   * Inits tree.
   */
  void init(const std::array<unsigned long, 3> &cellsPerDimension) {
    _cellsPerDimension = cellsPerDimension;
    if (_cellsPerDimension[0] != _cellsPerDimension[1] or _cellsPerDimension[1] != _cellsPerDimension[2]) {
      utils::ExceptionHandler::exception("Octree: Cells per dimension varies");
    }
    if ((_cellsPerDimension[0] & (_cellsPerDimension[0] - 1)) != 0) {
      utils::ExceptionHandler::exception("Octree: Cells per dimension is not a power of 2");
    }

    unsigned int exp = 0ul;
    while (_cellsPerDimension[0] != (1ul << exp)) {
      ++exp;
    }
    root = new internal::OctreeExternalNode<Particle, ParticleCell>(
        *_cells, exp, ArrayMath::mulScalar(ArrayMath::sub(_boxMax, _boxMin), 0.5), 0);
  }

  /**
   * Updates the tree.
   */
  void update() { root = root->update(*_cells); }

  /**
   * Returns whether update() would change the tree.
   * @return true = update needed; false = no update needed
   */
  bool isUpdateNeeded() const { return root->isUpdateNeeded(); }

  void apply(std::function<void(internal::OctreeNode<Particle, ParticleCell> &)> func,
             internal::ExecutionPolicy policy) {
    root->apply(func, policy);
  }

  /**
   * Sets min. number of elements inside of each node.
   * @param minElements Min. number of elements.
   */
  static void setMinElements(const unsigned long minElements) {
    internal::OctreeInternalNode<Particle, ParticleCell>::setMinElements(minElements);
  }

  /**
   * Sets max. number of elements inside of each node.
   * @param maxElements Max. number of elements.
   */
  static void setMaxElements(const unsigned long maxElements) {
    internal::OctreeExternalNode<Particle, ParticleCell>::setMaxElements(maxElements);
  }

 private:
  const std::array<double, 3> _boxMin;
  const std::array<double, 3> _boxMax;

  std::array<unsigned long, 3> _cellsPerDimension;

  internal::OctreeNode<Particle, ParticleCell> *root;
  std::vector<ParticleCell> *_cells;
};

}  // namespace autopas