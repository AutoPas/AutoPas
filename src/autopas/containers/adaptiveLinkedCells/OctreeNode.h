/**
 * @file OctreeNode.h
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include <array>
#include <memory>

namespace autopas {

/**
 * Class representing a single node in an octree.
 */
template <class Particle, class ParticleCell>
class OctreeNode {
 public:
  /**
   * Constructor for OctreeNode
   * @param level
   */
  OctreeNode(const unsigned int level) : _level(level) {}

  /**
   * Destructor of OctreeNode.
   */
  virtual ~OctreeNode() = default;

  /**
   * Get the containing cell of a specified position.
   * @param pos the position for which the cell is needed
   * @return cell at the given position
   */
  virtual ParticleCell &getContainingCell(const std::array<double, 3> &pos) const = 0;

  /**
   * Return the number of stored elements.
   * @return
   */
  virtual size_t getSize() const = 0;

  /**
   * Updates the tree.
   * @return new child.
   */
  // virtual std::unique_ptr<OctreeNode<ParticleCell>> update() = 0;

  /**
   * Returns whether update() would change the tree.
   * @return
   */
  virtual bool isUpdateNeeded() const = 0;

 private:
  const unsigned int _level;
  static std::vector<ParticleCell> *_cells;
};

}  // namespace autopas