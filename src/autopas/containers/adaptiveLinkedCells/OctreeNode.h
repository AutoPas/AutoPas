/**
 * @file OctreeNode.h
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include <array>
#include <memory>
#include <optional>

namespace autopas {
namespace internal {

/**
 * Class representing a single node in an octree.
 */
template <class Particle, class ParticleCell>
class OctreeNode {
 public:
  /**
   * Constructor for OctreeNode
   * @param level
   * @param index
   */
  OctreeNode(const unsigned int level, const unsigned int index) : _level(level), _index(index) {}

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
  virtual OctreeNode<Particle, ParticleCell> *update(std::vector<ParticleCell> &cells) = 0;

  /**
   * Returns whether update() would change the tree.
   * @return
   */
  virtual bool isUpdateNeeded() const = 0;

  /**
   * Return the base index.
   * @return
   */
  virtual size_t getIndex() const { return _index; }

 protected:
  const unsigned int _level;
  const unsigned int _index;
};

}  // namespace internal
}  // namespace autopas