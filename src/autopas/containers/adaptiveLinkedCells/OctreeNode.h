/**
 * @file OctreeNode.h
 * @author C.Menges
 * @date 20.06.2019
 */

#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <optional>

#include "autopas/utils/ArrayMath.h"

namespace autopas {
namespace internal {

enum class ExecutionPolicy { seq, par };

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
  OctreeNode(OctreeNode<Particle, ParticleCell> *parent, const unsigned int index, const std::array<double, 3> &center)
      : _parent(parent), _level((parent == nullptr) ? _maxExp : parent->_level - 3), _index(index), _center(center) {}

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

  virtual void apply(std::function<void(OctreeNode<Particle, ParticleCell> &)> func, ExecutionPolicy policy) = 0;

  /**
   * Get a string representation of the OctreeNode
   * @return string representation
   */
  virtual operator std::string() const = 0;

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

  /**
   * Return the base index.
   * @return
   */
  virtual const std::array<double, 3> &getCenter() const { return _center; }

  /**
   * Updates all Neighbors. This is necessary after splitting or combination of nodes.
   */
  virtual void updateNeigbors() = 0;

  /**
   * Return the base index.
   * @return
   */
  virtual size_t getLevel() const { return _level; }

  virtual std::vector<Particle> getOutliers() = 0;

  virtual bool isLeaf() const = 0;

  static void setMaxExp(size_t exp) { _maxExp = exp; }

 protected:
  OctreeNode<Particle, ParticleCell> *_parent;
  const size_t _level;
  const size_t _index;
  const std::array<double, 3> _center;
  static size_t _maxExp;
};

template <class Particle, class ParticleCell>
size_t OctreeNode<Particle, ParticleCell>::_maxExp = 0ul;

}  // namespace internal
}  // namespace autopas