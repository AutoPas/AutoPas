/**
 * @file OctreeNodeInterface.h
 *
 * @author Johannes Spies
 * @date 15.04.2021
 */
#pragma once

#include <array>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/octree/OctreeDirection.h"
#include "autopas/utils/inBox.h"

namespace autopas {
template <typename Particle>
class OctreeLeafNode;
template <typename Particle>
class OctreeInnerNode;

/**
 * The base class that provides the necessary function definitions that can be applied to an octree.
 *
 * @tparam Particle
 */
template <class Particle>
class OctreeNodeInterface {
 public:
  OctreeNodeInterface(std::array<double, 3> boxMin, std::array<double, 3> boxMax, OctreeNodeInterface<Particle> *parent)
      : _boxMin(boxMin), _boxMax(boxMax), _parent(parent) {}

  /**
   * Insert a particle into the octree.
   * @param p The particle to insert
   * @return The subtree below the current node that now contains the particle
   */
  virtual void insert(std::unique_ptr<OctreeNodeInterface<Particle>> &ref, Particle p) = 0;

  /**
   * Put all particles that are below this node into the vector.
   * @param ps A reference to the vector that should contain the particles after the operation
   */
  virtual void appendAllParticles(std::vector<Particle> &ps) = 0;

  /**
   * Put the min/max corner coordinates of every leaf into the vector.
   * @param boxes A reference to the vector that should contain pairs of the min/max corner coordinates
   */
  virtual void appendAllLeafBoxes(std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> &boxes) = 0;

  /**
   * Delete the entire tree below this node.
   */
  virtual void clearChildren(std::unique_ptr<OctreeNodeInterface<Particle>> &ref) = 0;

  /**
   * @copydoc CellBasedParticleContainer::getNumParticles()
   */
  virtual unsigned int getNumParticles() = 0;

  /**
   * Check if the node is a leaf or an inner node. The function exists for debugging.
   * @return true iff the node is a leaf, false otherwise.
   */
  virtual bool hasChildren() = 0;

  /**
   * Get a child by its index from the node.
   * @param index The index of the child. Must be between 0 and 7 inclusive.
   * @return A pointer to the child.
   */
  virtual OctreeNodeInterface<Particle> *getChild(int index) = 0;

  virtual std::optional<OctreeNodeInterface<Particle> *> getGreaterParentAlongAxis(
      int axis, int dir, OctreeNodeInterface<Particle> *embedded) = 0;

  virtual std::vector<OctreeNodeInterface<Particle> *> findTouchingLeaves(int axis, int dir,
                                                                          OctreeNodeInterface<Particle> *embedded) = 0;

  bool isPositive(int dir) {
    if (dir == 0) {
      throw std::runtime_error("[OctreeNodeInterface.h] Unable to obtain direction from 0");
    }
    return dir > 0;
  }

  bool isNegative(int dir) {
    if (dir == 0) {
      throw std::runtime_error("[OctreeNodeInterface.h] Unable to obtain direction from 0");
    }
    return dir < 0;
  }

  static int oppositeDirection(int dir) {
    switch(dir) {
      case 1: return -1;
      case -1: return 1;
      default: throw std::runtime_error("[OctreeNodeInterface.h] dir must be -1 or 1.");
    }
  }

  std::vector<OctreeNodeInterface<Particle> *> getNeighborsAlongAxis(int axis, int dir) {
    std::vector<OctreeNodeInterface<Particle> *> result = {};
    if(auto greaterParentOptional = getGreaterParentAlongAxis(axis, dir, this)) {
      OctreeNodeInterface<Particle> *greaterParent = *greaterParentOptional;
      if(!greaterParent->hasChildren()) {
        throw std::runtime_error("[OctreeNodeInterface.h] The greater parent must not be a leaf if it exists.");
      }

      for(int i = 0; i < 8; ++i) {
        auto *child = greaterParent->getChild(i);
        if(enclosesVolumeWithOtherOnAxis(axis, child)) {
          int downDir = oppositeDirection(dir);
          result = child->findTouchingLeaves(axis, downDir, this);
          break;
        }
      }
    }
    return result;
  }

  std::vector<OctreeNodeInterface<Particle> *> getNeighbors() {
    std::vector<OctreeNodeInterface<Particle> *> result = {};
    for (int axis = 0; axis < 3; ++axis) {
      for (int dir = -1; dir <= 1; dir += 2) {
        auto neighborsAlongAxis = getNeighborsAlongAxis(axis, dir);
        result.insert(result.end(), neighborsAlongAxis.begin(), neighborsAlongAxis.end());
      }
    }
    return result;
  }

  /**
   * Check if a 3d point is inside the node's axis aligned bounding box. (Set by the boxMin and boxMax fields.)
   * @param node The possible enclosing node
   * @param point The node to test
   * @return true if the point is inside the node's bounding box and false otherwise
   */
  bool isInside(std::array<double, 3> point) {
    using namespace autopas::utils;
    return inBox(point, _boxMin, _boxMax);
  }

  static bool volumeExistsOnAxis(int axis, std::array<double, 3> aMin, std::array<double, 3> aMax,
                                 std::array<double, 3> bMin, std::array<double, 3> bMax) {
    bool o1 = aMin[axis] < bMax[axis];
    bool o2 = bMin[axis] < aMax[axis];
    return o1 && o2;
  }

  bool enclosesVolumeWithOtherOnAxis(int axis, OctreeNodeInterface<Particle> *other) {
    return volumeExistsOnAxis(axis, this->getBoxMin(), this->getBoxMax(), other->getBoxMin(), other->getBoxMax());
  }

  /**
   * Check if the node's axis aligned bounding box overlaps with the given axis aligned bounding box.
   * @param otherMin The minimum coordinate of the other box
   * @param otherMax The maximum coordinate of the other box
   * @return true iff the overlapping volume is non-negative
   */
  bool overlapsBox(std::array<double, 3> otherMin, std::array<double, 3> otherMax) {
    bool result = true;
    for (auto d = 0; d < 3; ++d) {
      result &= (this->_boxMin[d] <= otherMax[d]) && (this->_boxMax[d] >= otherMin[d]);
    }
    return result;
  }

  void setBoxMin(std::array<double, 3> boxMin) { _boxMin = boxMin; }
  void setBoxMax(std::array<double, 3> boxMax) { _boxMax = boxMax; }
  std::array<double, 3> getBoxMin() { return _boxMin; }
  std::array<double, 3> getBoxMax() { return _boxMax; }

 protected:
  bool hasParent() { return _parent != nullptr; }

  OctreeNodeInterface<Particle> *_parent;
  std::array<double, 3> _boxMin, _boxMax;
};
}  // namespace autopas
