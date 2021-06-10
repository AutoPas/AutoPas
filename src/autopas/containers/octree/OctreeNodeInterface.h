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

  // To make clang happy.
  virtual ~OctreeNodeInterface() = default;

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

  virtual void appendAllLeaves(std::vector<OctreeLeafNode<Particle> *> &leaves) = 0;

  /**
   * Delete the entire tree below this node.
   */
  virtual void clearChildren(std::unique_ptr<OctreeNodeInterface<Particle>> &ref) = 0;

  /**
   * @copydoc CellBasedParticleContainer::getNumParticles()
   */
  virtual unsigned int getNumParticles() = 0;

  virtual OctreeNodeInterface<Particle> *SON(Octant O) = 0;

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
    switch (dir) {
      case 1:
        return -1;
      case -1:
        return 1;
      default:
        throw std::runtime_error("[OctreeNodeInterface.h] dir must be -1 or 1.");
    }
  }

  std::vector<OctreeNodeInterface<Particle> *> getNeighborsAlongAxis(int axis, int dir) {
    std::vector<OctreeNodeInterface<Particle> *> result = {};
    if (auto greaterParentOptional = getGreaterParentAlongAxis(axis, dir, this)) {
      OctreeNodeInterface<Particle> *greaterParent = *greaterParentOptional;
      if (!greaterParent->hasChildren()) {
        throw std::runtime_error("[OctreeNodeInterface.h] The greater parent must not be a leaf if it exists.");
      }

      for (int i = 0; i < 8; ++i) {
        auto *child = greaterParent->getChild(i);
        if (enclosesVolumeWithOtherOnAxis(axis, child)) {
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

  OctreeNodeInterface<Particle> *EQ_FACE_NEIGHBOR(Face I) {
    OctreeNodeInterface<Particle> *param, *P = this;
    if (ADJ(I, SONTYPE(P))) {
      param = FATHER(P)->EQ_FACE_NEIGHBOR(I);
    } else {
      param = FATHER(P);
    }
    return param->SON(REFLECT(I, SONTYPE(P)));
  }

  OctreeNodeInterface<Particle> *EQ_EDGE_NEIGHBOR(Edge I) {
    OctreeNodeInterface<Particle> *param, *P = this;
    if (ADJ(I, SONTYPE(P))) {
      param = FATHER(P)->EQ_EDGE_NEIGHBOR(I);
    } else if (COMMON_FACE(I, SONTYPE(P)) != O) {
      param = FATHER(P)->EQ_FACE_NEIGHBOR(COMMON_FACE(I, SONTYPE(P)));
    } else {
      param = FATHER(P);
    }
    return param->SON(REFLECT(I, SONTYPE(P)));
  }

  OctreeNodeInterface<Particle> *EQ_VERTEX_NEIGHBOR(Vertex I) {
    OctreeNodeInterface<Particle> *param, *P = this;
    if (ADJ(I, SONTYPE(P))) {
      param = FATHER(P)->EQ_VERTEX_NEIGHBOR(I);
    } else if (COMMON_EDGE(I, SONTYPE(P)) != OO) {
      param = FATHER(P)->EQ_EDGE_NEIGHBOR(COMMON_EDGE(I, SONTYPE(P)));
    } else if (COMMON_FACE(I, SONTYPE(P)) != O) {
      param = FATHER(P)->EQ_FACE_NEIGHBOR(COMMON_FACE(I, SONTYPE(P)));
    } else {
      param = FATHER(P);
    }
    return param->SON(REFLECT(I, SONTYPE(P)));
  }

  OctreeNodeInterface<Particle> *GTEQ_FACE_NEIGHBOR(Face I);
  OctreeNodeInterface<Particle> *GTEQ_EDGE_NEIGHBOR(Edge I);
  OctreeNodeInterface<Particle> *GTEQ_VERTEX_NEIGHBOR(Vertex I);

  std::vector<OctreeNodeInterface<Particle> *> getLeavesFromDirections(std::vector<Vertex> directions) {
    std::vector<OctreeNodeInterface<Particle> *> result;
    // TODO: Make this virtual and handled by the concrete subclasses
    if (hasChildren()) {
      for (auto d : directions) {
        int childIndex = vertexToIndex(d);
        auto child = getChild(childIndex);
        auto childLeaves = child->getLeavesFromDirections(directions);
        result.insert(result.end(), childLeaves.begin(), childLeaves.end());
      }
    } else {
      result.push_back(this);
    }
    return result;
  }

  std::vector<OctreeNodeInterface<Particle> *> getNeighborLeaves(Any direction) {
    auto opposite = getOppositeDirection(direction);
    auto directions = getAllowedDirections(opposite);
    auto neighborLeaves = getLeavesFromDirections(directions);
    return neighborLeaves;
  }

  void setBoxMin(std::array<double, 3> boxMin) { _boxMin = boxMin; }
  void setBoxMax(std::array<double, 3> boxMax) { _boxMax = boxMax; }
  std::array<double, 3> getBoxMin() { return _boxMin; }
  std::array<double, 3> getBoxMax() { return _boxMax; }

  OctreeNodeInterface<Particle> *getParent() { return _parent; }

 protected:
  bool hasParent() { return _parent != nullptr; }

  OctreeNodeInterface<Particle> *_parent;
  std::array<double, 3> _boxMin, _boxMax;
};

template <class Particle>
inline bool GRAY(OctreeNodeInterface<Particle> *node) {
  // According to Samet: "All non-leaf nodes are said to be GRAY"
  return node->hasChildren();
}

/**
 * Get the parent node of an arbitrary octree node.
 *
 * @tparam Particle
 * @param node
 * @return The parent of the given node if the node is not the root node, otherwise nullptr.
 */
template <class Particle>
inline OctreeNodeInterface<Particle> *FATHER(OctreeNodeInterface<Particle> *node) {
  return node->getParent();
}

/**
 * Get the octant in which a given node can be found in the parent.
 *
 * @tparam Particle
 * @param node The node to check
 * @return The octant in which the node is in the parent, OOO if the node either does not have a parent or could not be
 * found in the parent.
 */
template <class Particle>
static Octant SONTYPE(OctreeNodeInterface<Particle> *node) {
  Octant result = OOO;
  if (FATHER(node)) {
    for (Vertex *test = VERTICES(); *test != OOO; ++test) {
      if (FATHER(node)->SON(*test) == node) {
        result = *test;
        break;
      }
    }
    assert(result != OOO);
  }
  return result;
}

template <class Particle>
OctreeNodeInterface<Particle> *OctreeNodeInterface<Particle>::GTEQ_FACE_NEIGHBOR(Face I) {
  // Check precondition
  if (!contains(getFaces(), O, I)) {
    throw std::runtime_error("[OctreeNodeInterface.h] Received invalid face.");
  }

  auto null = [](OctreeNodeInterface<Particle> *T) { return T == nullptr; };

  // Find a common ancestor
  OctreeNodeInterface<Particle> *Q, *P = this;
  if ((not null(FATHER(P))) and ADJ(I, SONTYPE(P))) {
    Q = FATHER(P)->GTEQ_FACE_NEIGHBOR(I);
  } else {
    Q = FATHER(P);
  }

  if ((not null(Q)) and GRAY(Q)) {
    // Follow the reflected path to locate the neighbor
    return Q->SON(REFLECT(I, SONTYPE(P)));
  } else {
    return Q;
  }
}

template <class Particle>
OctreeNodeInterface<Particle> *OctreeNodeInterface<Particle>::GTEQ_EDGE_NEIGHBOR(Edge I) {
  // Check precondition
  if (!contains(getEdges(), OO, I)) {
    throw std::runtime_error("[OctreeNodeInterface.h] Received invalid edge.");
  }

  auto null = [](OctreeNodeInterface<Particle> *T) { return T == nullptr; };

  // Find a common ancestor
  OctreeNodeInterface<Particle> *Q, *P = this;
  if (null(FATHER(P))) {
    Q = nullptr;
  } else if (ADJ(I, SONTYPE(P))) {
    Q = FATHER(P)->GTEQ_EDGE_NEIGHBOR(I);
  } else if (Face common = COMMON_FACE(I, SONTYPE(P)); common != O) {
    Q = FATHER(P)->GTEQ_FACE_NEIGHBOR(common);
  } else {
    Q = FATHER(P);
  }

  if ((not null(Q)) and GRAY(Q)) {
    // Follow opposite path to locate the neighbor
    return Q->SON(REFLECT(I, SONTYPE(P)));
  } else {
    return Q;
  }
}

template <class Particle>
OctreeNodeInterface<Particle> *OctreeNodeInterface<Particle>::GTEQ_VERTEX_NEIGHBOR(Vertex I) {
  // Check precondition
  if (!contains(VERTICES(), OOO, I)) {
    throw std::runtime_error("[OctreeNodeInterface.h] Received invalid vertex.");
  }

  auto null = [](OctreeNodeInterface<Particle> *T) { return T == nullptr; };

  // Find a common ancestor
  OctreeNodeInterface<Particle> *Q, *P = this;
  if (null(FATHER(P))) {
    Q = nullptr;
  } else if (ADJ(I, SONTYPE(P))) {
    Q = FATHER(P)->GTEQ_VERTEX_NEIGHBOR(I);
  } else if (Edge commonEdge = COMMON_EDGE(I, SONTYPE(P)); commonEdge != OO) {
    Q = FATHER(P)->GTEQ_EDGE_NEIGHBOR(commonEdge);
  } else if (Face commonFace = COMMON_FACE(I, SONTYPE(P)); commonFace != O) {
    Q = FATHER(P)->GTEQ_FACE_NEIGHBOR(commonFace);
  } else {
    Q = FATHER(P);
  }

  if ((not null(Q)) and GRAY(Q)) {
    // Follow opposite path to locate the neighbor
    return Q->SON(REFLECT(I, SONTYPE(P)));
  } else {
    return Q;
  }
}
}  // namespace autopas
