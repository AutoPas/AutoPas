/**
 * @file OctreeNodeInterface.h
 *
 * @author Johannes Spies
 * @date 15.04.2021
 */
#pragma once

#include <array>
#include <set>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/octree/OctreeDirection.h"
#include "autopas/utils/inBox.h"

namespace autopas {
/**
 * This forward declaration is required since the `OctreeNodeInterface` provides a method to gather all
 * `OctreeLeafNode`s.
 * @tparam Particle
 */
template <typename Particle>
class OctreeLeafNode;

/**
 * The base class that provides the necessary function definitions that can be applied to an octree.
 *
 * @tparam Particle
 */
template <class Particle>
class OctreeNodeInterface {
 public:
  /**
   * Create an octree node interface by initializing the given fields.
   * @param boxMin The minimum coordinate of the enclosing box
   * @param boxMax The maximum coordinate of the enclosing box
   * @param parent A pointer to the parent of this node, may be nullptr for the root.
   * @param treeSplitThreshold Maximum number of particles inside a leaf before it tries to split itself
   * @param interactionLength The minimum distance at which a force is considered nonzero, cutoff+skin.
   * @param cellSizeFactor The cell size factor
   */
  OctreeNodeInterface(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                      OctreeNodeInterface<Particle> *parent, const int unsigned treeSplitThreshold,
                      const double interactionLength, const double cellSizeFactor)
      : _boxMin(boxMin),
        _boxMax(boxMax),
        _parent(parent),
        _treeSplitThreshold(treeSplitThreshold),
        _interactionLength(interactionLength),
        _cellSizeFactor(cellSizeFactor) {}

  /** To make clang happy. */
  virtual ~OctreeNodeInterface() = default;

  /** Default copy constructor */
  OctreeNodeInterface(const OctreeNodeInterface<Particle> &) = default;

  /**
   * Insert a particle into the octree.
   * @param p The particle to insert
   * @return A std::unique_ptr to a newly created subtree or nullptr if the subtree did not change
   */
  virtual std::unique_ptr<OctreeNodeInterface<Particle>> insert(const Particle &p) = 0;

  /**
   * Put all particles that are below this node into the vector.
   * @param ps A reference to the vector that should contain the particles after the operation
   */
  virtual void collectAllParticles(std::vector<Particle *> &ps) const = 0;

  /**
   * Put the min/max corner coordinates of every leaf into the vector.
   * @param boxes A reference to the vector that should contain pairs of the min/max corner coordinates
   */
  virtual void appendAllLeafBoxes(
      std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> &boxes) const = 0;

  /**
   * Put all leaves below this subtree into a given list.
   * @param leaves A reference to the vector that should contain pointers to the leaves
   */
  virtual void appendAllLeaves(std::vector<OctreeLeafNode<Particle> *> &leaves) const = 0;

  /**
   * Delete the entire tree below this node.
   * @param ref A reference that contains this node
   */
  virtual void clearChildren(std::unique_ptr<OctreeNodeInterface<Particle>> &ref) = 0;

  /**
   * @copydoc CellBasedParticleContainer::getNumberOfParticles()
   */
  virtual unsigned int getNumberOfParticles() = 0;

  /**
   * Get a child node of this node (if there are children) given a specific octant using the spacial structure of the
   * stored children.
   * @param O The octant
   * @return A pointer to a child node
   */
  virtual OctreeNodeInterface<Particle> *SON(octree::Octant O) = 0;

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

  /**
   * Find all leaves below this subtree that are in the given range.
   * @param min The minimum coordinate in 3D space of the query area
   * @param max The maximum coordinate in 3D space of the query area
   * @return A set of all leaf nodes that are in the query region
   */
  virtual std::set<OctreeLeafNode<Particle> *> getLeavesInRange(std::array<double, 3> min,
                                                                std::array<double, 3> max) = 0;

  /**
   * Check if a 3d point is inside the node's axis aligned bounding box. (Set by the boxMin and boxMax fields.)
   * @param point The node to test
   * @return true if the point is inside the node's bounding box and false otherwise
   */
  bool isInside(const std::array<double, 3> &point) {
    using namespace autopas::utils;
    return inBox(point, _boxMin, _boxMax);
  }

  /**
   * Check if the volume enclosed by two boxes a and b is nonzero on a specific axis.
   * @param axis An axis index (0, 1 or 2)
   * @param aMin The minimum coordinate of a's volume
   * @param aMax The maximum coordinate of a's volume
   * @param bMin The minimum coordinate of b's volume
   * @param bMax The maximum coordinate of b's volume
   * @return true iff the enclosed volume is greater than zero, false if the enclosed volume is equal to zero.
   */
  static bool volumeExistsOnAxis(const int axis, const std::array<double, 3> &aMin, const std::array<double, 3> &aMax,
                                 const std::array<double, 3> &bMin, const std::array<double, 3> &bMax) {
    const bool o1 = aMin[axis] < bMax[axis];
    const bool o2 = bMin[axis] < aMax[axis];
    return o1 and o2;
  }

  /**
   * Check if an octree node's box encloses volume with another octree node's box on a specific axis.
   * @param axis An axis index (0, 1 or 2)
   * @param other The octree node to check against.
   * @return true iff the enclosed volume is greater than zero, false if the enclosed volume is equal to zero.
   */
  bool enclosesVolumeWithOtherOnAxis(const int axis, const OctreeNodeInterface<Particle> *other) {
    return volumeExistsOnAxis(axis, this->getBoxMin(), this->getBoxMax(), other->getBoxMin(), other->getBoxMax());
  }

  /**
   * Check if the node's axis aligned bounding box overlaps with the given axis aligned bounding box.
   * @param otherMin The minimum coordinate of the other box
   * @param otherMax The maximum coordinate of the other box
   * @return true iff the overlapping volume is non-negative
   */
  bool overlapsBox(const std::array<double, 3> &otherMin, const std::array<double, 3> &otherMax) {
    bool result = true;
    for (auto d = 0; d < 3; ++d) {
      result &= (this->_boxMin[d] <= otherMax[d]) and (this->_boxMax[d] >= otherMin[d]);
    }
    return result;
  }

  /**
   * Get the enclosed volume between two boxes a and b.
   * @param aMin The minimum coordinate of box a
   * @param aMax The maximum coordinate of box b
   * @param bMin The minimum coordinate of box a
   * @param bMax The maximum coordinate of box b
   * @return The enclosed volume or zero if the boxes do not overlap
   */
  static double getEnclosedVolumeWith(const std::array<double, 3> &aMin, const std::array<double, 3> &aMax,
                                      const std::array<double, 3> &bMin, const std::array<double, 3> &bMax) {
    auto product = 1.0;
    int count = 0;
    for (auto d = 0; d < 3; ++d) {
      auto minOnAxis = std::max(aMin[d], bMin[d]);
      auto maxOnAxis = std::min(aMax[d], bMax[d]);
      auto dim = (maxOnAxis - minOnAxis);
      if (dim > 0) {
        ++count;
      } else {
        break;
      }
      product *= dim;
    }
    auto result = count == 3 ? product : 0;
    return result;
  }

  /**
   * Calculate the overlap volume between the node's axis aligned bounding box and the given box.
   * @param otherMin The minimum coordinate of the other box
   * @param otherMax The maximum coordinate of the other box
   * @return The volume enclosed by the two boxes
   */
  double getEnclosedVolumeWith(const std::array<double, 3> &otherMin, const std::array<double, 3> &otherMax) {
    return getEnclosedVolumeWith(this->getBoxMin(), this->getBoxMax(), otherMin, otherMax);
  }

  /**
   * Find a node (via the pointer structure) that is of equal size of the current node's bounding box, according to the
   * Samet paper. This function searches along a face.
   * @param I The face in which direction the search should find a node
   * @return An octree node
   */
  OctreeNodeInterface<Particle> *EQ_FACE_NEIGHBOR(const octree::Face I) {
    OctreeNodeInterface<Particle> *param, *P = this;
    if (ADJ(I, SONTYPE(P))) {
      param = FATHER(P)->EQ_FACE_NEIGHBOR(I);
    } else {
      param = FATHER(P);
    }
    return param->SON(REFLECT(I, SONTYPE(P)));
  }

  /**
   * Find a node (via the pointer structure) that is of equal size of the current node's bounding box, according to the
   * Samet paper. This function searches along an edge.
   * @param I The edge in which direction the search should find a node
   * @return An octree node
   */
  OctreeNodeInterface<Particle> *EQ_EDGE_NEIGHBOR(const octree::Edge I) {
    OctreeNodeInterface<Particle> *param, *P = this;
    if (ADJ(I, SONTYPE(P))) {
      param = FATHER(P)->EQ_EDGE_NEIGHBOR(I);
    } else if (COMMON_FACE(I, SONTYPE(P)) != octree::O) {
      param = FATHER(P)->EQ_FACE_NEIGHBOR(COMMON_FACE(I, SONTYPE(P)));
    } else {
      param = FATHER(P);
    }
    return param->SON(REFLECT(I, SONTYPE(P)));
  }

  /**
   * Find a node (via the pointer structure) that is of equal size of the current node's bounding box, according to the
   * Samet paper. This function searches along a vertex.
   * @param I The face in which direction the search should find a node
   * @return An octree node
   */
  OctreeNodeInterface<Particle> *EQ_VERTEX_NEIGHBOR(const octree::Vertex I) {
    OctreeNodeInterface<Particle> *param, *P = this;
    if (ADJ(I, SONTYPE(P))) {
      param = FATHER(P)->EQ_VERTEX_NEIGHBOR(I);
    } else if (COMMON_EDGE(I, SONTYPE(P)) != octree::OO) {
      param = FATHER(P)->EQ_EDGE_NEIGHBOR(COMMON_EDGE(I, SONTYPE(P)));
    } else if (COMMON_FACE(I, SONTYPE(P)) != octree::O) {
      param = FATHER(P)->EQ_FACE_NEIGHBOR(COMMON_FACE(I, SONTYPE(P)));
    } else {
      param = FATHER(P);
    }
    return param->SON(REFLECT(I, SONTYPE(P)));
  }

  /**
   * Find a node (via the pointer structure) that is of greater than or equal to the size of the current node's bounding
   * box, according to the Samet paper. This function searches along a face.
   * @param I The face in which direction the search should find a node
   * @return An octree node
   */
  OctreeNodeInterface<Particle> *GTEQ_FACE_NEIGHBOR(octree::Face I);

  /**
   * Find a node (via the pointer structure) that is of greater than or equal to the size of the current node's bounding
   * box, according to the Samet paper. This function searches along an edge.
   * @param I The edge in which direction the search should find a node
   * @return An octree node
   */
  OctreeNodeInterface<Particle> *GTEQ_EDGE_NEIGHBOR(octree::Edge I);

  /**
   * Find a node (via the pointer structure) that is of greater than or equal to the size of the current node's bounding
   * box, according to the Samet paper. This function searches along a vertex.
   * @param I The vertex in which direction the search should find a node
   * @return An octree node
   */
  OctreeNodeInterface<Particle> *GTEQ_VERTEX_NEIGHBOR(octree::Vertex I);

  /**
   * Find all leaf nodes along a list of given directions.
   * @param directions A list of allowed directions for traversal.
   * @return A list of leaf nodes
   */
  virtual std::vector<OctreeLeafNode<Particle> *> getLeavesFromDirections(
      const std::vector<octree::Vertex> &directions) = 0;

  /**
   * This function combines all required functions when traversing down a subtree of the octree and finding all leaves.
   * @param direction The "original" direction. The leaves will be found along the opposite direction.
   * @return A list of leaf nodes
   */
  std::vector<OctreeLeafNode<Particle> *> getNeighborLeaves(const octree::Any direction) {
    auto opposite = octree::getOppositeDirection(direction);
    auto directions = octree::getAllowedDirections(opposite);
    auto neighborLeaves = getLeavesFromDirections(directions);
    return neighborLeaves;
  }

  /**
   * Get the neighbor leaves in all directions
   * @return A set of (unique) neighboring leaves
   */
  std::set<OctreeLeafNode<Particle> *> getNeighborLeaves() {
    std::set<OctreeLeafNode<Particle> *> result;

    // Get all face neighbors
    for (auto face : octree::Tables::faces) {
      OctreeNodeInterface<Particle> *neighbor = GTEQ_FACE_NEIGHBOR(face);
      if (neighbor) {
        auto leaves = neighbor->getNeighborLeaves(face);
        result.insert(leaves.begin(), leaves.end());
      }
    }

    // Get all edge neighbors
    for (auto edge : octree::Tables::edges) {
      OctreeNodeInterface<Particle> *neighbor = GTEQ_EDGE_NEIGHBOR(edge);
      if (neighbor) {
        auto leaves = neighbor->getNeighborLeaves(edge);
        result.insert(leaves.begin(), leaves.end());
      }
    }

    // Get all face neighbors
    for (auto vertex : octree::Tables::vertices) {
      OctreeNodeInterface<Particle> *neighbor = GTEQ_VERTEX_NEIGHBOR(vertex);
      if (neighbor) {
        auto leaves = neighbor->getNeighborLeaves(vertex);
        result.insert(leaves.begin(), leaves.end());
      }
    }

    return result;
  }

  /**
   * Set the minimum coordinate of the enclosing box.
   * @param boxMin A point in 3D space
   */
  void setBoxMin(const std::array<double, 3> &boxMin) { _boxMin = boxMin; }

  /**
   * Set the maximum coordinate of the enclosing box.
   * @param boxMax A point in 3D space
   */
  void setBoxMax(const std::array<double, 3> &boxMax) { _boxMax = boxMax; }

  /**
   * Get the minimum coordinate of the enclosing box.
   * @return A point in 3D space
   */
  [[nodiscard]] std::array<double, 3> getBoxMin() const { return _boxMin; }

  /**
   * Get the maximum coordinate of the enclosing box.
   * @return A point in 3D space
   */
  [[nodiscard]] std::array<double, 3> getBoxMax() const { return _boxMax; }

  /**
   * Get the parent node of this node.
   * @return A pointer to the parent, can be nullptr.
   */
  OctreeNodeInterface<Particle> *getParent() const { return _parent; }

 protected:
  /**
   * Check if this is not the root node.
   * @return true iff there does not exist a parent for this node.
   */
  bool hasParent() { return _parent != nullptr; }

  /**
   * A pointer to the parent node. Can be nullptr (iff this is the root node).
   */
  OctreeNodeInterface<Particle> *_parent;

  /**
   * The min coordinate of the enclosed volume.
   */
  std::array<double, 3> _boxMin;

  /**
   * The max coordinate of the enclosed volume.
   */
  std::array<double, 3> _boxMax;

  /**
   * Maximum number of particles inside a leaf node before the leaf tries to split itself.
   */
  int unsigned _treeSplitThreshold;

  /**
   * The minimum distance at which a force is considered nonzero, cutoff+skin.
   */
  double _interactionLength;

  /**
   * The cell size factor for this node
   */
  double _cellSizeFactor;
};

/**
 * Check if a node is an inner node.
 * @tparam Particle The particle type used in this container
 * @param node A pointer to a node
 * @return true if the node has children, false otherwise.
 */
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
inline OctreeNodeInterface<Particle> *FATHER(const OctreeNodeInterface<Particle> *node) {
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
static octree::Octant SONTYPE(const OctreeNodeInterface<Particle> *node) {
  octree::Octant result = octree::OOO;
  if (FATHER(node)) {
    for (octree::Vertex test : octree::Tables::vertices) {
      if (FATHER(node)->SON(test) == node) {
        result = test;
        break;
      }
    }
    if (result == octree::OOO) {
      throw std::runtime_error("[OctreeNodeInterface::SONTYPE()] Unable to determine SONTYPE");
    }
  }
  return result;
}

template <class Particle>
OctreeNodeInterface<Particle> *OctreeNodeInterface<Particle>::GTEQ_FACE_NEIGHBOR(const octree::Face I) {
  // Check precondition
  if (not isFace(I)) {
    throw std::runtime_error("[OctreeNodeInterface::GTEQ_FACE_NEIGHBOR()] Received invalid face.");
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
OctreeNodeInterface<Particle> *OctreeNodeInterface<Particle>::GTEQ_EDGE_NEIGHBOR(const octree::Edge I) {
  // Check precondition
  if (not isEdge(I)) {
    throw std::runtime_error("[OctreeNodeInterface::GTEQ_EDGE_NEIGHBOR()] Received invalid edge.");
  }

  auto null = [](OctreeNodeInterface<Particle> *T) { return T == nullptr; };

  // Find a common ancestor
  OctreeNodeInterface<Particle> *Q, *P = this;
  if (null(FATHER(P))) {
    Q = nullptr;
  } else if (ADJ(I, SONTYPE(P))) {
    Q = FATHER(P)->GTEQ_EDGE_NEIGHBOR(I);
  } else if (octree::Face common = COMMON_FACE(I, SONTYPE(P)); common != octree::O) {
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
OctreeNodeInterface<Particle> *OctreeNodeInterface<Particle>::GTEQ_VERTEX_NEIGHBOR(const octree::Vertex I) {
  // Check precondition
  if (not isVertex(I)) {
    throw std::runtime_error("[OctreeNodeInterface::GTEQ_VERTEX_NEIGHBOR()] Received invalid vertex.");
  }

  // Find a common ancestor
  OctreeNodeInterface<Particle> *Q, *P = this;
  if (not FATHER(P)) {
    Q = nullptr;
  } else if (ADJ(I, SONTYPE(P))) {
    Q = FATHER(P)->GTEQ_VERTEX_NEIGHBOR(I);
  } else if (octree::Edge commonEdge = COMMON_EDGE(I, SONTYPE(P)); commonEdge != octree::OO) {
    Q = FATHER(P)->GTEQ_EDGE_NEIGHBOR(commonEdge);
  } else if (octree::Face commonFace = COMMON_FACE(I, SONTYPE(P)); commonFace != octree::O) {
    Q = FATHER(P)->GTEQ_FACE_NEIGHBOR(commonFace);
  } else {
    Q = FATHER(P);
  }

  if (Q and GRAY(Q)) {
    // Follow opposite path to locate the neighbor
    return Q->SON(REFLECT(I, SONTYPE(P)));
  } else {
    return Q;
  }
}
}  // namespace autopas
