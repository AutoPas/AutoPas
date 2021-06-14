/**
 * @file OctreeInnerNode.h
 *
 * @author Johannes Spies
 * @date 15.04.2021
 */

#pragma once

#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/inBox.h"

namespace autopas {
/**
 * Inner nodes of the octree data structure. An inner node always points to eight children, which can either be leaves
 * or inner nodes as well.
 *
 * @tparam Particle
 */
template <class Particle>
class OctreeInnerNode : public OctreeNodeInterface<Particle> {
 public:
  /**
   * Create an octree inner node that points to eight leaves.
   * @param boxMin The min coordinate of the octree box
   * @param boxMax The max coordinate of the octree box
   * @param parent A pointer to the parent node. Should be nullptr for root nodes.
   */
  OctreeInnerNode(std::array<double, 3> boxMin, std::array<double, 3> boxMax, OctreeNodeInterface<Particle> *parent)
      : OctreeNodeInterface<Particle>(boxMin, boxMax, parent) {
    using namespace autopas::utils;

    // The inner node is initialized with 8 leaves.
    auto center = ArrayMath::mulScalar(ArrayMath::add(boxMin, boxMax), 0.5);
    for (auto i = 0; i < _children.size(); ++i) {
      // Subdivide the bounding box of the parent.
      std::array<double, 3> newBoxMin = {};
      std::array<double, 3> newBoxMax = {};
      for (auto d = 0; d < 3; ++d) {
#if 1
        auto mask = 4 >> d;
#else
        auto mask = 1 << d;
#endif
        newBoxMin[d] = !(i & mask) ? boxMin[d] : center[d];
        newBoxMax[d] = !(i & mask) ? center[d] : boxMax[d];
      }

      // Assign new leaves as the children.
      _children[i] = std::make_unique<OctreeLeafNode<Particle>>(newBoxMin, newBoxMax, this);
    }
  }

  /**
   * @copydoc OctreeNodeInterface::insert()
   */
  void insert(std::unique_ptr<OctreeNodeInterface<Particle>> &ref, Particle p) override {
    if (!this->isInside(p.getR())) {
      throw std::runtime_error("[OctreeInnerNode.h] Attempting to insert particle that is not inside this node");
    }

    // Find a child to insert the particle into.
    for (auto &child : _children) {
      if (child->isInside(p.getR())) {
        child->insert(child, p);
        break;
      }
    }
  }

  /**
   * @copydoc OctreeNodeInterface::appendAllParticles()
   */
  void appendAllParticles(std::vector<Particle> &ps) override {
    // An inner node does not contain particles, traverse down to the children.
    for (auto &child : _children) {
      child->appendAllParticles(ps);
    }
  }

  /**
   * @copydoc OctreeNodeInterface::appendAllLeafBoxes()
   */
  void appendAllLeafBoxes(std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> &boxes) override {
    for (auto &child : _children) {
      child->appendAllLeafBoxes(boxes);
    }
  }

  /**
   * @copydoc OctreeNodeInterface::clearChildren()
   */
  void clearChildren(std::unique_ptr<OctreeNodeInterface<Particle>> &ref) override {
    for (auto &child : _children) {
      child->clearChildren(child);
    }

    std::unique_ptr<OctreeLeafNode<Particle>> newLeaf =
        std::make_unique<OctreeLeafNode<Particle>>(this->getBoxMin(), this->getBoxMax(), this->_parent);
    ref = std::move(newLeaf);
  }

  /**
   * @copydoc OctreeNodeInterface::getNumParticles()
   */
  unsigned int getNumParticles() override {
    unsigned int result = 0;
    for (auto &child : _children) {
      result += child->getNumParticles();
    }
    return result;
  }

  /**
   * @copydoc OctreeNodeInterface::hasChildren()
   */
  bool hasChildren() override { return true; }

  /**
   * @copydoc OctreeNodeInterface::getChild()
   */
  OctreeNodeInterface<Particle> *getChild(int index) override { return _children[index].get(); }

#if 0
  std::optional<OctreeNodeInterface<Particle> *> getGreaterParentAlongAxis(
      int axis, int dir, OctreeNodeInterface<Particle> *embedded) override {
    auto posDir = (dir > 0);
    auto negDir = (dir < 0);
    if (!posDir && !negDir) {
      throw std::runtime_error("[OctreeInnerNode] dir must be -1 or 1.");
    }

    // TODO(johannes): Check if those two hold for every case.
    auto smallerOnAxis = (embedded->getBoxMax()[axis] < this->getBoxMax()[axis]);
    auto biggerOnAxis = (embedded->getBoxMin()[axis] > this->getBoxMin()[axis]);

    std::optional<OctreeNodeInterface<Particle> *> result = std::nullopt;
    if ((posDir && smallerOnAxis) || (negDir && biggerOnAxis)) {
      result = std::make_optional(this);
    } else if (this->hasParent()) {
      // See if this node's parent includes a neighbor
      result = this->_parent->getGreaterParentAlongAxis(axis, dir, embedded);
    }
    return result;
  }

  bool isNearestInDirection(int childIndex, int axis, int dir) {
    int dirFlag = (childIndex >> axis) & 1;
    return (dirFlag == 0) == this->isNegative(dir);
  }

  std::vector<OctreeNodeInterface<Particle> *> findTouchingLeaves(int axis, int dir,
                                                                  OctreeNodeInterface<Particle> *embedded) override {
    std::vector<OctreeNodeInterface<Particle> *> result;
    for (int i = 0; i < _children.size(); ++i) {
      auto &child = _children[i];

      // Check if the child is in the right search direction
      if (isNearestInDirection(i, axis, dir)) {
        // Check if the child overlaps with the embedded node on the axis that are not the search axis
        int otherAxis1 = (axis + 1) % 3;
        bool volumeOnAxis1 = OctreeNodeInterface<Particle>::volumeExistsOnAxis(
            otherAxis1, child->getBoxMin(), child->getBoxMax(), embedded->getBoxMin(), embedded->getBoxMax());

        int otherAxis2 = (axis + 2) % 3;
        bool volumeOnAxis2 = OctreeNodeInterface<Particle>::volumeExistsOnAxis(
            otherAxis2, child->getBoxMin(), child->getBoxMax(), embedded->getBoxMin(), embedded->getBoxMax());

        if (volumeOnAxis1 && volumeOnAxis2) {
          auto otherTouching = child->findTouchingLeaves(axis, dir, embedded);
          result.insert(result.end(), otherTouching.begin(), otherTouching.end());
        }
      }
    }

    return result;
  }
#endif

  std::vector<OctreeNodeInterface<Particle> *> getLeavesFromDirections(std::vector<Vertex> directions) override {
    std::vector<OctreeNodeInterface<Particle> *> result;
    // Only take the children that are allowed (i.e. those which are in the given directions list)
    for (auto d : directions) {
      int childIndex = vertexToIndex(d);
      auto child = getChild(childIndex);
      auto childLeaves = child->getLeavesFromDirections(directions);
      result.insert(result.end(), childLeaves.begin(), childLeaves.end());
    }
    return result;
  }

  OctreeNodeInterface<Particle> *SON(Octant octant) override {
    // convert the Octant to a flat child index
    auto flat = vertexToIndex(octant);
    return _children[flat].get();
  }

  void appendAllLeaves(std::vector<OctreeLeafNode<Particle> *> &leaves) override {
    for (auto &child : _children) {
      child->appendAllLeaves(leaves);
    }
  }

 private:
  /**
   * Each inner node of an octree can contain exactly 8 children.
   */
  std::array<std::unique_ptr<OctreeNodeInterface<Particle>>, 8> _children;
};
}  // namespace autopas
