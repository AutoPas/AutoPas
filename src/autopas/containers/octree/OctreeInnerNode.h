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
   * @param treeSplitThreshold Maximum number of particles inside a leaf before it tries to split itself
   * @param interactionLength The minimum distance at which a force is considered nonzero, cutoff+skin.
   * @param cellSizeFactor The cell size factor
   */
  OctreeInnerNode(std::array<double, 3> boxMin, std::array<double, 3> boxMax, OctreeNodeInterface<Particle> *parent,
                  int unsigned treeSplitThreshold, double interactionLength, double cellSizeFactor)
      : OctreeNodeInterface<Particle>(boxMin, boxMax, parent, treeSplitThreshold, interactionLength, cellSizeFactor) {
    using namespace autopas::utils;

    // The inner node is initialized with 8 leaves.
    auto center = ArrayMath::mulScalar(ArrayMath::add(boxMin, boxMax), 0.5);
    for (auto i = 0; i < _children.size(); ++i) {
      // Subdivide the bounding box of the parent.
      std::array<double, 3> newBoxMin = {};
      std::array<double, 3> newBoxMax = {};
      for (auto d = 0; d < 3; ++d) {
        auto mask = 4 >> d;
        newBoxMin[d] = !(i & mask) ? boxMin[d] : center[d];
        newBoxMax[d] = !(i & mask) ? center[d] : boxMax[d];
      }

      // Assign new leaves as the children.
      _children[i] = std::make_unique<OctreeLeafNode<Particle>>(newBoxMin, newBoxMax, this, treeSplitThreshold,
                                                                interactionLength, cellSizeFactor);
    }
  }

  /**
   * Copy all children from the other octree into this octree. (Create a new, copied subtree)
   * @param other The other octree (to copy from)
   */
  OctreeInnerNode(const OctreeInnerNode<Particle> &other)
      : OctreeNodeInterface<Particle>(other._boxMin, other._boxMax, other._parent, other._treeSplitThreshold,
                                      other._interactionLength, other._cellSizeFactor) {
    for (auto i = 0; i < other._children.size(); ++i) {
      auto *otherChild = other._children[i].get();
      if (otherChild->hasChildren()) {
        _children[i] = std::make_unique<OctreeInnerNode<Particle>>((OctreeInnerNode<Particle> &)*other._children[i]);
      } else {
        _children[i] = std::make_unique<OctreeLeafNode<Particle>>((OctreeLeafNode<Particle> &)*other._children[i]);
      }
    }
  }

  /**
   * @copydoc OctreeNodeInterface::insert()
   */
  std::unique_ptr<OctreeNodeInterface<Particle>> insert(Particle p) override {
    if (!this->isInside(p.getR())) {
      // The exception is suppressed for AllContainersTests#testParticleAdding
      // throw std::runtime_error("[OctreeInnerNode.h] Attempting to insert particle that is not inside this node");
    }

    // Find a child to insert the particle into.
    for (auto &child : _children) {
      if (child->isInside(p.getR())) {
        auto ret = child->insert(p);
        if (ret) child = std::move(ret);
        break;
      }
    }

    return nullptr;
  }

  /**
   * @copydoc OctreeNodeInterface::appendAllParticles()
   */
  void appendAllParticles(std::vector<Particle *> &ps) override {
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

    std::unique_ptr<OctreeLeafNode<Particle>> newLeaf = std::make_unique<OctreeLeafNode<Particle>>(
        this->getBoxMin(), this->getBoxMax(), this->_parent, this->_treeSplitThreshold, this->_interactionLength,
        this->_cellSizeFactor);
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

  std::vector<OctreeLeafNode<Particle> *> getLeavesFromDirections(std::vector<Vertex> directions) override {
    std::vector<OctreeLeafNode<Particle> *> result;
    // Only take the children that are allowed (i.e. those which are in the given directions list)
    for (auto d : directions) {
      int childIndex = vertexToIndex(d);
      if (childIndex == -1) {
        throw std::runtime_error("[OctreeInnerNode.h] Calculated invalid child index");
      }

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

  std::set<OctreeLeafNode<Particle> *> getLeavesInRange(std::array<double, 3> min, std::array<double, 3> max) override {
    std::set<OctreeLeafNode<Particle> *> result;
    for (auto &child : _children) {
      double vol = child->getEnclosedVolumeWith(min, max);
      // Prevent iteration of the subtree if it is unnecessary
      if (vol > 0.0) {
        auto leaves = child->getLeavesInRange(min, max);
        result.insert(leaves.begin(), leaves.end());
      }
    }
    return result;
  }

 private:
  /**
   * Each inner node of an octree can contain exactly 8 children.
   */
  std::array<std::unique_ptr<OctreeNodeInterface<Particle>>, 8> _children;
};
}  // namespace autopas
