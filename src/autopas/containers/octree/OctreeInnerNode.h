/**
 * @file OctreeInnerNode.h
 *
 * @author Johannes Spies
 * @date 15.04.2021
 */

#pragma once

#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/containers/octree/OctreeStaticNodeSelector.h"
#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/inBox.h"

namespace autopas {
/**
 * Inner nodes of the octree data structure. An inner node always points to eight children, which can either be leaves
 * or inner nodes as well.
 *
 * @tparam Particle_T
 */
template <class Particle_T>
class OctreeInnerNode : public OctreeNodeInterface<Particle_T> {
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
  OctreeInnerNode(const std::array<double, 3> &boxMin, const std::array<double, 3> &boxMax,
                  OctreeNodeInterface<Particle_T> *parent, int unsigned treeSplitThreshold, double interactionLength,
                  double cellSizeFactor)
      : OctreeNodeInterface<Particle_T>(boxMin, boxMax, parent, treeSplitThreshold, interactionLength, cellSizeFactor) {
    using namespace autopas::utils;
    using namespace autopas::utils::ArrayMath::literals;

    // The inner node is initialized with 8 leaves.
    const auto center = (boxMin + boxMax) * 0.5;
    for (auto i = 0; i < _children.size(); ++i) {
      // Subdivide the bounding box of the parent.
      std::array<double, 3> newBoxMin = {};
      std::array<double, 3> newBoxMax = {};
      for (auto d = 0; d < 3; ++d) {
        // `i`, `d` and `mask` are used to generate minimum and maximum coordinates for every octree leaf's bounding
        // box. `i` represents a 3-bit wide number, where each bit corresponds to an axis. The following table
        // visualizes this layout:
        //
        //          <-- msb   lsb -->
        // +------++-----+---+---+---+
        // | bit  || ... | 2 | 1 | 0 |
        // +------+------+---+---+---+
        // | axis ||       x | y | z |
        // +------++---------+---+---+
        //
        // When an axis bit is not set, a leaf is expected to range from the minimum coordinate to the center of the
        // enclosing box on the respective axis `d`. If the axis bit is set, the leaf's region should start at the
        // center coordinate and end at the maximum coordinate of the enclosing box.
        //
        // `mask` is used to extract the individual axis components from `i`: Using the shift operation, `mask` can
        // either become 0b100, 0b010 or 0b001 since `d` ranges from 0 to 2 inclusively. This covers all available axis
        // in 3 dimensions.
        const auto mask = 4 >> d;
        newBoxMin[d] = (not(i & mask)) ? boxMin[d] : center[d];
        newBoxMax[d] = (not(i & mask)) ? center[d] : boxMax[d];
      }

      // Assign new leaves as the children.
      _children[i] = std::make_unique<OctreeLeafNode<Particle_T>>(newBoxMin, newBoxMax, this, treeSplitThreshold,
                                                                  interactionLength, cellSizeFactor);
    }
  }

  /**
   * Copy all children from the other octree into this octree. (Create a new, copied subtree)
   * @param other The other octree (to copy from)
   */
  OctreeInnerNode(const OctreeInnerNode<Particle_T> &other)
      : OctreeNodeInterface<Particle_T>(other._boxMin, other._boxMax, other._parent, other._treeSplitThreshold,
                                        other._interactionLength, other._cellSizeFactor) {
    for (auto i = 0; i < other._children.size(); ++i) {
      auto *otherChild = other._children[i].get();
      if (otherChild->hasChildren()) {
        _children[i] =
            std::make_unique<OctreeInnerNode<Particle_T>>((OctreeInnerNode<Particle_T> &)*other._children[i]);
      } else {
        _children[i] = std::make_unique<OctreeLeafNode<Particle_T>>((OctreeLeafNode<Particle_T> &)*other._children[i]);
      }
    }
  }

  /**
   * @copydoc OctreeNodeInterface::insert()
   */
  std::unique_ptr<OctreeNodeInterface<Particle_T>> insert(const Particle_T &p) override {
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

  bool deleteParticle(Particle_T &particle) override {
    for (auto &child : _children) {
      if (child->isInside(particle.getR())) {
        return child->deleteParticle(particle);
      }
    }
    // clang-format off
    utils::ExceptionHandler::exception(
        "Particle not found in this node!"
        "\nBoxMin: " + utils::ArrayUtils::to_string(this->_boxMin) +
        "\nBoxMax: " + utils::ArrayUtils::to_string(this->_boxMax) +
        "\n" + particle.toString());
    // clang-format on
    return false;
  }

  /**
   * @copydoc OctreeNodeInterface::collectAllParticles()
   */
  void collectAllParticles(std::vector<Particle_T *> &ps) const override {
    // An inner node does not contain particles, traverse down to the children.
    for (auto &child : _children) {
      child->collectAllParticles(ps);
    }
  }

  /**
   * @copydoc OctreeNodeInterface::appendAllLeafBoxes()
   */
  void appendAllLeafBoxes(std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> &boxes) const override {
    for (auto &child : _children) {
      child->appendAllLeafBoxes(boxes);
    }
  }

  /**
   * @copydoc OctreeNodeInterface::clearChildren()
   */
  void clearChildren(std::unique_ptr<OctreeNodeInterface<Particle_T>> &ref) override {
    for (auto &child : _children) {
      child->clearChildren(child);
    }

    std::unique_ptr<OctreeLeafNode<Particle_T>> newLeaf = std::make_unique<OctreeLeafNode<Particle_T>>(
        this->getBoxMin(), this->getBoxMax(), this->_parent, this->_treeSplitThreshold, this->_interactionLength,
        this->_cellSizeFactor);
    ref = std::move(newLeaf);
  }

  /**
   * @copydoc OctreeNodeInterface::size()
   */
  size_t size() const override {
    unsigned int result = 0;
    for (const auto &child : _children) {
      result += child->size();
    }
    return result;
  }

  /**
   * @copydoc OctreeNodeInterface::getNumberOfParticles()
   */
  size_t getNumberOfParticles(IteratorBehavior behavior) const override {
    unsigned int result = 0;
    for (const auto &child : _children) {
      result += child->getNumberOfParticles(behavior);
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
  OctreeNodeInterface<Particle_T> *getChild(int index) override { return _children[index].get(); }

  std::vector<OctreeLeafNode<Particle_T> *> getLeavesFromDirections(
      const std::vector<octree::Vertex> &directions) override {
    std::vector<OctreeLeafNode<Particle_T> *> result;
    // Only take the children that are allowed (i.e. those which are in the given directions list)
    for (auto d : directions) {
      int childIndex = vertexToIndex(d);
      if (childIndex < 0) {
        throw std::runtime_error("[OctreeInnerNode::getLeavesFromDirections()] Calculated invalid child index");
      }

      auto child = getChild(childIndex);
      auto childLeaves = child->getLeavesFromDirections(directions);
      result.insert(result.end(), childLeaves.begin(), childLeaves.end());
    }
    return result;
  }

  OctreeNodeInterface<Particle_T> *SON(octree::Octant octant) override {
    // convert the Octant to a flat child index
    auto flat = vertexToIndex(octant);
    return _children[flat].get();
  }

  void appendAllLeaves(std::vector<OctreeLeafNode<Particle_T> *> &leaves) const override {
    for (auto &child : _children) {
      child->appendAllLeaves(leaves);
    }
  }

  std::set<OctreeLeafNode<Particle_T> *> getLeavesInRange(const std::array<double, 3> &min,
                                                          const std::array<double, 3> &max) override {
    std::set<OctreeLeafNode<Particle_T> *> result;
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

  /**
   * @copydoc OctreeNodeWrapper::forEach()
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda) {
    for (auto &child : _children) {
      withStaticNodeType(child, [&](auto nodePtr) { nodePtr->forEach(forEachLambda); });
    }
  }

  /**
   * @copydoc OctreeNodeWrapper::reduce()
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result) {
    for (auto &child : _children) {
      withStaticNodeType(child, [&](auto nodePtr) { nodePtr->reduce(reduceLambda, result); });
    }
  }

  /**
   * Apply the forEach lambda to each particle in the region.
   *
   * @tparam Lambda Function type
   * @param forEachLambda Function to apply
   * @param lowerCorner Lower corner of region
   * @param higherCorner Higher corner of region
   * @param behavior Parameter is only there to reuse functionality already implemented in FullParticleCell, should be
   * set to ownedOrHalo
   */
  template <typename Lambda>
  void forEach(Lambda forEachLambda, const std::array<double, 3> &lowerCorner,
               const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    for (auto &child : _children) {
      double vol = child->getEnclosedVolumeWith(lowerCorner, higherCorner);
      if (vol > 0.0)
        withStaticNodeType(child,
                           [&](auto nodePtr) { nodePtr->forEach(forEachLambda, lowerCorner, higherCorner, behavior); });
    }
  }

  /**
   * Apply the reduce lambda to each particle in the region.
   *
   * @tparam Lambda Function type
   * @tparam A Initial value type
   * @param reduceLambda Function to apply
   * @param result Initial value
   * @param lowerCorner Lower corner of region
   * @param higherCorner Higher corner of region
   * @param behavior Parameter is only there to reuse functionality already implemented in FullParticleCell, should be
   * set to ownedOrHalo
   */
  template <typename Lambda, typename A>
  void reduce(Lambda reduceLambda, A &result, const std::array<double, 3> &lowerCorner,
              const std::array<double, 3> &higherCorner, IteratorBehavior behavior) {
    for (auto &child : _children) {
      double vol = child->getEnclosedVolumeWith(lowerCorner, higherCorner);
      if (vol > 0.0)
        withStaticNodeType(
            child, [&](auto nodePtr) { nodePtr->reduce(reduceLambda, result, lowerCorner, higherCorner, behavior); });
    }
  }

 private:
  /**
   * Each inner node of an octree can contain exactly 8 children.
   */
  std::array<std::unique_ptr<OctreeNodeInterface<Particle_T>>, 8> _children;
};
}  // namespace autopas
