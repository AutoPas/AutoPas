/**
 * @file OctreeLeafNode.h
 *
 * @author Johannes Spies
 * @date 15.04.2021
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/octree/OctreeInnerNode.h"
#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {
template <typename Particle>
class Octree;

/**
 * An octree leaf node. This class utilizes the FullParticleCell to store the actual particles.
 *
 * @tparam Particle
 */
template <typename Particle>
class OctreeLeafNode : public OctreeNodeInterface<Particle>, public FullParticleCell<Particle> {
 public:
  /**
   * Create an empty octree leaf node
   * @param boxMin The min coordinate of the octree box
   * @param boxMax The max coordinate of the octree box
   * @param parent A pointer to the parent node. Should be nullptr for root nodes.
   * @param treeSplitThreshold Maximum number of particles inside a leaf before it tries to split itself
   */
  OctreeLeafNode(std::array<double, 3> boxMin, std::array<double, 3> boxMax, OctreeNodeInterface<Particle> *parent,
                 int unsigned treeSplitThreshold)
      : OctreeNodeInterface<Particle>(boxMin, boxMax, parent, treeSplitThreshold),
        FullParticleCell<Particle>(utils::ArrayMath::sub(boxMax, boxMin)) {}

  /**
   * @copydoc OctreeNodeInterface::insert()
   */
  void insert(std::unique_ptr<OctreeNodeInterface<Particle>> &ref, Particle p) override {
    if (!this->isInside(p.getR())) {
      throw std::runtime_error("[OctreeLeafNode.h] Attempting to insert particle that is not inside this node");
    }

    // TODO(johannes): Check if the size of the new leaves is >= cellSizeFactor*interactionLength
    if (this->_particles.size() < this->_treeSplitThreshold) {
      this->_particles.push_back(p);
    } else {
      std::unique_ptr<OctreeNodeInterface<Particle>> newInner = std::make_unique<OctreeInnerNode<Particle>>(
          this->getBoxMin(), this->getBoxMax(), this->_parent, this->_treeSplitThreshold);
      newInner->insert(newInner, p);
      for (auto cachedParticle : this->_particles) {
        newInner->insert(newInner, cachedParticle);
      }

      // Set the reference of the parent to this leaf to the new inner node.
      ref = std::move(newInner);
    }

    // TODO: Take a look at this later
    // return std::make_unique<OctreeNodeInterface<Particle>>(*this);
  }

  /**
   * @copydoc OctreeNodeInterface::appendAllParticles()
   */
  void appendAllParticles(std::vector<Particle> &ps) override {
    ps.insert(ps.end(), this->_particles.begin(), this->_particles.end());
  }

  /**
   * @copydoc OctreeNodeInterface::appendAllLeafBoxes()
   */
  void appendAllLeafBoxes(std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> &boxes) override {
    auto minMax = std::make_pair(this->getBoxMin(), this->getBoxMax());
    boxes.push_back(minMax);
  }

  /**
   * @copydoc OctreeNodeInterface::clearChildren()
   */
  void clearChildren(std::unique_ptr<OctreeNodeInterface<Particle>> &ref) override { this->_particles.clear(); }

  /**
   * @copydoc OctreeNodeInterface::getNumParticles()
   */
  unsigned int getNumParticles() override { return this->_particles.size(); }

  /**
   * @copydoc OctreeNodeInterface::hasChildren()
   */
  bool hasChildren() override { return false; }

  /**
   * @copydoc OctreeNodeInterface::getChild()
   */
  OctreeNodeInterface<Particle> *getChild(int index) override {
    throw std::runtime_error("[OctreeLeafNode] Unable to return child by index in leaf");
  }

#if 0
  std::optional<OctreeNodeInterface<Particle> *> getGreaterParentAlongAxis(
      int axis, int dir, OctreeNodeInterface<Particle> *embedded) override {
    std::optional<OctreeNodeInterface<Particle> *> result = std::nullopt;
    if (this->hasParent()) {
      result = this->_parent->getGreaterParentAlongAxis(axis, dir, embedded);
    }
    return result;
  }

  std::vector<OctreeNodeInterface<Particle> *> findTouchingLeaves(int axis, int dir,
                                                                  OctreeNodeInterface<Particle> *embedded) override {
    int otherAxis1 = (axis + 1) % 3;
    int otherAxis2 = (axis + 2) % 3;
    if (!OctreeNodeInterface<Particle>::volumeExistsOnAxis(otherAxis1, this->getBoxMin(), this->getBoxMax(),
                                                           embedded->getBoxMin(), embedded->getBoxMax()) ||
        !OctreeNodeInterface<Particle>::volumeExistsOnAxis(otherAxis2, this->getBoxMin(), this->getBoxMax(),
                                                           embedded->getBoxMin(), embedded->getBoxMax())) {
      throw std::runtime_error("[OctreeLeafNode] Leaf does not overlap with requested box on necessary axis.");
    }

    std::vector<OctreeNodeInterface<Particle> *> result;
    result.push_back(this);
    return result;
  }
#endif

  std::vector<OctreeLeafNode<Particle> *> getLeavesFromDirections(std::vector<Vertex> directions) override {
    return {this};
  }

  OctreeNodeInterface<Particle> *SON(Octant O) override { throw std::runtime_error("Unable to get SON of leaf node"); }

  void appendAllLeaves(std::vector<OctreeLeafNode<Particle> *> &leaves) override { leaves.push_back(this); }

  std::set<OctreeNodeInterface<Particle> *> getLeavesInRange(std::array<double, 3> min,
                                                             std::array<double, 3> max) override {
    // TODO(johannes): Check what counts as zero
    if (this->getEnclosedVolumeWith(min, max) < 0.0001f) {
      throw std::runtime_error("[OctreeLeafNode.h] The given region does not overlap with this leaf");
    }
    return {this};
  }

  /**
   * Clear the list that contains all neighbor nodes that have already been processed.
   */
  void clearAlreadyProcessedList() { _alreadyProcessed.clear(); }

  /**
   * Check if a node has already been processed.
   * @param other A pointer to another node
   * @return True if the other node is in the already processed list
   */
  bool alreadyProcessed(OctreeLeafNode<Particle> *other) {
    return std::find(_alreadyProcessed.begin(), _alreadyProcessed.end(), other) != _alreadyProcessed.end();
  }

  /**
   * Put a neighbor in the already processed list.
   * @param other The node to put in
   */
  void markAlreadyProcessed(OctreeLeafNode<Particle> *other) { _alreadyProcessed.push_back(other); }

 private:
  /**
   * The list that contains the neighbors that have already been processed in one traversal run.
   */
  std::vector<OctreeLeafNode<Particle> *> _alreadyProcessed;
};
}  // namespace autopas
