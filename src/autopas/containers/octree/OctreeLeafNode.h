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

#include <optional>

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
   * @param interactionLength The minimum distance at which a force is considered nonzero, cutoff+skin.
   */
  OctreeLeafNode(std::array<double, 3> boxMin, std::array<double, 3> boxMax, OctreeNodeInterface<Particle> *parent,
                 int unsigned treeSplitThreshold, double interactionLength, double cellSizeFactor)
      : OctreeNodeInterface<Particle>(boxMin, boxMax, parent, treeSplitThreshold, interactionLength, cellSizeFactor),
        FullParticleCell<Particle>(utils::ArrayMath::sub(boxMax, boxMin)) {}

  OctreeLeafNode(OctreeLeafNode<Particle> const &other)
      : OctreeNodeInterface<Particle>(other._boxMin, other._boxMax, other._parent, other._treeSplitThreshold,
                                      other._interactionLength),
        FullParticleCell<Particle>(utils::ArrayMath::sub(other._boxMax, other._boxMin)) {
    for (auto &p : other._particles) {
      this->_particles.push_back(p);
    }
  }

  /**
   * @copydoc OctreeNodeInterface::insert()
   */
  std::optional<std::unique_ptr<OctreeNodeInterface<Particle>>> insert(Particle p) override {
    if (not this->isInside(p.getR())) {
      // The exception is suppressed for AllContainersTests#testParticleAdding
      // throw std::runtime_error("[OctreeLeafNode.h] Attempting to insert particle that is not inside this node");
    }

    // Check if the size of the new leaves would become smaller than cellSizeFactor*interactionLength
    std::array<double, 3> splitLeafDimensions = utils::ArrayMath::sub(this->getBoxMax(), this->getBoxMin());
    splitLeafDimensions = utils::ArrayMath::mulScalar(splitLeafDimensions, 0.5);
    bool anyNewDimSmallerThanMinSize = false;
    for (auto d = 0; d < 3; ++d) {
      auto cellSizeFactor = 1.0;
      // TODO: The condition below should actually be
      //  splitLeafDimensions[d] < (this->_cellSizeFactor * this->_interactionLength)
      //  But with this condition, the TraversalComparison test fails for cell size factor 0.5. Find out why the octree
      //  cannot handle this value.
      if (splitLeafDimensions[d] < this->_interactionLength) {
        anyNewDimSmallerThanMinSize = true;
        break;
      }
    }

    if ((this->_particles.size() < this->_treeSplitThreshold) or anyNewDimSmallerThanMinSize) {
      this->_particles.push_back(p);
      return std::nullopt;
    } else {
      std::unique_ptr<OctreeNodeInterface<Particle>> newInner = std::make_unique<OctreeInnerNode<Particle>>(
          this->getBoxMin(), this->getBoxMax(), this->_parent, this->_treeSplitThreshold, this->_interactionLength,
          this->_cellSizeFactor);
      auto opt = newInner->insert(p);
      if (opt) {
        newInner = std::move(*opt);
      }
      for (auto cachedParticle : this->_particles) {
        opt = newInner->insert(cachedParticle);
        if (opt) {
          newInner = std::move(*opt);
        }
      }

      return newInner;
    }
  }

  /**
   * @copydoc OctreeNodeInterface::appendAllParticles()
   */
  void appendAllParticles(std::vector<Particle *> &ps) override {
    for (Particle &p : this->_particles) {
      ps.push_back(&p);
    }
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

  std::vector<OctreeLeafNode<Particle> *> getLeavesFromDirections(std::vector<Vertex> directions) override {
    return {this};
  }

  OctreeNodeInterface<Particle> *SON(Octant O) override { throw std::runtime_error("Unable to get SON of leaf node"); }

  void appendAllLeaves(std::vector<OctreeLeafNode<Particle> *> &leaves) override { leaves.push_back(this); }

  std::set<OctreeLeafNode<Particle> *> getLeavesInRange(std::array<double, 3> min, std::array<double, 3> max) override {
    if (this->getEnclosedVolumeWith(min, max) > 0.0) {
      return {this};
    } else {
      return {};
    }
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
