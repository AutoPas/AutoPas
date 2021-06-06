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

template <typename Particle>
class OctreeLeafNode : public OctreeNodeInterface<Particle>, public FullParticleCell<Particle> {
 public:
  OctreeLeafNode(std::array<double, 3> boxMin, std::array<double, 3> boxMax, OctreeNodeInterface<Particle> *parent)
      : OctreeNodeInterface<Particle>(boxMin, boxMax, parent),
        FullParticleCell<Particle>(utils::ArrayMath::sub(boxMax, boxMin)) {}

  /**
   * @copydoc OctreeNodeInterface::insert()
   */
  void insert(std::unique_ptr<OctreeNodeInterface<Particle>> &ref, Particle p) override {
    if (!this->isInside(p.getR())) {
      throw std::runtime_error("[OctreeLeafNode.h] Attempting to insert particle that is not inside this node");
    }

    // TODO(johannes): Make this constant tunable or move it to a better suited location
    const int unsigned maxParticlesInLeaf = 4;

    if (this->_particles.size() < maxParticlesInLeaf) {
      this->_particles.push_back(p);
    } else {
      std::unique_ptr<OctreeNodeInterface<Particle>> newInner =
          std::make_unique<OctreeInnerNode<Particle>>(this->getBoxMin(), this->getBoxMax(), this->_parent);
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

  OctreeNodeInterface<Particle> *SON(Octant O) override {
    throw std::runtime_error("Unable to get SON of leaf node");
  }

  void appendAllLeaves(std::vector<OctreeLeafNode<Particle> *> &leaves) override {
    leaves.push_back(this);
  }
};
}  // namespace autopas
