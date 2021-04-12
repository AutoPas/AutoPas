/**
 * @file Octree.h
 *
 * @author Johannes Spies
 * @date 09.04.2021
 */

#pragma once

#include "autopas/cells/FullParticleCell.h"
#include "autopas/containers/CellBasedParticleContainer.h"
#include <cstdio>

namespace autopas {

using Position = std::array<double, 3>;

static inline Position operator+(Position a, Position b) {
  Position result = {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
  return result;
}

static inline Position operator*(double s, Position p) {
  Position result = {s*p[0], s*p[1], s*p[2]};
  return result;
}

// TODO(johannes): Documentation
template <class Particle>
class Octree : public CellBasedParticleContainer<FullParticleCell<Particle>> {
 public:
  using ParticleCell = FullParticleCell<Particle>;
  using ParticleType = typename ParticleCell::ParticleType;

  Octree(Position boxMin, Position boxMax, const double cutoff,
         const double skin) : CellBasedParticleContainer<ParticleCell>(boxMin, boxMax, cutoff, skin), _root(nullptr) {
    printf("Johannes' Octree()\n");
    initOctree(boxMin, boxMax);
  }

  [[nodiscard]] std::vector<ParticleType> updateContainer() override {
    printf("Johannes' Octree::updateContainer\n");
    auto result = std::vector<ParticleType>();
    return result;
  }

  void iteratePairwise(TraversalInterface *traversal) override {
    printf("Johannes' Octree::iteratePairwise\n");
  }

  /**
   * @copydoc ParticleContainerInterface::getContainerType()
   */
  [[nodiscard]] ContainerOption getContainerType() const override { return ContainerOption::octree; }

  /**
   * @copydoc ParticleContainerInterface::getParticleCellTypeEnum()
   */
  [[nodiscard]] CellType getParticleCellTypeEnum() override { return CellType::FullParticleCell; }

  /**
   * @copydoc ParticleContainerInterface::addParticleImpl()
   */
  void addParticleImpl(const ParticleType &p) override {
    //printf("Johannes' Octree::addParticleImpl\n");
    Node::insert(_root, p);
  }

  /**
   * @copydoc ParticleContainerInterface::addHaloParticleImpl()
   */
  void addHaloParticleImpl(const ParticleType &haloParticle) override {
    printf("Johannes' Octree::addHaloParticleImpl\n");
  }

  /**
   * @copydoc ParticleContainerInterface::updateHaloParticle()
   */
  bool updateHaloParticle(const ParticleType &haloParticle) override {
    printf("Johannes' Octree::updateHaloParticle\n");
    return true;
  }

  void deleteHaloParticles() override {
    printf("Johannes' Octree::deleteHaloParticles\n");
  }

  void rebuildNeighborLists(TraversalInterface *traversal) override {
    printf("Johannes' Octree::rebuildNeighborLists\n");
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, true> begin(IteratorBehavior behavior) override {
    printf("Johannes' Octree::begin<..., true>\n");
    return ParticleIteratorWrapper<ParticleType, true>();
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, false> begin(IteratorBehavior behavior) const override {
    printf("Johannes' Octree::begin<..., false>\n");
    return ParticleIteratorWrapper<ParticleType, false>();
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, true> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior) override {
    printf("Johannes' Octree::getRegionIterator<..., true>\n");
    return ParticleIteratorWrapper<ParticleType, true>();
  }

  [[nodiscard]] ParticleIteratorWrapper<ParticleType, false> getRegionIterator(
      const std::array<double, 3> &lowerCorner, const std::array<double, 3> &higherCorner,
      IteratorBehavior behavior) const override {
    printf("Johannes' Octree::getRegionIterator<..., false>\n");
    return ParticleIteratorWrapper<ParticleType, false>();
  }

  /**
   * @copydoc ParticleContainerInterface::getTraversalSelectorInfo()
   */
  [[nodiscard]] TraversalSelectorInfo getTraversalSelectorInfo() const override {
    // TODO(johannes): Figure out what these values should be
    std::array<unsigned long, 3> dims = {1, 1, 1};
    std::array<double, 3> cellLength = {1, 1, 1};
    return TraversalSelectorInfo(dims, 0.0, cellLength, 1);
  }

 private:
  struct Node {
    /** Node type that distinguishes inner nodes from leaf nodes.*/
    enum Type {
      INNER, LEAF,
    } type;

    /** Axis aligned bounding box */
    Position boxMin, boxMax;

    /** Since the node can either be a leaf or an inner node, the memory can be used depending on the current type. */
#if 0
    union {
      // Fields for type == Type::INNER
      struct {
        Node *children[8];
      };

      // Fields for type == Type::LEAF
      struct {
        //std::vector<const ParticleType &> particles;
        //std::vector<int> particles;
        //int particles[MAX_NUMBER_OF_PARTICLES_PER_NODE];
        std::array<const ParticleType &, MAX_NUMBER_OF_PARTICLES_PER_NODE> particles;
        int unsigned particlesUsed;
      };
    };
#else
    Node *children[8];
    std::vector<Particle> particles;
#endif

    explicit Node(Type type, Position nodeMinCorner, Position nodeMaxCorner)
        : type(type),
          boxMin(nodeMinCorner),
          boxMax(nodeMaxCorner) {}

    // TODO(johannes): This functionality may already exist
    /**
     * Check if a 3d point is inside the node's axis aligned bounding box. (Set by the boxMin and boxMax fields.)
     * @param node The possible enclosing node
     * @param point The node to test
     * @return true if the point is inside the node's bounding box and false otherwise
     */
    static bool isInside(struct Node *node, Position point) {
      bool result = true;
      // Check if the point is within the bounding box for each dimension.
      for(auto d = 0; d < 3; d++) {
        result &= (node->boxMin[d] <= point[d]) && (point[d] <= node->boxMax[d]);
      }
      return result;
    }

    /**
     * Check if the node's axis aligned bounding box overlaps with the given axis aligned bounding box.
     * @param otherMin The minimum coordinate of the other box
     * @param otherMax The maximum coordinate of the other box
     * @return true iff the overlapping volume is non-negative
     */
    bool overlapsBox(Position otherMin, Position otherMax) {
      bool result = true;
      for(auto d = 0; d < 3; ++d) {
        result &= (boxMin[d] <= otherMax[d]) && (boxMax[d] >= otherMin[d]);
      }
      return result;
    }

    static void insert(Node *parent, Particle p) {
      assert(parent != nullptr); // initOctree() should have been called already
      const int unsigned maxParticlesInLeaf = 16;

      if (isInside(parent, p.getR())) {
        if (parent->type == Node::Type::LEAF) {
          if (parent->particles.size() < maxParticlesInLeaf) {
            parent->particles.push_back(p);
            //parent->particles[parent->particlesUsed++] = p;
          } else {
            // Create a new subdivision
            auto parentMin = parent->boxMin;
            auto parentCenter = 0.5 * (parent->boxMin + parent->boxMax);
            auto parentMax = parent->boxMax;
            for (auto i = 0; i < 8; i++) {
              Position boxMin = {};
              Position boxMax = {};

              // Subdivide the bounding box of the parent.
              for (auto d = 0; d < 3; ++d) {
                auto mask = 1 << d;
                boxMin[d] = !(i & mask) ? parentMin[d] : parentCenter[d];
                boxMax[d] = !(i & mask) ? parentCenter[d] : parentMax[d];
              }

              auto *newChild = new Node(Node::Type::LEAF, boxMin, boxMax);

              parent->children[i] = newChild;
            }
            parent->type = Node::Type::INNER;

            for(auto oldParticle : parent->particles) {
              insert(parent, oldParticle);
            }
            parent->particles.clear();

            // Insert the particle
            insert(parent, p);
          }
        } else {
          // Decide in which child the particle should be inserted.
          for (auto i = 0; i < 8; i++) {
            auto *child = parent->children[i];
            if (isInside(child, p.getR())) {
              insert(child, p);
              break;  // The particle should only be inserted into one box to avoid duplicate force calculations
            }
          }
        }
      }
    }
  };

  Node *_root;

  void initOctree(Position boxMin, Position boxMax) {
    _root = new Node(Node::Type::LEAF, boxMin, boxMax);
  }
};

} // namespace autopas
