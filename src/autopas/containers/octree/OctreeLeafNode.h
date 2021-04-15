/**
 * @file OctreeLeafNode.h
 *
 * @author Johannes Spies
 * @date 15.04.2021
 */
#pragma once

#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/containers/octree/OctreeInnerNode.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {
    template<typename Particle>
    class OctreeLeafNode : public OctreeNodeInterface<Particle> {
    public:
        OctreeLeafNode(std::array<double, 3> boxMin, std::array<double, 3> boxMax)
                : OctreeNodeInterface<Particle>(boxMin, boxMax) {}

        /**
         * @copydoc OctreeNodeInterface::insert()
         */
        OctreeNodeInterface<Particle> *insert(Particle p) override {
            using namespace autopas::utils;

            // TODO(johannes): Make this constant tunable or move it to a better suited location
            const int unsigned maxParticlesInLeaf = 16;

            OctreeNodeInterface<Particle> *result;
            if (particles.size() < maxParticlesInLeaf) {
                particles.push_back(p);
                result = this;
            } else {
                // Create a new subdivision based on the lower, center and maximum coordinate
                auto parentMin = this->getBoxMin();
                auto parentCenter = ArrayMath::mulScalar(ArrayMath::add(this->getBoxMin(), this->getBoxMax()), 0.5);
                auto parentMax = this->getBoxMax();

                std::array<OctreeNodeInterface<Particle> *, 8> newChildren;
                for (auto i = 0; i < 8; i++) {
                    std::array<double, 3> newBoxMin = {};
                    std::array<double, 3> newBoxMax = {};

                    // Subdivide the bounding box of the parent.
                    for (auto d = 0; d < 3; ++d) {
                        auto mask = 1 << d;
                        newBoxMin[d] = !(i & mask) ? parentMin[d] : parentCenter[d];
                        newBoxMax[d] = !(i & mask) ? parentCenter[d] : parentMax[d];
                    }

                    newChildren[i] = new OctreeLeafNode<Particle>(newBoxMin, newBoxMax);
                }
                result = new OctreeInnerNode<Particle>(parentMin, parentMax, newChildren);

                // Put all particles from the leaf inside the new nodes
                for (auto oldParticle : particles) {
                    result = result->insert(oldParticle);
                }
                particles.clear();

                // Insert the new particle
                result = result->insert(p);
            }

            return result;
        }

        /**
         * @copydoc OctreeNodeInterface::appendAllParticles()
         */
        void appendAllParticles(std::vector<Particle> &ps) override {
            ps.insert(ps.end(), particles.begin(), particles.end());
        }

        /**
         * @copydoc OctreeNodeInterface::clearChildren()
         */
        OctreeNodeInterface<Particle> *clearChildren() override {
            return this;
        }

    private:
        std::vector<Particle> particles;
    };
} // namespace autopas
