/**
 * @file OctreeInnerNode.h
 *
 * @author Johannes Spies
 * @date 15.04.2021
 */
#pragma once

#include <variant>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/inBox.h"
#include "autopas/containers/octree/OctreeNodeInterface.h"

namespace autopas {
    template<class Particle>
    class OctreeInnerNode : public OctreeNodeInterface<Particle> {
    public:
        OctreeInnerNode(std::array<double, 3> boxMin, std::array<double, 3> boxMax)
                : OctreeNodeInterface<Particle>(boxMin, boxMax) {
            using namespace autopas::utils;

            // The inner node is initialized with 8 leaves.
            auto center = ArrayMath::mulScalar(ArrayMath::add(boxMin, boxMax), 0.5);
            for(auto i = 0; i < _children.size(); ++i) {
                // Subdivide the bounding box of the parent.
                std::array<double, 3> newBoxMin = {};
                std::array<double, 3> newBoxMax = {};
                for (auto d = 0; d < 3; ++d) {
                    auto mask = 1 << d;
                    newBoxMin[d] = !(i & mask) ? boxMin[d] : center[d];
                    newBoxMax[d] = !(i & mask) ? center[d] : boxMax[d];
                }

                // Assign new leaves as the children.
                _children[i] = std::make_unique<OctreeLeafNode<Particle>>(newBoxMin, newBoxMax);
            }
        }

        /**
         * @copydoc OctreeNodeInterface::insert()
         */
        void insert(std::unique_ptr<OctreeNodeInterface<Particle>> &ref, Particle p) override {
            if(!this->isInside(p.getR())) {
                throw std::runtime_error("Attempting to insert particle that is not inside this node");
            }

            // Find a child to insert the particle into.
            for(auto &child : _children) {
                if(child->isInside(p.getR())) {
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
        void clearChildren(std::unique_ptr<OctreeNodeInterface<Particle>> &ref) {
            for (auto &child : _children) {
                child->clearChildren(child);
            }

            std::unique_ptr<OctreeLeafNode<Particle>> newLeaf = std::make_unique<OctreeLeafNode<Particle>>(this->getBoxMin(), this->getBoxMax());
            ref = std::move(newLeaf);
        }

        /**
         * @copydoc OctreeNodeInterface::getNumParticles()
         */
        unsigned int getNumParticles() override {
            unsigned int result = 0;
            for(auto &child : _children) {
                result += child->getNumParticles();
            }
            return result;
        }

    private:
        /**
         * Each inner node of an octree can contain exactly 8 children.
         */
        std::array<std::unique_ptr<OctreeNodeInterface<Particle>>, 8> _children;
    };
} // namespace autopas
