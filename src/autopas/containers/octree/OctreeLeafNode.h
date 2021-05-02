/**
 * @file OctreeLeafNode.h
 *
 * @author Johannes Spies
 * @date 15.04.2021
 */
#pragma once

#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/containers/octree/OctreeInnerNode.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/utils/ArrayMath.h"

namespace autopas {
    template<typename Particle>
    class Octree;

    template<typename Particle>
    class OctreeLeafNode : public OctreeNodeInterface<Particle>, public FullParticleCell<Particle> {
    public:
        OctreeLeafNode(std::array<double, 3> boxMin, std::array<double, 3> boxMax)
                : OctreeNodeInterface<Particle>(boxMin, boxMax),
                        FullParticleCell<Particle>(utils::ArrayMath::sub(boxMax, boxMin)) {}

        /**
         * @copydoc OctreeNodeInterface::insert()
         */
        void insert(std::unique_ptr<OctreeNodeInterface<Particle>> &ref, Particle p) override {
            if(!this->isInside(p.getR())) {
                throw std::runtime_error("Attempting to insert particle that is not inside this node");
            }

            // TODO(johannes): Make this constant tunable or move it to a better suited location
            const int unsigned maxParticlesInLeaf = 16;

            if(this->_particles.size() < maxParticlesInLeaf) {
                this->_particles.push_back(p);
            } else {
                std::unique_ptr<OctreeNodeInterface<Particle>> newInner =
                        std::make_unique<OctreeInnerNode<Particle>>(this->getBoxMin(), this->getBoxMax());
                newInner->insert(newInner, p);
                for(auto cachedParticle : this->_particles) {
                    newInner->insert(newInner, cachedParticle);
                }

                // Set the reference of the parent to this leaf to the new inner node.
                ref = std::move(newInner);
            }
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
        void clearChildren(std::unique_ptr<OctreeNodeInterface<Particle>> &ref) override {
            // Nothing to do for a leaf.
        }

        /**
         * @copydoc OctreeNodeInterface::getNumParticles()
         */
        unsigned int getNumParticles() override {
            return this->_particles.size();
        }
    };
} // namespace autopas
