/**
 * @file OctreeNodeInterface.h
 *
 * @author Johannes Spies
 * @date 15.04.2021
 */
#pragma once

#include <array>
#include <vector>

#include "autopas/utils/inBox.h"
#include "autopas/cells/FullParticleCell.h"

namespace autopas {
    template<typename Particle>
    class OctreeLeafNode;

    /**
     * The base class that provides the necessary function definitions that can be applied to an octree.
     *
     * @tparam Particle
     */
    template<class Particle>
    class OctreeNodeInterface {
    public:
        OctreeNodeInterface(std::array<double, 3> boxMin, std::array<double, 3> boxMax)
                : _boxMin(boxMin), _boxMax(boxMax) {}

        /**
         * Insert a particle into the octree.
         * @param p The particle to insert
         * @return The subtree below the current node that now contains the particle
         */
        virtual void insert(std::unique_ptr<OctreeNodeInterface<Particle>> &ref, Particle p) = 0;

        /**
         * Put all particles that are below this node into the vector.
         * @param ps A reference to the vector that should contain the particles after the operation
         */
        virtual void appendAllParticles(std::vector<Particle> &ps) = 0;

        /**
         * Put the min/max corner coordinates of every leaf into the vector.
         * @param boxes A reference to the vector that should contain pairs of the min/max corner coordinates
         */
         virtual void appendAllLeafBoxes(std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> &boxes) = 0;

        /**
         * Delete the entire tree below this node.
         * @return A leaf node.
         */
        virtual void clearChildren(std::unique_ptr<OctreeNodeInterface<Particle>> &ref) = 0;

        /**
         * @copydoc CellBasedParticleContainer::getNumParticles()
         */
        virtual unsigned int getNumParticles() = 0;

        /**
         * Check if a 3d point is inside the node's axis aligned bounding box. (Set by the boxMin and boxMax fields.)
         * @param node The possible enclosing node
         * @param point The node to test
         * @return true if the point is inside the node's bounding box and false otherwise
         */
        bool isInside(std::array<double, 3> point) {
            using namespace autopas::utils;
            return inBox(point, _boxMin, _boxMax);
        }

        /**
         * Check if the node's axis aligned bounding box overlaps with the given axis aligned bounding box.
         * @param otherMin The minimum coordinate of the other box
         * @param otherMax The maximum coordinate of the other box
         * @return true iff the overlapping volume is non-negative
         */
        bool overlapsBox(std::array<double, 3> otherMin, std::array<double, 3> otherMax) {
            bool result = true;
            for (auto d = 0; d < 3; ++d) {
                result &= (this->_boxMin[d] <= otherMax[d]) && (this->_boxMax[d] >= otherMin[d]);
            }
            return result;
        }

        void setBoxMin(std::array<double, 3> boxMin) {_boxMin = boxMin;}
        void setBoxMax(std::array<double, 3> boxMax) {_boxMax = boxMax;}
        std::array<double, 3> getBoxMin() {return _boxMin;}
        std::array<double, 3> getBoxMax() {return _boxMax;}

    protected:
        std::array<double, 3> _boxMin, _boxMax;
    };
} // namespace autpas
