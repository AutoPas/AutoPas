/**
 * @file OctreeTests.cpp
 * @author Johannes Spies
 * @date 15.04.2021
 */

#include "OctreeTests.h"
#include "autopas/containers/octree/Octree.h"
#include "autopas/particles/Particle.h"

TEST_F(OctreeTest, testDummy) {
    using namespace autopas;

    std::array<double, 3> min = {0, 0, 0}, max = {2, 2, 2};
    Octree<ParticleFP64> tree(min, max, 0.001f, 0.1f);
}