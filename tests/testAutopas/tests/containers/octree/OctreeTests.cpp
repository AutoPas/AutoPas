/**
 * @file OctreeTests.cpp
 * @author Johannes Spies
 * @date 15.04.2021
 */

#include "OctreeTests.h"

#include "autopas/containers/octree/Octree.h"
#include "autopas/containers/octree/OctreeDirection.h"
#include "autopas/particles/Particle.h"

TEST_F(OctreeTest, testDummy) {
  using namespace autopas;

  std::array<double, 3> min = {0, 0, 0}, max = {2, 2, 2};
  Octree<ParticleFP64> tree(min, max, 0.001f, 0.1f);
}

/**
 * This test checks (pairwise) whether the boxes split on one layer are correctly aligned with themselves.
 */
TEST_F(OctreeTest, testDebugIndexing) {
  using namespace autopas;

  // insert more than 4 dummy particles to trigger splitting of the root leaf.

  std::array<double, 3> min = {0, 0, 0}, max = {1, 1, 1};
  std::unique_ptr<OctreeNodeInterface<ParticleFP64>> root = std::make_unique<OctreeLeafNode<ParticleFP64>>(min, max, nullptr);
  // Add some dummy particles that split the nodes
  int dummyParticleCount = 8;
  for (int i = 0; i < 8; ++i) {
    root->insert(root, ParticleFP64({0.01f, (double)i / (double)dummyParticleCount, 0.01f}, {0, 0, 0}, 0));
  }

  int axisPairs[3][4][2] = {
      {
          // x axis cases
          {0b000, 0b001},
          {0b010, 0b011},
          {0b100, 0b101},
          {0b110, 0b111},
      },
      {
          // y axis cases
          {0b000, 0b010},
          {0b001, 0b011},
          {0b100, 0b110},
          {0b101, 0b111},
      },
      {
          // z axis cases
          {0b000, 0b100},
          {0b001, 0b101},
          {0b010, 0b110},
          {0b011, 0b111},
      },
  };

  ASSERT_TRUE(root->hasChildren());
  for (int axis = 0; axis < 3; ++axis) {
    for (int option = 0; option < 4; ++option) {
      int box1Index = axisPairs[axis][option][0];
      int box2Index = axisPairs[axis][option][1];

      OctreeNodeInterface<ParticleFP64> *box1 = root->getChild(box1Index);
      OctreeNodeInterface<ParticleFP64> *box2 = root->getChild(box2Index);
      ASSERT_NE(box1, box2);

      auto box1Min = box1->getBoxMin();
      auto box1Max = box1->getBoxMax();
      auto box2Min = box2->getBoxMin();
      auto box2Max = box2->getBoxMax();
      ASSERT_GE(box1Max[axis], box2Min[axis]);
      ASSERT_LT(box1Min[axis], box2Min[axis]);
    }
  }
}

/**
 * This is a giant test case that simulates an entire octree run in one route.
 */
TEST_F(OctreeTest, testParentFinding) {
  using namespace autopas;

  using Position = std::array<double, 3>;

  std::array<double, 3> min = {0, 0, 0}, max = {1, 1, 1};
  std::unique_ptr<OctreeNodeInterface<ParticleFP64>> root = std::make_unique<OctreeLeafNode<ParticleFP64>>(min, max, nullptr);

  // Create some particle positions that are in a specific cell
  auto maxParticlesPerCell = 4;
  auto dummyParticleCount = 4*maxParticlesPerCell;
  std::vector<Position> positions;
  for (int i = 0; i < dummyParticleCount; ++i) {
    positions.push_back({(double)i/(double)dummyParticleCount, 0.1, 0.1});
  }

  // Insert the particles
  for (int i = 0; i < dummyParticleCount; ++i) {
    root->insert(root, ParticleFP64(positions[i], {0, 0, 0}, 0));
  }

  // Get the split cells
  ASSERT_TRUE(root->hasChildren());
  OctreeNodeInterface<ParticleFP64> *child0 = root->getChild(0b000);
  ASSERT_TRUE(child0->hasChildren());
  OctreeNodeInterface<ParticleFP64> *child01 = child0->getChild(0b001);
  ASSERT_FALSE(child01->hasChildren());
  ASSERT_EQ(child01->getNumParticles(), 4);

  // Check if we get the node we wanted from getGreaterParentAlongAxis
  if(auto greaterParentChild01Optional = child01->getGreaterParentAlongAxis(0, 1, child01)) {
    auto greaterParentChild01 = *greaterParentChild01Optional;
    ASSERT_EQ(greaterParentChild01, root.get());
    ASSERT_TRUE(greaterParentChild01->hasChildren());

    auto phaseTwoStartNode = greaterParentChild01->getChild(0b001);
    auto touching = phaseTwoStartNode->findTouchingLeaves(0, -1, child01);

    // The elements returned from the find touching routine should all be leaves.
    for(auto *elem : touching) {
      ASSERT_FALSE(elem->hasChildren());
    }

    ASSERT_EQ(touching.size(), 1);
    auto methodTouching = child01->getNeighborsAlongAxis(0, 1);
    //ASSERT_EQ(touching[0], methodTouching[0]);
  } else {
    FAIL();
  }
}

TEST_F(OctreeTest, testDirectionIndexing) {
  using namespace autopas;

  ASSERT_EQ(NEG_X, getDirection(0, false));
  ASSERT_EQ(POS_X, getDirection(0, true));
  ASSERT_EQ(NEG_Y, getDirection(1, false));
  ASSERT_EQ(POS_Y, getDirection(1, true));
  ASSERT_EQ(NEG_Z, getDirection(2, false));
  ASSERT_EQ(POS_Z, getDirection(2, true));
}