/**
 * @file OctreeTests.cpp
 * @author Johannes Spies
 * @date 15.04.2021
 */

#include "OctreeTests.h"

#include <cstdio>

#include "autopas/containers/octree/Octree.h"
#include "autopas/containers/octree/OctreeDirection.h"
#include "autopas/containers/octree/OctreeNodeInterface.h"
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
  std::unique_ptr<OctreeNodeInterface<ParticleFP64>> root =
      std::make_unique<OctreeLeafNode<ParticleFP64>>(min, max, nullptr, 4);
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

#if 0
/**
 * This is a giant test case that simulates an entire octree run in one route.
 */
TEST_F(OctreeTest, testParentFinding) {
  using namespace autopas;

  using Position = std::array<double, 3>;

  std::array<double, 3> min = {0, 0, 0}, max = {1, 1, 1};
  std::unique_ptr<OctreeNodeInterface<ParticleFP64>> root =
      std::make_unique<OctreeLeafNode<ParticleFP64>>(min, max, nullptr);

  // Create some particle positions that are in a specific cell
  auto maxParticlesPerCell = 4;
  auto dummyParticleCount = 4 * maxParticlesPerCell;
  std::vector<Position> positions;
  for (int i = 0; i < dummyParticleCount; ++i) {
    positions.push_back({(double)i / (double)dummyParticleCount, 0.1, 0.1});
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
  if (auto greaterParentChild01Optional = child01->getGreaterParentAlongAxis(0, 1, child01)) {
    auto greaterParentChild01 = *greaterParentChild01Optional;
    ASSERT_EQ(greaterParentChild01, root.get());
    ASSERT_TRUE(greaterParentChild01->hasChildren());

    auto phaseTwoStartNode = greaterParentChild01->getChild(0b001);
    auto touching = phaseTwoStartNode->findTouchingLeaves(0, -1, child01);

    // The elements returned from the find touching routine should all be leaves.
    for (auto *elem : touching) {
      ASSERT_FALSE(elem->hasChildren());
    }

    ASSERT_EQ(touching.size(), 1);
    auto methodTouching = child01->getNeighborsAlongAxis(0, 1);
    // ASSERT_EQ(touching[0], methodTouching[0]);
  } else {
    FAIL();
  }
}
#endif

TEST_F(OctreeTest, testDirectionIndexing) {
  using namespace autopas;

  ASSERT_EQ(NEG_X, getDirection(0, false));
  ASSERT_EQ(POS_X, getDirection(0, true));
  ASSERT_EQ(NEG_Y, getDirection(1, false));
  ASSERT_EQ(POS_Y, getDirection(1, true));
  ASSERT_EQ(NEG_Z, getDirection(2, false));
  ASSERT_EQ(POS_Z, getDirection(2, true));
}

static bool isOdd(int n) { return (n % 2) == 1; }

/**
 * Evaluate the plane equation
 * <n, x-a> = 0
 * for n=planeNormal, a=pointOnPlane, x=test and <.,.> denotes the default inner product.
 *
 * @param planeNormal The normal vector of the plane
 * @param pointOnPlane An arbitrary point on the plane
 * @param test A point to put in the plane vector
 * @return A value close to zero if test lays in the plane
 */
static double planeEquation(std::array<double, 3> planeNormal, std::array<double, 3> pointOnPlane,
                            std::array<double, 3> test) {
  double inner = 0;
  for (int d = 0; d < 3; ++d) {
    inner += planeNormal[d] * (test[d] - pointOnPlane[d]);
  }
  return inner;
}

/**
 * Get the axis that the face is perpendicular to.
 *
 * @param f The face
 * @return An axis index: 0, 1 or 2
 */
static int getAxis(autopas::Face f) { return (f - 1) >> 1; }

/**
 * This testcase checks if the indexing of octree node children works as intended.
 */
TEST_F(OctreeTest, testChildIndexing) {
  using namespace autopas;

  // Create an inner node that is split once.
  std::array<double, 3> min = {0, 0, 0}, max = {1, 1, 1};
  OctreeInnerNode<ParticleFP64> inner(min, max, nullptr, 16);

  // Get the center of the node
  std::array<double, 3> center = utils::ArrayMath::mulScalar(utils::ArrayMath::add(min, max), 0.5);

  // Check if the child indexing works
  using Node = OctreeNodeInterface<ParticleFP64>;

  // LDB
  int ldbIndex = vertexToIndex(LDB);
  Node *ldb = inner.getChild(ldbIndex);
  ASSERT_EQ(ldb->getBoxMin(), min);
  ASSERT_EQ(ldb->getBoxMax(), center);

  // LDF
  int ldfIndex = vertexToIndex(LDF);
  Node *ldf = inner.getChild(ldfIndex);
  ASSERT_EQ(ldf->getBoxMin()[0], min[0]);
  ASSERT_EQ(ldf->getBoxMin()[1], min[1]);
  ASSERT_EQ(ldf->getBoxMin()[2], center[2]);
  ASSERT_EQ(ldf->getBoxMax()[0], center[0]);
  ASSERT_EQ(ldf->getBoxMax()[1], center[1]);
  ASSERT_EQ(ldf->getBoxMax()[2], max[2]);

  // LUB
  int lubIndex = vertexToIndex(LUB);
  Node *lub = inner.getChild(lubIndex);
  ASSERT_EQ(lub->getBoxMin()[0], min[0]);
  ASSERT_EQ(lub->getBoxMin()[1], center[1]);
  ASSERT_EQ(lub->getBoxMin()[2], min[2]);
  ASSERT_EQ(lub->getBoxMax()[0], center[0]);
  ASSERT_EQ(lub->getBoxMax()[1], max[1]);
  ASSERT_EQ(lub->getBoxMax()[2], center[2]);

  // LUF
  int lufIndex = vertexToIndex(LUF);
  Node *luf = inner.getChild(lufIndex);
  ASSERT_EQ(luf->getBoxMin()[0], min[0]);
  ASSERT_EQ(luf->getBoxMin()[1], center[1]);
  ASSERT_EQ(luf->getBoxMin()[2], center[2]);
  ASSERT_EQ(luf->getBoxMax()[0], center[0]);
  ASSERT_EQ(luf->getBoxMax()[1], max[1]);
  ASSERT_EQ(luf->getBoxMax()[2], max[2]);

  // RDB
  int rdbIndex = vertexToIndex(RDB);
  Node *rdb = inner.getChild(rdbIndex);
  ASSERT_EQ(rdb->getBoxMin()[0], center[0]);
  ASSERT_EQ(rdb->getBoxMin()[1], min[1]);
  ASSERT_EQ(rdb->getBoxMin()[2], min[2]);
  ASSERT_EQ(rdb->getBoxMax()[0], max[0]);
  ASSERT_EQ(rdb->getBoxMax()[1], center[1]);
  ASSERT_EQ(rdb->getBoxMax()[2], center[2]);

  // RDF
  int rdfIndex = vertexToIndex(RDF);
  Node *rdf = inner.getChild(rdfIndex);
  ASSERT_EQ(rdf->getBoxMin()[0], center[0]);
  ASSERT_EQ(rdf->getBoxMin()[1], min[1]);
  ASSERT_EQ(rdf->getBoxMin()[2], center[2]);
  ASSERT_EQ(rdf->getBoxMax()[0], max[0]);
  ASSERT_EQ(rdf->getBoxMax()[1], center[1]);
  ASSERT_EQ(rdf->getBoxMax()[2], max[2]);

  // RUB
  int rubIndex = vertexToIndex(RUB);
  Node *rub = inner.getChild(rubIndex);
  ASSERT_EQ(rub->getBoxMin()[0], center[0]);
  ASSERT_EQ(rub->getBoxMin()[1], center[1]);
  ASSERT_EQ(rub->getBoxMin()[2], min[2]);
  ASSERT_EQ(rub->getBoxMax()[0], max[0]);
  ASSERT_EQ(rub->getBoxMax()[1], max[1]);
  ASSERT_EQ(rub->getBoxMax()[2], center[2]);

  // RUF
  int rufIndex = vertexToIndex(RUF);
  Node *ruf = inner.getChild(rufIndex);
  ASSERT_EQ(ruf->getBoxMin(), center);
  ASSERT_EQ(ruf->getBoxMax(), max);
}

template <typename Particle>
static void verifyFaceNeighbor(autopas::Face face, autopas::OctreeNodeInterface<Particle> *node,
                               autopas::OctreeNodeInterface<Particle> *neighbor) {
  using namespace autopas;

  int touchAxis = getAxis(face);
  int touchAxisDir = isOdd(face) ? -1 : 1;
  int oppositeDir = -1 * touchAxisDir;

  auto nodePointInclInFacePlane = (touchAxisDir == -1) ? node->getBoxMin() : node->getBoxMax();
  auto neighborPointInclInFacePlane = (oppositeDir == -1) ? neighbor->getBoxMin() : neighbor->getBoxMax();

  std::array<double, 3> normal = {0, 0, 0};
  normal[touchAxis] = 1;

  // Check if the face touches the other cube
  ASSERT_DOUBLE_EQ(planeEquation(normal, nodePointInclInFacePlane, neighborPointInclInFacePlane), 0.0);
  ASSERT_DOUBLE_EQ(planeEquation(normal, neighborPointInclInFacePlane, nodePointInclInFacePlane), 0.0);

  int overlapAxis1 = (touchAxis + 1) % 3;
  int overlapAxis2 = (overlapAxis1 + 1) % 3;

  // Check if the neighbor has volume with the leaf in the direction of the face
  ASSERT_TRUE(node->enclosesVolumeWithOtherOnAxis(overlapAxis1, neighbor));
  ASSERT_TRUE(node->enclosesVolumeWithOtherOnAxis(overlapAxis2, neighbor));
}

/**
 * Create a 3D point in a given box
 * @param min The minimum coordinate of the box
 * @param max The maximum coordinate of the box
 * @return A poorly random-sampled point in the box
 */
static std::array<double, 3> random3D(std::array<double, 3> min, std::array<double, 3> max) {
  using namespace autopas;

  // Create a random point in the [0,0,0] to [1,1,1] cube
  std::array<double, 3> randomPosition = {(double)rand() / (double)RAND_MAX, (double)rand() / (double)RAND_MAX,
                                          (double)rand() / (double)RAND_MAX};

  // Map in the given space
  std::array<double, 3> dim = utils::ArrayMath::sub(max, min);
  randomPosition = utils::ArrayMath::mul(randomPosition, dim);  // Map [0,1] to [0,width on axis]
  randomPosition = utils::ArrayMath::add(randomPosition, min);

  return randomPosition;
}

/**
 * Create an octree filled with particles whose positions are created by a very bad random number generator
 * @param rootRef A reference in which the newly generated octree will be saved
 * @param seed A seed value for the very bad RNG
 * @param min The minimum coordinate of the octree's bounding box
 * @param max The maximum coordinate of the octree's bounding box
 * @param randomParticleCount How many particles should be spawned
 */
static void createRandomOctree(std::unique_ptr<autopas::OctreeNodeInterface<autopas::ParticleFP64>> &rootRef, int seed,
                               std::array<double, 3> min, std::array<double, 3> max, int randomParticleCount) {
  using namespace autopas;
  std::unique_ptr<autopas::OctreeNodeInterface<autopas::ParticleFP64>> root =
      std::make_unique<OctreeLeafNode<ParticleFP64>>(min, max, nullptr, 16);
  srand(seed);
  for (int particleIndex = 0; particleIndex < randomParticleCount; ++particleIndex) {
    auto randomPosition = random3D(min, max);
    root->insert(root, ParticleFP64(randomPosition, {0, 0, 0}, 0));
  }
  rootRef = std::move(root);
}

/**
 * This testcase
 * 1. instantiates an octree data structure (not the octree container itself),
 * 2. puts a fixed number of randomly distributed particles in it,
 * 3. obtains all leaves in the octree and
 * 4. checks for each leaf whether the neighbor finding routines return sane nodes for each available direction.
 * TODO: Write when a neighbor is considered sane
 */
TEST_F(OctreeTest, testNeighborLocator) {
  using namespace autopas;

  // Create an octree with a random particle configuration.
  std::unique_ptr<OctreeNodeInterface<ParticleFP64>> root;
  createRandomOctree(root, 1234, {0, 0, 0}, {1, 1, 1}, 1000);

  // Find all leaves
  std::vector<OctreeLeafNode<ParticleFP64> *> leaves;
  root->appendAllLeaves(leaves);
  // Log the leaves
  FILE *out = OctreeLogger::leavesToJSON(fopen("leaves.json", "w"), leaves);
  if (out) {
    fclose(out);
  }

  // Check the properties of each found neighbor for each leaf.
  int leafIndex = 0;
  for (auto leaf : leaves) {
    // Check for each available face if the returned neighbor is valid.
    for (Face *face = getFaces(); *face != O; ++face) {
      auto neighbor = leaf->GTEQ_FACE_NEIGHBOR(*face);
      if (neighbor != nullptr) {
        verifyFaceNeighbor(*face, leaf, neighbor);

        auto neighborLeaves = neighbor->getNeighborLeaves(*face);

        for (auto neighborLeaf : neighborLeaves) {
          ASSERT_NE(neighborLeaf, nullptr);
          verifyFaceNeighbor(*face, leaf, neighborLeaf);
        }
      } else {
        // TODO(johannes): The only case in which it is allowed for the GTEQ_FACE_NEIGHBOR method to return nullptr is
        //  when the requested face is also on a face of the enclosing min/max cube. This can be checked in this branch.
      }
    }

    // Check for each available edge if the returned neighbor is valid.
    for (Edge *edge = getEdges(); *edge != OO; ++edge) {
      auto neighbor = leaf->GTEQ_EDGE_NEIGHBOR(*edge);
      if (neighbor != nullptr) {
      } else {
        // TODO(johannes): Check if the leaf is touching the border in this direction
      }
    }

    // Check for each available vertex if the returned neighbor is valid.
    for (Vertex *vertex = VERTICES(); *vertex != OOO; ++vertex) {
      auto neighbor = leaf->GTEQ_VERTEX_NEIGHBOR(*vertex);
      if (neighbor != nullptr) {
        // Get all corners of the leaf and the neighbor

        // See if there matches only one
      } else {
        // TODO(johannes): Check if the leaf is touching the border in this direction
      }
    }

    ++leafIndex;
  }
}

/**
 * Check if the getEnclosedVolumeWith function works as expected
 */
TEST_F(OctreeTest, testOverlapVolume) {
  using namespace autopas;
  ASSERT_DOUBLE_EQ(OctreeNodeInterface<ParticleFP64>::getEnclosedVolumeWith({0, 0, 0}, {1, 1, 1}, {1, 0, 0}, {2, 1, 1}),
                   0);
  ASSERT_DOUBLE_EQ(OctreeNodeInterface<ParticleFP64>::getEnclosedVolumeWith({0, 0, 0}, {1, 1, 1}, {0, 0, 0}, {1, 1, 1}),
                   1);
  ASSERT_DOUBLE_EQ(
      OctreeNodeInterface<ParticleFP64>::getEnclosedVolumeWith({0, 0, 0}, {1, 1, 1}, {0.5, 0.5, 0.5}, {1, 1, 1}),
      0.5 * 0.5 * 0.5);
}

/**
 * Create a octree filled with particles randomly distributed across a fixed box. Then, for each particle it is checked
 * whether the range neighbor finding algorithm only finds boxes that are valid.
 */
TEST_F(OctreeTest, testRangeNeighborFinding) {
  using namespace autopas;

  // Create an octree with a random particle configuration.
  std::unique_ptr<OctreeNodeInterface<ParticleFP64>> root;
  std::array<double, 3> min = {0, 0, 0}, max = {1, 1, 1};
  createRandomOctree(root, 1234, min, max, 1000);

  // Iterate all particles
  std::vector<OctreeLeafNode<ParticleFP64> *> leaves;
  root->appendAllLeaves(leaves);

  // Pick random points and see if the found boxes contain the point
  for (int pointIndex = 0; pointIndex < 1000; ++pointIndex) {
    auto randomPosition = random3D(min, max);

    for (double pseudoInteractionLength : {0.1, 1.0}) {
      auto pseudoBoxMin = utils::ArrayMath::subScalar(randomPosition, pseudoInteractionLength);
      auto pseudoBoxMax = utils::ArrayMath::addScalar(randomPosition, pseudoInteractionLength);

      auto containingBoxes = root->getLeavesInRange(pseudoBoxMin, pseudoBoxMax);
      for (OctreeLeafNode<ParticleFP64> *leaf : containingBoxes) {
        ASSERT_GE(leaf->getEnclosedVolumeWith(pseudoBoxMin, pseudoBoxMax), 0);
      }
    }
  }
}
