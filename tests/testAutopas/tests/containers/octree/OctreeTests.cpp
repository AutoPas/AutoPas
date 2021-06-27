/**
 * @file OctreeTests.cpp
 * @author Johannes Spies
 * @date 15.04.2021
 */

#include "OctreeTests.h"

#include <autopas/molecularDynamics/LJFunctor.h>
#include <autopas/selectors/ContainerSelector.h>
#include <autopas/utils/StaticCellSelector.h>
#include <testingHelpers/commonTypedefs.h>

#include <cstdio>
#include <string>
#include <vector>

#include "autopas/containers/octree/Octree.h"
#include "autopas/containers/octree/OctreeDirection.h"
#include "autopas/containers/octree/OctreeNodeInterface.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/particles/Particle.h"

using ::testing::_;

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
      std::make_unique<OctreeLeafNode<ParticleFP64>>(min, max, nullptr, 4, 1);
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
  OctreeInnerNode<ParticleFP64> inner(min, max, nullptr, 16, 1);

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
 * Get a particle from a random distribution between min and max
 * @param min The minimum coordinate of the cube from which the position is sampled
 * @param max The maximum coordinate of the cube from which the position is sampled
 * @return A particle with no speed and id 0
 */
static autopas::ParticleFP64 getRandomlyDistributedParticle(std::array<double, 3> min, std::array<double, 3> max) {
  auto randomPosition = random3D(min, max);
  return autopas::ParticleFP64(randomPosition, {0, 0, 0}, 0);
}

/**
 * Create an octree filled with particles whose positions are created by a very bad random number generator
 * @param rootRef A reference in which the newly generated octree will be saved
 * @param min The minimum coordinate of the octree's bounding box
 * @param max The maximum coordinate of the octree's bounding box
 * @param randomParticleCount How many particles should be spawned
 */
static void createRandomOctree(std::unique_ptr<autopas::OctreeNodeInterface<autopas::ParticleFP64>> &rootRef,
                               std::array<double, 3> min, std::array<double, 3> max, int randomParticleCount) {
  using namespace autopas;
  std::unique_ptr<OctreeNodeInterface<autopas::ParticleFP64>> root =
      std::make_unique<OctreeLeafNode<ParticleFP64>>(min, max, nullptr, 16, 1);
  for (int particleIndex = 0; particleIndex < randomParticleCount; ++particleIndex) {
    auto randomParticle = getRandomlyDistributedParticle(min, max);
    root->insert(root, randomParticle);
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
  srand(1234);
  createRandomOctree(root, {0, 0, 0}, {1, 1, 1}, 1000);

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
  srand(1234);
  createRandomOctree(root, min, max, 1000);

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

/**
 * Test whether the leaf splitting behavior is correct for a leaf that cannot split.
 */
TEST_F(OctreeTest, testUnableToSplit) {
  using namespace autopas;

  // Create a small octree with a high interaction length
  std::array<double, 3> min = {}, max = {1, 1, 1};
  int unsigned treeSplitThreshold = 4;
  std::unique_ptr<OctreeNodeInterface<ParticleFP64>> root =
      std::make_unique<OctreeLeafNode<ParticleFP64>>(min, max, nullptr, treeSplitThreshold, 1.0);
  ASSERT_FALSE(root->hasChildren());

  // Insert particles
  srand(1234);
  for (int unsigned i = 0; i < 2 * treeSplitThreshold; ++i) {
    auto particle = getRandomlyDistributedParticle(min, max);
    root->insert(root, particle);

    // The node should never split because of the interaction length
    ASSERT_FALSE(root->hasChildren());
  }
}

/**
 * Test whether the leaf splitting behavior is correct for a leaf that can split.
 */
TEST_F(OctreeTest, testAbleToSplit) {
  using namespace autopas;

  // Create a small octree with a high interaction length
  std::array<double, 3> min = {}, max = {1, 1, 1};
  int unsigned treeSplitThreshold = 4;
  std::unique_ptr<OctreeNodeInterface<ParticleFP64>> root =
      std::make_unique<OctreeLeafNode<ParticleFP64>>(min, max, nullptr, treeSplitThreshold, .5);
  ASSERT_FALSE(root->hasChildren());

  // Insert particles, the first should not cause the octree to split
  srand(1234);
  for (int unsigned i = 0; i < treeSplitThreshold; ++i) {
    auto particle = getRandomlyDistributedParticle(min, max);
    root->insert(root, particle);
    ASSERT_FALSE(root->hasChildren());
  }

  // These should cause the octree to split
  for (int unsigned i = 0; i < treeSplitThreshold; ++i) {
    auto particle = getRandomlyDistributedParticle(min, max);
    root->insert(root, particle);
    ASSERT_TRUE(root->hasChildren());
  }
}

std::tuple<OctreeTest::Vector3DList, OctreeTest::Vector3DList,
           std::vector<std::tuple<unsigned long, unsigned long, double>>>
OctreeTest::calculateForcesAndPairs(autopas::ContainerOption containerOption, autopas::TraversalOption traversalOption,
                                    autopas::DataLayoutOption dataLayoutOption, autopas::Newton3Option newton3Option,
                                    size_t numParticles, size_t numHaloParticles, std::array<double, 3> boxMax,
                                    double cellSizeFactor, bool doSlightShift) {
  using namespace autopas;

  std::array<double, 3> _boxMin{0, 0, 0};
  double _cutoff{1.};
  double _cutoffsquare = _cutoff * _cutoff;
  static constexpr double _eps{1.};
  static constexpr double _sig{1.};

  // Construct container
  ContainerSelector<Molecule> selector{_boxMin, boxMax, _cutoff};
  double skin = _cutoff * 0.1;
  double interactionLength = _cutoff + skin;
  selector.selectContainer(
      containerOption, autopas::ContainerSelectorInfo{cellSizeFactor, skin, 32, autopas::LoadEstimatorOption::none});
  auto container = selector.getCurrentContainer();

  // Create a functor that is able to calculate forces
  autopas::LJFunctor<Molecule, true /*applyShift*/, false /*useMixing*/, autopas::FunctorN3Modes::Both,
                     false /*calculateGlobals*/>
      ljFunctor{_cutoff};
  ljFunctor.setParticleProperties(_eps * 24, _sig * _sig);

  // Initialize the RNG in order to be able to initialize both containers with the same particles
  srand(1234);

  // Fill a smaller portion of the octree region with particles
  auto boxCenter = utils::ArrayMath::add(_boxMin, boxMax);
  boxCenter = utils::ArrayMath::mulScalar(boxCenter, 0.5);
  for (int unsigned i = 0; i < numParticles; ++i) {
    auto position = random3D(_boxMin, boxMax);
    auto particle = Molecule(position, {0, 0, 0}, i);
    container->addParticle(particle);
  }

  // Fill the area around the octree with halo particles
  auto haloMin = utils::ArrayMath::subScalar(_boxMin, interactionLength);
  auto haloMax = utils::ArrayMath::addScalar(boxMax, interactionLength);
  for (int unsigned i = 0; i < numHaloParticles;) {
    auto position = random3D(haloMin, haloMax);
    if (!utils::inBox(position, _boxMin, boxMax)) {
      auto particle = Molecule(position, {0, 0, 0}, i);
      container->addHaloParticle(particle);
      ++i;
    }
  }

  // Obtain a compatible traversal
  auto traversal =
      autopas::utils::withStaticCellType<Molecule>(container->getParticleCellTypeEnum(), [&](auto particleCellDummy) {
        return autopas::TraversalSelector<decltype(particleCellDummy)>::generateTraversal(
            traversalOption, mockFunctor, container->getTraversalSelectorInfo(), dataLayoutOption, newton3Option);
      });

  if (doSlightShift) {
    // TODO(johannes): Shift the particles
  }

  // Specify the behavior that should be executed for each particle pair
  int unsigned numPairs = 0;
  std::vector<std::tuple<unsigned long, unsigned long, double>> particlePairs;
  bool useNewton3 = newton3Option == Newton3Option::enabled;
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, useNewton3))
      .Times(testing::AtLeast(1))
      .WillRepeatedly(testing::WithArgs<0, 1>([&](auto &i, auto &j) {
        ++numPairs;

        // Store the particle pair interaction if it is within cutoff range
        auto dr = utils::ArrayMath::sub(i.getR(), j.getR());
        double dr2 = utils::ArrayMath::dot(dr, dr);
        if (dr2 <= _cutoffsquare) {
          particlePairs.template emplace_back(i.getID(), j.getID(), dr2);
        }

        // Do, what the LJ functor would do
        ljFunctor.AoSFunctor(i, j, useNewton3);
      }));

  // non useNewton3 variant should not happen
  EXPECT_CALL(mockFunctor, AoSFunctor(_, _, not useNewton3)).Times(0);

  // Perform the traversal
  mockFunctor.initTraversal();
  container->iteratePairwise(traversal.get());
  mockFunctor.endTraversal(newton3Option);

  printf("Johannes' AoS functor called %d times\n", numPairs);

  // Obtain all calculated forces
  std::vector<std::array<double, 3>> forces(numParticles), positions(numParticles);
  for (auto it = container->begin(autopas::IteratorBehavior::owned); it.isValid(); ++it) {
    EXPECT_TRUE(it->isOwned());
    auto f = it->getF();
    auto r = it->getR();
    auto id = it->getID();
    forces.at(id) = f;
    positions.at(id) = r;
  }

  return std::make_tuple(positions, forces, particlePairs);
}

TEST_P(OctreeTest, testCustomParticleDistribution) {
  using namespace autopas;

  auto [boxMax, numParticles, numHaloParticles] = GetParam();

  // Stolen from the traversalTest
  constexpr double rel_err_tolerance = 1.0e-10;
  constexpr double rel_err_tolerance_globals = 1.0e-12;

  // Obtain a starting configuration
  auto containerOption = ContainerOption::octree;
  auto traversalOption = TraversalOption::ot_c01;
  auto dataLayoutOption = DataLayoutOption::aos;
  auto newton3Option = Newton3Option::disabled;
  auto cellSizeFactor = 1;
  auto doSlightShift = false;
  auto particleDeletionPosition =
      false;  // TODO: Change this to DeletionPosition::never or something more reasonable in general and put it back in

  // Calculate the forces using the octree
  auto [calculatedPositions, calculatedForces, calculatedPairs] =
      calculateForcesAndPairs(containerOption, traversalOption, dataLayoutOption, newton3Option, numParticles,
                              numHaloParticles, boxMax, cellSizeFactor, doSlightShift);

  // Calculate the forces using the reference implementation
  auto [referencePositions, referenceForces, referencePairs] = calculateForcesAndPairs(
      autopas::ContainerOption::linkedCells, autopas::TraversalOption::lc_c08, autopas::DataLayoutOption::aos,
      autopas::Newton3Option::enabled, numParticles, numHaloParticles, boxMax, cellSizeFactor, doSlightShift);

  // Calculate which pairs are in the set difference between the reference pairs and the calculated pairs
  std::sort(calculatedPairs.begin(), calculatedPairs.end());
  std::sort(referencePairs.begin(), referencePairs.end());
  std::vector<std::tuple<long unsigned, long unsigned, double>> diff;
  std::set_difference(referencePairs.begin(), referencePairs.end(), calculatedPairs.begin(), calculatedPairs.end(),
                      std::inserter(diff, diff.begin()));
  EXPECT_TRUE(diff.empty());

  // Check whether the forces generated by the octree match the reference forces
  for (size_t i = 0; i < numParticles; ++i) {
    // Check if calculated forces match
    for (unsigned int d = 0; d < 3; ++d) {
      double calculatedForce = calculatedForces[i][d];
      double referenceForce = referenceForces[i][d];
      EXPECT_NEAR(calculatedForce, referenceForce, std::fabs(calculatedForce * rel_err_tolerance))
          << "#p" << numParticles << " #hp" << numHaloParticles << " Particle id: " << i;
    }

    // Check if returned particle positions match
    for (unsigned int d = 0; d < 3; ++d) {
      double calculatedPosition = calculatedPositions[i][d];
      double referencePosition = referencePositions[i][d];
      EXPECT_NEAR(calculatedPosition, referencePosition, std::fabs(calculatedPosition * rel_err_tolerance))
          << "#p" << numParticles << " #hp" << numHaloParticles << " Particle id: " << i;
    }

    for (std::tuple<long unsigned, long unsigned, double> &tup : diff) {
      // Get the index of the missing particle
      long unsigned missing;
      if (std::get<0>(tup) == i) {
        missing = std::get<1>(tup);
      } else if (std::get<1>(tup) == i) {
        missing = std::get<0>(tup);
      } else {
        continue;
      }

      double dr = std::sqrt(std::get<2>(tup));
      if (missing < numParticles) {
        // The missing particle is another owned particle
        auto p1 = calculatedPositions[i];
        auto p2 = calculatedPositions[missing];
        printf("Interaction between particles %lu<->%lu (positions: %f, %f, %f, %f, %f, %f, distance: %f) is missing\n",
               i, missing, p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], dr);
      } else {
        // The missing particle is a halo particle
        printf("Interaction between particles %lu<->%lu (distance: %f) is missing\n", i, missing, dr);
      }
    }
  }
}

/**
 * Lambda to generate a readable string out of the parameters of this test.
 */
static auto toString = [](const auto &info) {
  auto [boxMax, numParticles, numHaloParticles] = info.param;
  std::stringstream resStream;
  resStream << boxMax[0] << "_" << boxMax[1] << "_" << boxMax[2] << "_NP" << std::to_string(numParticles) << "_NH"
            << std::to_string(numHaloParticles);
  std::string res = resStream.str();
  std::replace(res.begin(), res.end(), '-', '_');
  std::replace(res.begin(), res.end(), '.', '_');
  return res;
};

static std::vector<GeneratorSpec> getTestParams() {
  std::vector<GeneratorSpec> result = {};
  for (auto boxMax : std::vector<std::array<double, 3>>{{3.0, 3.0, 3.0}, {10.0, 10.0, 10.0}, {20.0, 20.0, 20.0}}) {
    for (auto numParticles : {40, 50, 60, 70, 80}) {
      for (auto numHaloParticles : {0, 1, 50}) {
        result.emplace_back(boxMax, numParticles, numHaloParticles);
      }
    }
  }
  return result;
}

INSTANTIATE_TEST_SUITE_P(OctreeTestGenerated, OctreeTest, ::testing::ValuesIn(getTestParams()), toString);
