/**
 * @file C08TraversalTest.cpp
 * @author F. Gratl
 * @date 24.05.18
 */

#include "C08TraversalTest.h"
#include "testingHelpers/commonTypedefs.h"

using ::testing::_;
using ::testing::AtLeast;

void testC08Traversal(const std::array<size_t, 3>& edgeLength) {
  MFunctor functor;
  std::vector<FPCell> cells;
  cells.resize(edgeLength[0] * edgeLength[1] * edgeLength[2]);
  autopas::Particle defaultParticle;

  GridGenerator::fillWithParticles(cells, {edgeLength[0], edgeLength[1], edgeLength[2]}, defaultParticle);
#ifdef AUTOPAS_OPENMP
  int numThreadsBefore = omp_get_max_threads();
  omp_set_num_threads(4);
#endif
  autopas::C08Traversal<FPCell, MFunctor, false, true> c08Traversal({edgeLength[0], edgeLength[1], edgeLength[2]},
                                                                    &functor);

  // every particle interacts with 13 others. Last layer of each dim is covered
  // by previous interactions
  EXPECT_CALL(functor, AoSFunctor(_, _, true))
      .Times((edgeLength[0] - 1) * (edgeLength[1] - 1) * (edgeLength[2] - 1) * 13);
  c08Traversal.traverseCellPairs(cells);
#ifdef AUTOPAS_OPENMP
  omp_set_num_threads(numThreadsBefore);
#endif
}

TEST_F(C08TraversalTest, testTraversal10x10x10) {
  std::array<size_t, 3> edgeLength = {10, 10, 10};
  testC08Traversal(edgeLength);
}

TEST_F(C08TraversalTest, testTraversal2x2x2) {
  std::array<size_t, 3> edgeLength = {2, 2, 2};
  testC08Traversal(edgeLength);
}

TEST_F(C08TraversalTest, testTraversal3x3x3) {
  std::array<size_t, 3> edgeLength = {3, 3, 3};
  testC08Traversal(edgeLength);
}

TEST_F(C08TraversalTest, testTraversal2x3x4) {
  std::array<size_t, 3> edgeLength = {2, 3, 4};
  testC08Traversal(edgeLength);
}

TEST_F(C08TraversalTest, testTraversal7x8x9) {
  std::array<size_t, 3> edgeLength = {7, 8, 9};
  testC08Traversal(edgeLength);
}

template <autopas::BlackBoxTraversalOption blackBoxTraversalOption>
class C08BasedTraversalDummy : public autopas::C08BasedTraversal<FPCell, false, true, blackBoxTraversalOption> {
 public:
  explicit C08BasedTraversalDummy(const std::array<unsigned long, 3>& dims)
      : autopas::C08BasedTraversal<FPCell, false /*useSoA*/, true /*newton3*/, blackBoxTraversalOption>(dims) {}

  autopas::TraversalOptions getTraversalType() override {
    autopas::utils::ExceptionHandler::exception("not yet implemented.");
    return autopas::TraversalOptions::dummyTraversal;
  }

  template <typename LoopBody>
  inline void c08Traversal_public(LoopBody&& loopBody) {
    this->c08Traversal(loopBody);
  }
};

TEST_F(C08TraversalTest, testOuterTraversal) {
  constexpr std::array<size_t, 3> edgeLength = {7, 8, 9};
  C08BasedTraversalDummy<autopas::outer> traversal(edgeLength);
  std::array<std::array<std::array<int, edgeLength[2]>, edgeLength[1]>, edgeLength[0]> touchableArray = {};
  traversal.c08Traversal_public([&](unsigned long x, unsigned long y, unsigned long z) { touchableArray[x][y][z]++; });

  // the upper halo is NOT traversed => "x < edgeLength[0] - 1"
  for (unsigned int x = 0; x < edgeLength[0] - 1; ++x) {
    for (unsigned int y = 0; y < edgeLength[1] - 1; ++y) {
      for (unsigned int z = 0; z < edgeLength[2] - 1; ++z) {
        bool outside = false;
        if (x < 2 or x >= edgeLength[0] - 3) outside = true;
        if (y < 2 or y >= edgeLength[1] - 3) outside = true;
        if (z < 2 or z >= edgeLength[2] - 3) outside = true;
        EXPECT_EQ(touchableArray[x][y][z], outside ? 1 : 0) << "x: " << x << ", y: " << y << ", z: " << z;
      }
    }
  }
}

TEST_F(C08TraversalTest, testOuterTraversalColors) {
  // the idea of this test is to check whether the colors are used correctly, such that no race conditions can occur.

#if defined(AUTOPAS_OPENMP)
  int previousThreadCount = autopas::autopas_get_max_threads();
  omp_set_num_threads(1);
#endif

  constexpr std::array<size_t, 3> edgeLength = {7, 8, 9};
  C08BasedTraversalDummy<autopas::outer> traversal(edgeLength);
  std::array<std::array<std::array<int, edgeLength[2]>, edgeLength[1]>, edgeLength[0]> touchableArray = {};

  unsigned long color = 0;
  traversal.c08Traversal_public([&](unsigned long x, unsigned long y, unsigned long z) {
    if (x < 2 and y < 2 and z < 2) {
      // from the lower entries we can get the colors (at least if they are traversed first)
      color = autopas::utils::ThreeDimensionalMapping::threeToOneD(x, y, z, {2, 2, 2});
    }
    touchableArray[x][y][z] = color + 1;
  });

  // the upper halo is NOT traversed => "x < edgeLength[0] - 1"
  for (unsigned int x = 0; x < edgeLength[0] - 1; ++x) {
    for (unsigned int y = 0; y < edgeLength[1] - 1; ++y) {
      for (unsigned int z = 0; z < edgeLength[2] - 1; ++z) {
        unsigned long currentColor = touchableArray[x][y][z];
        bool outside = false;
        if (x < 2 or x >= edgeLength[0] - 3) outside = true;
        if (y < 2 or y >= edgeLength[1] - 3) outside = true;
        if (z < 2 or z >= edgeLength[2] - 3) outside = true;
        if (not outside) continue;

        for (unsigned int dx = 0; dx < 2; ++dx) {
          for (unsigned int dy = 0; dy < 2; ++dy) {
            for (unsigned int dz = 0; dz < 2; ++dz) {
              if (dx == 0 and dy == 0 and dz == 0) continue;
              int neighborColor = touchableArray[x + dx][y + dy][z + dz];
              EXPECT_NE(currentColor, neighborColor)
                  << "x: " << x << ", y: " << y << ", z: " << z << "; dx: " << dx << ", dy: " << dy << ", dz: " << dz;
            }
          }
        }
      }
    }
  }
#if defined(AUTOPAS_OPENMP)
  omp_set_num_threads(previousThreadCount);
#endif
}

TEST_F(C08TraversalTest, testInnerTraversal) {
  constexpr std::array<size_t, 3> edgeLength = {7, 8, 9};
  C08BasedTraversalDummy<autopas::inner> traversal(edgeLength);
  std::array<std::array<std::array<int, edgeLength[2]>, edgeLength[1]>, edgeLength[0]> touchableArray = {};
  traversal.c08Traversal_public([&](unsigned long x, unsigned long y, unsigned long z) { touchableArray[x][y][z]++; });
  for (unsigned int x = 0; x < edgeLength[0]; ++x) {
    for (unsigned int y = 0; y < edgeLength[1]; ++y) {
      for (unsigned int z = 0; z < edgeLength[2]; ++z) {
        bool outside = false;
        if (x < 2 or x >= edgeLength[0] - 3) outside = true;
        if (y < 2 or y >= edgeLength[1] - 3) outside = true;
        if (z < 2 or z >= edgeLength[2] - 3) outside = true;
        EXPECT_EQ(touchableArray[x][y][z], outside ? 0 : 1) << "x: " << x << ", y: " << y << ", z: " << z;
      }
    }
  }
}