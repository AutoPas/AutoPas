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
class C08BasedTraversalDummy
    : public autopas::C08BasedTraversal<FPCell, MFunctor, false, true, blackBoxTraversalOption> {
 public:
  explicit C08BasedTraversalDummy(const std::array<unsigned long, 3>& dims, MFunctor* pairwiseFunctor)
      : autopas::C08BasedTraversal<FPCell, MFunctor, false /*useSoA*/, true /*newton3*/, blackBoxTraversalOption>(
            dims, pairwiseFunctor) {}

  autopas::TraversalOptions getTraversalType() override {
    autopas::utils::ExceptionHandler::exception("not yet implemented.");
    return autopas::TraversalOptions::dummyTraversal;
  }
  FRIEND_TEST(C08TraversalTest, testOuterTraversal);
};

TEST_F(C08TraversalTest, testOuterTraversal) {
  constexpr std::array<size_t, 3> edgeLength = {7, 8, 9};
  MFunctor functor;
  C08BasedTraversalDummy<autopas::outer> traversal(edgeLength, &functor);
  std::array<std::array<std::array<int, edgeLength[2]>, edgeLength[1]>, edgeLength[0]> touchableArray = {};
  traversal.c08TraversalOuter([&](unsigned long x, unsigned long y, unsigned long z) { touchableArray[x][y][z]++; });
  for (unsigned int x = 0; x < edgeLength[0]; ++x) {
    for (unsigned int y = 0; y < edgeLength[1]; ++y) {
      for (unsigned int z = 0; z < edgeLength[2]; ++z) {
        bool shouldBeTrue = false;
        if (x < 3 or x >= edgeLength[0] - 2) shouldBeTrue = true;
        if (y < 3 or y >= edgeLength[1] - 2) shouldBeTrue = true;
        if (z < 3 or z >= edgeLength[2] - 2) shouldBeTrue = true;
        EXPECT_EQ(touchableArray[x][y][z], shouldBeTrue ? 1 : 0) << "x: " << x << ", y: " << y << ", z: " << z;
      }
    }
  }
}