/**
 * @file SortedSoAViewTest.cpp
 * @date 26.06.2026
 * @author hmeyran
 */

#include "SortedSoAViewTest.h"

#include "autopas/utils/SortedSoAView.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "testingHelpers/commonTypedefs.h"

using TestFunctor = LJFunctorType<>;
using TestView = autopas::SortedSoAView<Molecule, TestFunctor>;

/**
 * Fills a cell with n particles at positions (i, 0, 0) for i in [0, n).
 */
static FMCell makeLineCell(size_t n) {
  FMCell cell;
  for (size_t i = 0; i < n; ++i) {
    Molecule p({static_cast<double>(i), 0.0, 0.0}, {0.0, 0.0, 0.0}, i, 0);
    cell.addParticle(p);
  }
  return cell;
}

/**
 * Loads a cell's SoA buffer using a temporary LJFunctor instance.
 */
static void loadSoA(FMCell &cell) {
  constexpr double cutoff = 10.0;
  TestFunctor f(cutoff);
  f.setParticleProperties(24.0, 1.0);
  f.initTraversal();
  f.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize=*/false);
  f.endTraversal(false);
}

/**
 * Particles placed at x = 3, 1, 2 (unsorted). After constructing a SortedSoAView along {1,0,0}
 * the projIdx vector must be sorted ascending by x-projection.
 */
TEST_F(SortedSoAViewTest, testSortAscending) {
  FMCell cell;
  cell.addParticle(Molecule({3.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));
  cell.addParticle(Molecule({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0));
  cell.addParticle(Molecule({2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 2, 0));
  loadSoA(cell);

  autopas::SoA<Molecule::SoAArraysType> sortedSoa;
  std::vector<std::pair<double, size_t>> projIdx;

  TestView view(cell._particleSoABuffer.constructView(), {1.0, 0.0, 0.0}, sortedSoa, projIdx);

  ASSERT_EQ(projIdx.size(), 3);
  for (size_t i = 1; i < projIdx.size(); ++i) {
    EXPECT_LE(projIdx[i - 1].first, projIdx[i].first) << "projIdx not sorted at i=" << i;
  }
  // Projection values must match x-coordinates
  EXPECT_DOUBLE_EQ(projIdx[0].first, 1.0);
  EXPECT_DOUBLE_EQ(projIdx[1].first, 2.0);
  EXPECT_DOUBLE_EQ(projIdx[2].first, 3.0);
}

/**
 * After construction, x-positions in the sorted SoA must appear in ascending sorted order,
 * and force arrays must be zero-initialised (computed attributes start at 0).
 */
TEST_F(SortedSoAViewTest, testPackingAndZeroInit) {
  FMCell cell = makeLineCell(4);  // particles at x = 0, 1, 2, 3
  loadSoA(cell);

  autopas::SoA<Molecule::SoAArraysType> sortedSoa;
  std::vector<std::pair<double, size_t>> projIdx;

  TestView view(cell._particleSoABuffer.constructView(), {1.0, 0.0, 0.0}, sortedSoa, projIdx);

  auto sortedView = view.getView();
  ASSERT_EQ(sortedView.size(), 4);

  const auto *xPtr = sortedView.template begin<Molecule::AttributeNames::posX>();
  const auto *fxPtr = sortedView.template begin<Molecule::AttributeNames::forceX>();
  const auto *fyPtr = sortedView.template begin<Molecule::AttributeNames::forceY>();
  const auto *fzPtr = sortedView.template begin<Molecule::AttributeNames::forceZ>();

  for (size_t i = 0; i < sortedView.size(); ++i) {
    // Positions packed in sorted (ascending x) order
    EXPECT_DOUBLE_EQ(xPtr[i], static_cast<double>(i)) << "x-position mismatch at sorted index " << i;
    // Force arrays zero-initialised
    EXPECT_DOUBLE_EQ(fxPtr[i], 0.0) << "forceX not zero at sorted index " << i;
    EXPECT_DOUBLE_EQ(fyPtr[i], 0.0) << "forceY not zero at sorted index " << i;
    EXPECT_DOUBLE_EQ(fzPtr[i], 0.0) << "forceZ not zero at sorted index " << i;
  }
}

/**
 * Forces written to the sorted SoA buffer must be scattered back (+=) to the correct
 * original-index slots in the source SoA buffer after scatterBack().
 */
TEST_F(SortedSoAViewTest, testScatterBack) {
  // Particles at x = 2, 0, 1 — deliberately unsorted so original and sorted indices differ.
  FMCell cell;
  cell.addParticle(Molecule({2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0, 0));
  cell.addParticle(Molecule({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0));
  cell.addParticle(Molecule({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 2, 0));
  loadSoA(cell);

  autopas::SoA<Molecule::SoAArraysType> sortedSoa;
  std::vector<std::pair<double, size_t>> projIdx;

  TestView view(cell._particleSoABuffer.constructView(), {1.0, 0.0, 0.0}, sortedSoa, projIdx);

  // Write known sentinel forces into the sorted buffer: sorted slot i gets value (i + 1).
  auto sortedView = view.getView();
  auto *fxPtr = sortedView.template begin<Molecule::AttributeNames::forceX>();
  for (size_t i = 0; i < sortedView.size(); ++i) {
    fxPtr[i] = static_cast<double>(i + 1);
  }

  view.scatterBack();

  // After scatter-back, original slot projIdx[i].second must contain (i + 1).
  const auto *origFxPtr = cell._particleSoABuffer.template begin<Molecule::AttributeNames::forceX>();
  for (size_t i = 0; i < projIdx.size(); ++i) {
    EXPECT_DOUBLE_EQ(origFxPtr[projIdx[i].second], static_cast<double>(i + 1))
        << "scatter-back mismatch for sorted index " << i << " -> original index " << projIdx[i].second;
  }
}

/**
 * Constructing a SortedSoAView on an empty SoA and calling scatterBack() must not crash.
 */
TEST_F(SortedSoAViewTest, testEmptySoA) {
  FMCell cell;
  loadSoA(cell);

  autopas::SoA<Molecule::SoAArraysType> sortedSoa;
  std::vector<std::pair<double, size_t>> projIdx;

  EXPECT_NO_THROW({
    TestView view(cell._particleSoABuffer.constructView(), {1.0, 0.0, 0.0}, sortedSoa, projIdx);
    view.scatterBack();
  });
  EXPECT_EQ(projIdx.size(), 0);
}
