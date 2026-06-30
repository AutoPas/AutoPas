/**
 * @file CellFunctorSortingDataTest.cpp
 * @author hmeyran
 * @date 30.06.26
 *  Unit tests for CellFunctor::computeSortingData
 *
 *  Uses synthetic sorted projection vectors so expected bounds can be verified
 *  by inspection, independently of any force kernel.
 */

#include "CellFunctorSortingDataTest.h"

/**
 * Build a projIdx vector from a plain list of projection values.
 * Original indices are assigned in order (irrelevant for bound computation).
 */
std::vector<std::pair<double, size_t>> makeProjIdx(std::initializer_list<double> projs) {
  std::vector<std::pair<double, size_t>> v;
  size_t idx = 0;
  for (double p : projs) v.push_back({p, idx++});
  return v;
}

/**
 * startI skips i particles whose projection lies entirely below projJ[0] - cutoff.
 * With projI = {1, 4, 7}, projJ = {5, 8}, cutoff = 3:
 *   threshold = 5 - 3 = 2 → first i where projI > 2 is index 1 (projI[1] = 4).
 */
TEST_F(CellFunctorSortingDataTest, testStartI) {
  auto projI = makeProjIdx({1.0, 4.0, 7.0});
  auto projJ = makeProjIdx({5.0, 8.0});
  std::vector<size_t> maxIdx, minIdx;

  const auto data = _cf.computeSortingData(projI, projJ, maxIdx, minIdx);

  EXPECT_EQ(data.startI, 1u);
}

/**
 * maxIndex[i] is the first j where projJ[j] > projI[i] + cutoff.
 * With projI = {0, 5, 10}, projJ = {6, 7, 11}, cutoff = 3:
 *   i=1 (proj=5): first j > 8 → j=2 (proj=11),  maxIndex[1] = 2
 *   i=2 (proj=10): first j > 13 → none,          maxIndex[2] = 3 (= nJ)
 */
TEST_F(CellFunctorSortingDataTest, testMaxIndex) {
  auto projI = makeProjIdx({0.0, 5.0, 10.0});
  auto projJ = makeProjIdx({6.0, 7.0, 11.0});
  std::vector<size_t> maxIdx, minIdx;

  _cf.computeSortingData(projI, projJ, maxIdx, minIdx);

  // startI = 1 (projI[0]=0 <= threshold=3 so skipped)
  EXPECT_EQ(maxIdx[1], 2u);
  EXPECT_EQ(maxIdx[2], 3u);
}

/**
 * minIndex[i] is the first j where projJ[j] >= projI[i] - cutoff (left-side pruning).
 * With projI = {0, 5, 10}, projJ = {6, 7, 11}, cutoff = 3:
 *   i=1 (proj=5): first j >= 2 → j=0 (proj=6),  minIndex[1] = 0
 *   i=2 (proj=10): first j >= 7 → j=1 (proj=7), minIndex[2] = 1
 */
TEST_F(CellFunctorSortingDataTest, testMinIndex) {
  auto projI = makeProjIdx({0.0, 5.0, 10.0});
  auto projJ = makeProjIdx({6.0, 7.0, 11.0});
  std::vector<size_t> maxIdx, minIdx;

  _cf.computeSortingData(projI, projJ, maxIdx, minIdx);

  EXPECT_EQ(minIdx[1], 0u);
  EXPECT_EQ(minIdx[2], 1u);
}

/**
 * When projI lies entirely below projJ - cutoff no i can interact with any j.
 * projI = {0, 1}, projJ = {10, 11}, cutoff = 3 → threshold = 10-3=7,
 * both projI values <= 7, so startI = 2 = nI.
 */
TEST_F(CellFunctorSortingDataTest, testAllPruned) {
  auto projI = makeProjIdx({0.0, 1.0});
  auto projJ = makeProjIdx({10.0, 11.0});
  std::vector<size_t> maxIdx, minIdx;

  const auto data = _cf.computeSortingData(projI, projJ, maxIdx, minIdx);

  EXPECT_EQ(data.startI, 2u);
}
