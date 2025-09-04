/**
* @file PatternBenchmarkTest.cpp
 * @author J.Rief
 * @date 04.09.25
 */
#include <gtest/gtest.h>
#include "autopas/utils/PatternBenchmark.h"
/**
   *
   * Tests, if the Benchmark class can properly select optimal vectorization patterns from its benchmark results.
   */
TEST(PatternBenchmarkTest,getBenchmarkResultTest) {
  constexpr size_t benchmarkSize = autopas::PatternBenchmark::_benchmarkSize;
  std::array<autopas::VectorizationPatternOption::Value, benchmarkSize * benchmarkSize> benchmarkResults{};
  /* to test the getBenchmarkResult method we fill the patternBenchmark object with artificial optimal patterns
   * if you present the pattern benchmark results in a grid of size benchmarkSize x benchmarkSize and the x-axis is the
   * amount of particles in the first buffer and the y-axis is the amount of particles in the second buffer we fill the
   * artificial optimal patterns with: the pVecx1 pattern in the top right corner the p2xVecDiv2 pattern for the top
   * edge of the grid the pVecDiv2x2 pattern for the right edge of the grid the p1xVec pattern for the rest of the grid.
   * Then, we test if the getBenchmarkResult properly maps buffer-size-pairs outside of the benchmarkSize x
   * benchmarkSize radius to the correct corner or edge.
   */
  // fcs: number of particles in first SoABuffer
  // scs: number of particles in second SoABuffer
  for (size_t fcs = 1; fcs <= benchmarkSize; ++fcs) {
    for (size_t scs = 1; scs <= benchmarkSize; ++scs) {
      if (fcs == benchmarkSize and scs == benchmarkSize) {
        benchmarkResults[(fcs - 1) + benchmarkSize * (scs - 1)] = autopas::VectorizationPatternOption::Value::pVecx1;
      } else if (fcs == benchmarkSize) {
        benchmarkResults[(fcs - 1) + benchmarkSize * (scs - 1)] =
            autopas::VectorizationPatternOption::Value::p2xVecDiv2;
      } else if (scs == benchmarkSize) {
        benchmarkResults[(fcs - 1) + benchmarkSize * (scs - 1)] =
            autopas::VectorizationPatternOption::Value::pVecDiv2x2;
      } else {
        benchmarkResults[(fcs - 1) + benchmarkSize * (scs - 1)] = autopas::VectorizationPatternOption::Value::p1xVec;
      }
    }
  }
  autopas::PatternBenchmark patternBenchmark;
  patternBenchmark._optimalPatternsNewton3On = benchmarkResults;
  patternBenchmark._optimalPatternsNewton3Off = benchmarkResults;
  patternBenchmark._patternsCalculated = true;

  // bottom left corner

  EXPECT_EQ(patternBenchmark.getBenchmarkResult(1, 1, true), autopas::VectorizationPatternOption::Value::p1xVec);
  EXPECT_EQ(patternBenchmark.getBenchmarkResult(1, 1, false), autopas::VectorizationPatternOption::Value::p1xVec);


  //  right edge

  EXPECT_EQ(patternBenchmark.getBenchmarkResult(2 * benchmarkSize, 1, true),
            autopas::VectorizationPatternOption::Value::p2xVecDiv2);
  EXPECT_EQ(patternBenchmark.getBenchmarkResult(2 * benchmarkSize, 1, false),
            autopas::VectorizationPatternOption::Value::p2xVecDiv2);

  // upper edge

  EXPECT_EQ(patternBenchmark.getBenchmarkResult(1, 2 * benchmarkSize, true),
            autopas::VectorizationPatternOption::Value::pVecDiv2x2);
  EXPECT_EQ(patternBenchmark.getBenchmarkResult(1, 2 * benchmarkSize, false),
            autopas::VectorizationPatternOption::Value::pVecDiv2x2);

  //  top right corner

  EXPECT_EQ(patternBenchmark.getBenchmarkResult(2 * benchmarkSize, 2 * benchmarkSize, true),
            autopas::VectorizationPatternOption::Value::pVecx1);
  EXPECT_EQ(patternBenchmark.getBenchmarkResult(2 * benchmarkSize, 2 * benchmarkSize, false),
            autopas::VectorizationPatternOption::Value::pVecx1);

  // if at least one buffer has 0 particles than it should default to p1xVec pattern
  EXPECT_EQ(patternBenchmark.getBenchmarkResult(0, 1, true), autopas::VectorizationPatternOption::Value::p1xVec);
  EXPECT_EQ(patternBenchmark.getBenchmarkResult(0, 1, false), autopas::VectorizationPatternOption::Value::p1xVec);


  EXPECT_EQ(patternBenchmark.getBenchmarkResult(1, 0, true), autopas::VectorizationPatternOption::Value::p1xVec);
  EXPECT_EQ(patternBenchmark.getBenchmarkResult(1, 0, false), autopas::VectorizationPatternOption::Value::p1xVec);


  EXPECT_EQ(patternBenchmark.getBenchmarkResult(0, 0, true), autopas::VectorizationPatternOption::Value::p1xVec);
  EXPECT_EQ(patternBenchmark.getBenchmarkResult(0, 0, false), autopas::VectorizationPatternOption::Value::p1xVec);

}