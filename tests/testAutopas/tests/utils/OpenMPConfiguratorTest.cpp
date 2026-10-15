/**
 * @file OpenMPConfiguratorTest.cpp
 * @author S. J. Newcome
 * @date 28/09/2026
 */

#include "OpenMPConfiguratorTest.h"

#include <utility>
#include <vector>

#include "autopas/utils/OpenMPConfigurator.h"
#include "autopas/utils/WrapOpenMP.h"

using namespace autopas;

namespace {

/**
 * Saves the OpenMP schedule on construction and restores it on destruction, so that tests do not leave their schedule
 * behind for other tests.
 */
class ScheduleGuard final {
 public:
  /**
   * Saves the current OpenMP schedule.
   */
  ScheduleGuard() { autopas_get_schedule(&_kindBefore, &_chunkSizeBefore); }

  /**
   * Restores the saved OpenMP schedule.
   */
  ~ScheduleGuard() { autopas_set_schedule(_kindBefore, _chunkSizeBefore); }

  /**
   * delete copy constructor.
   */
  ScheduleGuard(const ScheduleGuard &) = delete;

  /**
   * delete copy assignment operator.
   * @return deleted, so not important
   */
  ScheduleGuard &operator=(const ScheduleGuard &) = delete;

 private:
  omp_sched_t _kindBefore{};
  int _chunkSizeBefore{};
};

}  // namespace

Configuration OpenMPConfiguratorTest::getConfiguration(OpenMPKindOption ompKind) {
  OpenMPChunkSize ompChunkSize = 4;
  if (ompKind == OpenMPKindOption::omp_auto) {
    ompChunkSize = 1;
  } else if (OpenMPKindOption::getAuto4OMPSelectionMethods().contains(ompKind)) {
    ompChunkSize = OpenMPChunkSize::golden;
  }
  return {ContainerOption::linkedCells,
          1.,
          TraversalOption::lc_c08,
          LoadEstimatorOption::none,
          DataLayoutOption::aos,
          Newton3Option::enabled,
          ompKind,
          ompChunkSize,
          InteractionTypeOption::pairwise,
          VectorizationPatternOption::NA};
}

/**
 * For the schedule kind given by the parameter:
 * - If Configuration::hasCompatibleValues() accepts it, checks that it can be set, the same way LogicHandler sets it,
 *   and that the OpenMP runtime then reports it back.
 * - Otherwise, checks that trying to set it fails.
 */
TEST_P(OpenMPConfiguratorTest, setScheduleTest) {
  const auto ompKind = GetParam();
  // Sanity check that only the OpenMP schedule kind can make the configuration incompatible.
  ASSERT_TRUE(getConfiguration(OpenMPKindOption::omp_static).hasCompatibleValues());

  const auto config = getConfiguration(ompKind);
  const OpenMPConfigurator ompConfig(config.ompKind, config.ompChunkSize, {}, autopas_get_max_threads());

  if (not config.hasCompatibleValues()) {
    // Kinds the OpenMP runtime does not support must never be set. Should one get past
    // Configuration::hasCompatibleValues(), setting it throws.
    EXPECT_ANY_THROW(autopas_set_schedule(ompConfig));
    return;
  }

#ifndef AUTOPAS_USE_OPENMP
  GTEST_SKIP() << "Without OpenMP, there is no schedule to set.";
#else
  const ScheduleGuard scheduleGuard;
  autopas_set_schedule(ompConfig);

  omp_sched_t kind{};
  int chunkSizeGot{};
  autopas_get_schedule(&kind, &chunkSizeGot);
  EXPECT_EQ(kind, ompConfig.getOMPKind());
  // OpenMP runtimes ignore the chunk size for auto.
  if (ompKind != OpenMPKindOption::omp_auto) {
    EXPECT_EQ(chunkSizeGot, ompConfig.getActualOMPChunkSize());
  }
#endif
}

/**
 * Checks that Configuration::hasCompatibleValues() only accepts the balanced chunk size with the static schedule kind.
 */
TEST_P(OpenMPConfiguratorTest, balancedOnlyWithStaticTest) {
  const auto ompKind = GetParam();
  auto config = getConfiguration(ompKind);
  config.ompChunkSize = OpenMPChunkSize::balanced;
  EXPECT_EQ(config.hasCompatibleValues(), ompKind == OpenMPKindOption::omp_static);
}

INSTANTIATE_TEST_SUITE_P(Generated, OpenMPConfiguratorTest, ::testing::ValuesIn(OpenMPKindOption::getAllOptions()),
                         [](const ::testing::TestParamInfo<OpenMPKindOption> &info) { return info.param.to_string(); });

/**
 * Checks that the OpenMP runtime itself rejects the schedule kinds above the schedule level AutoPas is built with,
 * i.e. that it is necessary to exclude them. As OpenMPConfigurator refuses to produce these kinds, they are passed to
 * the runtime directly.
 *
 * The runtimes handle invalid kinds differently:
 * - libgomp (used with GCC) ignores them, keeping the previous schedule.
 * - libomp (used with Clang) warns and falls back to static.
 * - libomp's static_steal is not checked at the clangLibomp level, as libomp knows it but does not set it correctly.
 */
TEST(OpenMPRuntimeTest, runtimeRejectsKindsAboveScheduleLevel) {
#ifndef AUTOPAS_USE_OPENMP
  GTEST_SKIP() << "Without OpenMP, there is no runtime to reject anything.";
#else
  // The first of LB4OMP's scheduling techniques, omp_sched_fsc in LB4OMP's omp.h.
  // 103 comes from LB4OMP itself.
  // All the kinds from this until kmp_sched_upper are LB4OMP kinds. We don't check these for maintainence purposes: if
  // a new kind gets added to LB4OMP, the test is still valid; if a new kind gets added to libomp, this will correctly
  // fail as 103 is expected only to work with LB4OMP builds and not standard libomp builds.
  const auto lb4ompFsc = static_cast<omp_sched_t>(103);

  std::vector<omp_sched_t> unsupportedKinds;
  switch (OpenMPKindOption::runtimeScheduleLevel) {
    case OpenMPKindOption::ScheduleLevel::gccLibgomp:
      unsupportedKinds = {omp_sched_trapezoidal, omp_sched_static_steal, lb4ompFsc};
      break;
    case OpenMPKindOption::ScheduleLevel::clangLibomp:
      // static-steal is incorrectly supported at this level. We exclude it from the test as libomp does not correctly
      // reject it.
    case OpenMPKindOption::ScheduleLevel::clangFixedLibomp:
      unsupportedKinds = {lb4ompFsc};
      break;
    case OpenMPKindOption::ScheduleLevel::clangLB4OMP:
      GTEST_SKIP() << "LB4OMP supports all schedule kinds.";
  }

  const ScheduleGuard scheduleGuard;
  for (const auto unsupportedKind : unsupportedKinds) {
    // Start from a schedule every runtime supports, so that the rejection is noticed however the runtime handles it.
    autopas_set_schedule(omp_sched_dynamic, 7);
    autopas_set_schedule(unsupportedKind, 4);

    omp_sched_t kind{};
    int chunkSize{};
    autopas_get_schedule(&kind, &chunkSize);
    EXPECT_NE(kind, unsupportedKind) << "The OpenMP runtime accepted the schedule kind " << unsupportedKind
                                     << ", although schedule level "
                                     << static_cast<int>(OpenMPKindOption::runtimeScheduleLevel) << " excludes it.";
  }
#endif
}

/**
 * Checks OpenMPConfigurator::fallsBackToOtherKindForAllLoops() for a single loop, right at the thresholds of libomp's
 * and LB4OMP's criteria. With libgomp, nothing falls back.
 */
TEST(OpenMPConfiguratorFallbackTest, fallsBackToOtherKindForSingleLoop) {
  const bool hasFallbacks = OpenMPKindOption::runtimeScheduleLevel != OpenMPKindOption::ScheduleLevel::gccLibgomp;

  // static_steal needs at least one chunk per thread: 29 iterations are 8 chunks of size 4, 28 iterations are 7.
  EXPECT_FALSE(OpenMPConfigurator(OpenMPKindOption::lb4omp_static_steal, 4, {29}, 8).fallsBackToOtherKindForAllLoops());
  EXPECT_EQ(OpenMPConfigurator(OpenMPKindOption::lb4omp_static_steal, 4, {28}, 8).fallsBackToOtherKindForAllLoops(),
            hasFallbacks);
  EXPECT_EQ(
      OpenMPConfigurator(OpenMPKindOption::lb4omp_static_steal, 4, {1000000}, 1).fallsBackToOtherKindForAllLoops(),
      hasFallbacks);

  // guided needs (2 * chunk + 1) * numThreads < loopCount: (2 * 4 + 1) * 8 = 72.
  EXPECT_FALSE(OpenMPConfigurator(OpenMPKindOption::omp_guided, 4, {73}, 8).fallsBackToOtherKindForAllLoops());
  EXPECT_EQ(OpenMPConfigurator(OpenMPKindOption::omp_guided, 4, {72}, 8).fallsBackToOtherKindForAllLoops(),
            hasFallbacks);
  EXPECT_EQ(OpenMPConfigurator(OpenMPKindOption::omp_guided, 4, {1000000}, 1).fallsBackToOtherKindForAllLoops(),
            hasFallbacks);

  // auto is guided with a chunk size of 1, whichever chunk size is configured: (2 * 1 + 1) * 8 = 24.
  EXPECT_FALSE(OpenMPConfigurator(OpenMPKindOption::omp_auto, 4, {25}, 8).fallsBackToOtherKindForAllLoops());
  EXPECT_EQ(OpenMPConfigurator(OpenMPKindOption::omp_auto, 4, {24}, 8).fallsBackToOtherKindForAllLoops(), hasFallbacks);

  // Other kinds never fall back.
  for (const OpenMPKindOption kind : {OpenMPKindOption::omp_static, OpenMPKindOption::omp_dynamic,
                                      OpenMPKindOption::llvm_trapezoidal, OpenMPKindOption::lb4omp_fac2a}) {
    EXPECT_FALSE(OpenMPConfigurator(kind, 1000, {1}, 8).fallsBackToOtherKindForAllLoops()) << kind.to_string();
    EXPECT_FALSE(OpenMPConfigurator(kind, 1000, {1}, 1).fallsBackToOtherKindForAllLoops()) << kind.to_string();
  }
}

/**
 * Checks that OpenMPConfigurator::fallsBackToOtherKindForAllLoops() only holds if every loop falls back.
 */
TEST(OpenMPConfiguratorFallbackTest, fallsBackToOtherKindForAllLoops) {
  const bool hasFallbacks = OpenMPKindOption::runtimeScheduleLevel != OpenMPKindOption::ScheduleLevel::gccLibgomp;

  // With 8 threads, guided with chunk size 4 falls back for loops of up to 72 iterations.
  EXPECT_EQ(OpenMPConfigurator(OpenMPKindOption::omp_guided, 4, {10, 72, 30}, 8).fallsBackToOtherKindForAllLoops(),
            hasFallbacks);
  EXPECT_FALSE(OpenMPConfigurator(OpenMPKindOption::omp_guided, 4, {10, 73, 30}, 8).fallsBackToOtherKindForAllLoops())
      << "One loop does not fall back.";
  EXPECT_FALSE(OpenMPConfigurator(OpenMPKindOption::omp_guided, 4, {}, 8).fallsBackToOtherKindForAllLoops())
      << "Without loops, nothing falls back.";
}

/**
 * Checks OpenMPConfigurator::getGoldenChunkSize() against values calculated by hand from LB4OMP's formula:
 * loopCount / (2^(floor(log2(loopCount / numThreads) * 0.618) + 1) * numThreads), with integer divisions.
 */
TEST(OpenMPConfiguratorChunkSizeTest, getGoldenChunkSize) {
  // floor(log2(1000) * 0.618) = floor(6.16) = 6, so 8000 / (2^7 * 8) = 7
  EXPECT_EQ(OpenMPConfigurator::getGoldenChunkSize(8000, 8), 7);
  // floor(log2(27) * 0.618) = floor(2.94) = 2, so 27 / (2^3 * 1) = 3
  EXPECT_EQ(OpenMPConfigurator::getGoldenChunkSize(27, 1), 3);
  // floor(log2(8) * 0.618) = floor(1.85) = 1, so 64 / (2^2 * 8) = 2
  EXPECT_EQ(OpenMPConfigurator::getGoldenChunkSize(64, 8), 2);
  // floor(log2(7) * 0.618) = floor(1.73) = 1, so 63 / (2^2 * 8) = 1
  EXPECT_EQ(OpenMPConfigurator::getGoldenChunkSize(63, 8), 1);
  // With fewer than two iterations per thread, LB4OMP's formula gives 0, which is replaced by 1.
  EXPECT_EQ(OpenMPConfigurator::getGoldenChunkSize(15, 8), 1);
  EXPECT_EQ(OpenMPConfigurator::getGoldenChunkSize(7, 8), 1);
  EXPECT_EQ(OpenMPConfigurator::getGoldenChunkSize(0, 8), 1);
}

/**
 * Checks that OpenMPConfigurator::getActualOMPChunkSize() converts every type of chunk size to what the OpenMP runtime
 * takes.
 */
TEST(OpenMPConfiguratorChunkSizeTest, getActualOMPChunkSize) {
  // With 8 threads, the golden chunk sizes of these loops are 7, 2 and 24 (= 100000 / (2^9 * 8)).
  const std::vector<size_t> loopCounts{8000, 64, 100000};

  EXPECT_EQ(OpenMPConfigurator(OpenMPKindOption::omp_dynamic, 4, loopCounts, 8).getActualOMPChunkSize(), 4);
  EXPECT_ANY_THROW(OpenMPConfigurator(OpenMPKindOption::omp_dynamic, 0, loopCounts, 8).getActualOMPChunkSize())
      << "A fixed chunk size must be at least 1.";
  EXPECT_ANY_THROW(OpenMPConfigurator(OpenMPKindOption::sequential, OpenMPChunkSize::notApplicable, loopCounts, 8)
                       .getActualOMPChunkSize())
      << "The not-applicable chunk size has no actual chunk size.";
  EXPECT_EQ(OpenMPConfigurator(OpenMPKindOption::omp_static, OpenMPChunkSize::balanced, loopCounts, 8)
                .getActualOMPChunkSize(),
            0);
  EXPECT_EQ(
      OpenMPConfigurator(OpenMPKindOption::omp_dynamic, OpenMPChunkSize::golden, loopCounts, 8).getActualOMPChunkSize(),
      2)
      << "The golden chunk size should be the smallest golden chunk size of all loops.";
  EXPECT_EQ(OpenMPConfigurator(OpenMPKindOption::omp_dynamic, OpenMPChunkSize::golden, {}, 8).getActualOMPChunkSize(),
            1)
      << "Without loop counts, the golden chunk size should be 1.";
  // For auto and Auto4OMP's selection methods, the kind determines the chunk size.
  const std::vector<std::pair<OpenMPKindOption, int>> autoKindsAndChunkSizes{
      {OpenMPKindOption::omp_auto, 1},
      {OpenMPKindOption::auto4omp_randomSel, 2},
      {OpenMPKindOption::auto4omp_exhaustiveSel, 3},
      {OpenMPKindOption::auto4omp_binarySearch, 4},
      {OpenMPKindOption::auto4omp_expertSel, 5},
  };
  for (const auto &[kind, expectedChunkSize] : autoKindsAndChunkSizes) {
    EXPECT_EQ(OpenMPConfigurator(kind, OpenMPChunkSize::golden, loopCounts, 8).getActualOMPChunkSize(),
              expectedChunkSize)
        << kind.to_string();
  }
}

/**
 * Checks that Configuration::hasCompatibleValues() only accepts auto with a chunk size of 1, and Auto4OMP's selection
 * methods with the golden chunk size. Kinds the OpenMP runtime does not support are skipped.
 */
TEST(OpenMPConfiguratorChunkSizeTest, autoKindsOnlyWithTheirChunkSize) {
  const std::vector<OpenMPChunkSize> chunkSizes{1, OpenMPChunkSize::golden, 4};

  const OpenMPKindOption autoKind{OpenMPKindOption::omp_auto};
  if (autoKind.isSupportedAtScheduleLevel(OpenMPKindOption::runtimeScheduleLevel)) {
    auto config = OpenMPConfiguratorTest::getConfiguration(autoKind);
    for (const auto chunkSize : chunkSizes) {
      config.ompChunkSize = chunkSize;
      EXPECT_EQ(config.hasCompatibleValues(), chunkSize == OpenMPChunkSize(1)) << "auto, " << chunkSize.to_string();
    }
  }

  for (const auto ompKind : OpenMPKindOption::getAuto4OMPSelectionMethods()) {
    if (not ompKind.isSupportedAtScheduleLevel(OpenMPKindOption::runtimeScheduleLevel)) {
      continue;
    }
    auto config = OpenMPConfiguratorTest::getConfiguration(ompKind);
    for (const auto chunkSize : chunkSizes) {
      config.ompChunkSize = chunkSize;
      EXPECT_EQ(config.hasCompatibleValues(), chunkSize == OpenMPChunkSize::golden)
          << ompKind.to_string() << ", " << chunkSize.to_string();
    }
  }
}

/**
 * Checks OpenMPConfigurator::fallsBackToChunkSize1() right at the thresholds for the golden and the balanced chunk
 * size.
 */
TEST(OpenMPConfiguratorFallbackTest, fallsBackToChunkSize1) {
  // With 8 threads, the golden chunk size is 1 for up to 63 iterations (see the test getGoldenChunkSize).
  EXPECT_TRUE(
      OpenMPConfigurator(OpenMPKindOption::omp_dynamic, OpenMPChunkSize::golden, {63}, 8).fallsBackToChunkSize1());
  EXPECT_FALSE(
      OpenMPConfigurator(OpenMPKindOption::omp_dynamic, OpenMPChunkSize::golden, {64}, 8).fallsBackToChunkSize1());
  EXPECT_TRUE(OpenMPConfigurator(OpenMPKindOption::omp_dynamic, OpenMPChunkSize::golden, {64, 63, 8000}, 8)
                  .fallsBackToChunkSize1());
  EXPECT_TRUE(
      OpenMPConfigurator(OpenMPKindOption::omp_dynamic, OpenMPChunkSize::golden, {}, 8).fallsBackToChunkSize1());

  // balanced is the same as a chunk size of 1 if no loop has more iterations than there are threads.
  EXPECT_TRUE(
      OpenMPConfigurator(OpenMPKindOption::omp_static, OpenMPChunkSize::balanced, {8, 3}, 8).fallsBackToChunkSize1());
  EXPECT_FALSE(
      OpenMPConfigurator(OpenMPKindOption::omp_static, OpenMPChunkSize::balanced, {8, 9}, 8).fallsBackToChunkSize1());
  EXPECT_FALSE(
      OpenMPConfigurator(OpenMPKindOption::omp_static, OpenMPChunkSize::balanced, {}, 8).fallsBackToChunkSize1())
      << "Without loop counts, this cannot be known.";

  // Fixed and not-applicable chunk sizes do not fall back.
  EXPECT_FALSE(OpenMPConfigurator(OpenMPKindOption::omp_static, 1, {8}, 8).fallsBackToChunkSize1());
  EXPECT_FALSE(OpenMPConfigurator(OpenMPKindOption::omp_dynamic, 4, {8}, 8).fallsBackToChunkSize1());
  EXPECT_FALSE(
      OpenMPConfigurator(OpenMPKindOption::sequential, OpenMPChunkSize::notApplicable, {8}, 8).fallsBackToChunkSize1());
}

/**
 * Checks that setting static with the balanced chunk size results in static without a chunk size, which the OpenMP
 * runtimes report as a chunk size of 0.
 */
TEST(OpenMPConfiguratorChunkSizeTest, setBalancedSchedule) {
#ifndef AUTOPAS_USE_OPENMP
  GTEST_SKIP() << "Without OpenMP, there is no schedule to set.";
#else
  const ScheduleGuard scheduleGuard;
  // Start from a schedule with a chunk size, so that it is noticed if the chunk size is not reset.
  autopas_set_schedule(omp_sched_static, 7);
  autopas_set_schedule(OpenMPConfigurator(OpenMPKindOption::omp_static, OpenMPChunkSize::balanced, {100}, 4));

  omp_sched_t kind{};
  int chunkSizeGot{};
  autopas_get_schedule(&kind, &chunkSizeGot);
  EXPECT_EQ(kind, omp_sched_static);
  EXPECT_EQ(chunkSizeGot, 0);
#endif
}