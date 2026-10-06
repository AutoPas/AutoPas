/**
 * @file OpenMPConfiguratorTest.cpp
 * @author S. J. Newcome
 * @date 28/09/2026
 */

#include "OpenMPConfiguratorTest.h"

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
  return {ContainerOption::linkedCells,
          1.,
          TraversalOption::lc_c08,
          LoadEstimatorOption::none,
          DataLayoutOption::aos,
          Newton3Option::enabled,
          ompKind,
          chunkSize,
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
  const OpenMPConfigurator ompConfig(config.ompKind, config.ompChunkSize);

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
    EXPECT_EQ(chunkSizeGot, ompConfig.getOMPChunkSize());
  }
#endif
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