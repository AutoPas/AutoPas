/**
 * @file RuntimeScheduledLoopCountsTest.cpp
 * @author S. J. Newcome
 * @date 08/10/2026
 */

#include "RuntimeScheduledLoopCountsTest.h"

#include <atomic>
#include <limits>
#include <optional>
#include <set>
#include <type_traits>
#include <vector>

#include "autopas/containers/CompatibleTraversals.h"
#include "autopas/tuning/selectors/ContainerSelector.h"
#include "autopas/tuning/selectors/TraversalSelector.h"
#include "autopas/utils/OpenMPConfigurator.h"
#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/checkFunctorType.h"
#include "autopas/utils/generators/UniformGenerator.h"
#include "testingHelpers/EmptyPairwiseFunctor.h"
#include "testingHelpers/EmptyTriwiseFunctor.h"
#include "testingHelpers/GenerateValidConfigurations.h"
#include "testingHelpers/NumThreadGuard.h"
#include "testingHelpers/commonTypedefs.h"

// This test asks the OpenMP runtime how many iterations its loops have, through the OpenMP tools interface (OMPT).
// For this, the test binary has to contain an OMPT tool, which is only compiled if:
// - Clang and libomp >= 19 are used.
// - The thread sanitizer is not used. An OpenMP runtime only loads one tool, and prefers one in the binary. With the
//   thread sanitizer, it has to load Archer instead, which tells the sanitizer about OpenMP's synchronization.
#if defined(AUTOPAS_USE_OPENMP) && defined(__clang__) && __has_include(<omp-tools.h>)
#if __clang_major__ >= 19 && !__has_feature(thread_sanitizer)
#define AUTOPAS_TEST_OMPT_LOOP_COUNTS
#include <omp-tools.h>
#endif
#endif

#ifdef AUTOPAS_TEST_OMPT_LOOP_COUNTS
// The OMPT tool. An OMPT tool consists of functions which the OpenMP runtime calls, so their names and signatures are
// prescribed by the OpenMP specification (see omp-tools.h). Hence, they have parameters this tool has no use for.
//
// As the runtime finds the tool by the name of ompt_start_tool(), it loads it for all tests of the binary. Outside of
// the test below, the only effect is that recordLoopCount() is called and returns immediately.
namespace {
/**
 * Whether the loops the OpenMP runtime reports are currently recorded.
 */
std::atomic<bool> recordLoops{false};

/**
 * Number of iterations of the recorded loops, in the order they were started.
 */
std::vector<size_t> recordedLoopCounts;

/**
 * The callback the runtime calls in every thread at the begin and at the end of every worksharing construct (loops,
 * single, sections, ...).
 *
 * Records the number of iterations of the loops which the runtime reports to this function, which the runtime calls
 * when an event occurs such as a worksharing construct beginning or ending. This only holds for schedules with a kind
 * other than static, dynamic, or guided. The test uses trapezoidal.
 *
 * @param workType The kind of worksharing construct, for loops including the schedule kind.
 * @param endpoint Whether the construct begins or ends.
 * @param count Number of iterations of the loop.
 */
void recordLoopCount(ompt_work_t workType, ompt_scope_endpoint_t endpoint, ompt_data_t * /*parallelData*/,
                     ompt_data_t * /*taskData*/, uint64_t count, const void * /*codeptr*/) {
  // Every thread of the team reports the loop, so only take the report of thread 0.
  if (recordLoops and workType == ompt_work_loop_other and endpoint == ompt_scope_begin and
      autopas::autopas_get_thread_num() == 0) {
    recordedLoopCounts.push_back(count);
  }
}

/**
 * The runtime calls this once when it loads the tool.
 *
 * Registers recordLoopCount(). The functions of the runtime a tool may call are not linkable symbols. Instead, the
 * runtime passes a lookup function, which returns them by name as generic function pointers that have to be cast.
 *
 * @param lookup Returns the functions of the runtime by name.
 * @return Non-zero, which tells the runtime to keep the tool active.
 */
int initializeTool(ompt_function_lookup_t lookup, int /*initialDeviceNum*/, ompt_data_t * /*toolData*/) {
  const auto setCallback = reinterpret_cast<ompt_set_callback_t>(lookup("ompt_set_callback"));
  // ompt_set_callback registers callbacks of all kinds, so it takes a generic function pointer as well.
  setCallback(ompt_callback_work, reinterpret_cast<ompt_callback_t>(&recordLoopCount));
  return 1;
}

/**
 * The runtime calls this when it shuts down, for tools to e.g. write their results.
 *
 * This tool has nothing to clean up. The function has to exist nonetheless, as the runtime calls it unconditionally.
 */
void finalizeTool(ompt_data_t * /*toolData*/) {}

}  // namespace

/**
 * Entry point of the OMPT tool. When the OpenMP runtime initializes, it looks for a function with exactly this name, so
 * it must not be in a namespace and must have C linkage.
 * @return The functions to initialize and to finalize the tool, and a value the runtime passes to both, which this tool
 * does not use.
 */
extern "C" ompt_start_tool_result_t *ompt_start_tool(unsigned int /*ompVersion*/, const char * /*runtimeVersion*/) {
  static ompt_start_tool_result_t result{&initializeTool, &finalizeTool, {0}};
  return &result;
}
#endif

/**
 * For an arbitrary valid configuration of every container, traversal and cell size factor, runs the traversal on a
 * filled container and checks that the loop counts it reports are those of the loops the OpenMP runtime scheduled.
 */
TEST_F(RuntimeScheduledLoopCountsTest, testLoopCountsMatchOpenMPRuntime) {
#ifndef AUTOPAS_TEST_OMPT_LOOP_COUNTS
  GTEST_SKIP() << "Needs the OpenMP tools interface (OMPT) of Clang and libomp >= 19, without the thread sanitizer.";
#else
  // Restores the OpenMP schedule at the end of the test.
  struct ScheduleGuard {
    ScheduleGuard() { autopas::autopas_get_schedule(&kind, &chunkSize); }
    ~ScheduleGuard() { autopas::autopas_set_schedule(kind, chunkSize); }
    omp_sched_t kind{};
    int chunkSize{};
  } scheduleGuard;
  // The kind recordLoopCount() recognizes the loops with schedule(runtime) by.
  autopas::autopas_set_schedule(autopas::OpenMPConfigurator(autopas::OpenMPKindOption::lb4omp_trapezoidal, 1));
  const NumThreadGuard numThreadGuard(2);

  // Runs the given code and returns the loop counts of the loops with schedule(runtime) in it.
  const auto measureLoopCounts = [](const auto &code) {
    recordedLoopCounts.clear();
    recordLoops = true;
    code();
    recordLoops = false;
    return recordedLoopCounts;
  };

  // Not every OpenMP runtime loads the tool, e.g. LB4OMP is built without OMPT. So check with a loop of known length.
  const auto probedLoopCounts = measureLoopCounts([]() {
    AUTOPAS_OPENMP(parallel for schedule(runtime))
    for (int i = 0; i < 7; ++i) {
    }
  });
  if (probedLoopCounts != std::vector<size_t>{7}) {
    GTEST_SKIP() << "The OpenMP runtime does not report the loops with schedule(runtime) through OMPT.";
  }

  constexpr std::array<double, 3> boxMin{0., 0., 0.};
  constexpr std::array<double, 3> boxMax{8., 8., 8.};
  constexpr double cutoff = 1.;
  constexpr double skin = 0.2;
  constexpr unsigned int clusterSize = 4;
  constexpr size_t sortingThreshold = std::numeric_limits<size_t>::max();
  constexpr size_t numParticles = 500;
  constexpr size_t numHaloParticles = 200;

  const auto static1OnlyTraversals = autopas::compatibleTraversals::allTraversalsSupportingOnlyStatic1Scheduling();
  std::set<autopas::TraversalOption> checkedTraversals;
  // Checks an arbitrary configuration of every container, traversal and cell size factor for the interaction type of
  // the given functor.
  const auto checkLoopCounts = [&](auto &functor) {
    using Functor = std::remove_reference_t<decltype(functor)>;
    const autopas::InteractionTypeOption interactionType = autopas::utils::isPairwiseFunctor<Functor>()
                                                               ? autopas::InteractionTypeOption::pairwise
                                                               : autopas::InteractionTypeOption::triwise;
    for (const auto &containerOption : autopas::ContainerOption::getAllOptions()) {
      for (const auto &traversalOption : autopas::TraversalOption::getAllOptions()) {
        for (const double cellSizeFactor : {0.5, 1., 1.5}) {
          // The loop counts do not depend on the other components of the configuration.
          const auto config =
              getArbitraryConfiguration(interactionType, containerOption, traversalOption, std::nullopt, std::nullopt,
                                        std::nullopt, cellSizeFactor, std::nullopt, std::nullopt, std::nullopt,
                                        /*throwIfNone*/ false);
          if (not config) {
            continue;
          }
          SCOPED_TRACE(config->toString());

          const autopas::ContainerSelectorInfo containerInfo{
              boxMin,      boxMax,           cutoff,           config->cellSizeFactor, skin,
              clusterSize, sortingThreshold, sortingThreshold, config->loadEstimator};
          auto container = autopas::ContainerSelector<Molecule>::generateContainer(config->container, containerInfo);
          autopas::generators::UniformGenerator::fillWithParticles(*container, Molecule({0., 0., 0.}, {0., 0., 0.}, 0),
                                                                   container->getBoxMin(), container->getBoxMax(),
                                                                   numParticles);
          autopas::generators::UniformGenerator::fillWithHaloParticles(
              *container, Molecule({0., 0., 0.}, {0., 0., 0.}, numParticles /*initial ID*/), container->getCutoff(),
              numHaloParticles);

          const auto traversal = autopas::TraversalSelector::generateTraversalFromConfig<Molecule, Functor>(
              *config, functor, container->getTraversalSelectorInfo());
          if (not traversal) {
            // The traversal is not applicable to the domain.
            continue;
          }

          // As in LogicHandler: from the container holding the particles, before the neighbor lists are rebuilt.
          const auto reportedLoopCounts = traversal->getRuntimeScheduledLoopCounts(
              container->getTraversalSelectorInfo(),
              container->getNumberOfParticles(autopas::IteratorBehavior::ownedOrHalo));

          container->rebuildNeighborLists(traversal.get());
          const auto measuredLoopCounts = measureLoopCounts([&]() { container->computeInteractions(traversal.get()); });

          if (not reportedLoopCounts) {
            if (not static1OnlyTraversals.contains(config->traversal)) {
              EXPECT_TRUE(measuredLoopCounts.empty())
                  << "The traversal has loops with schedule(runtime), but does not report loop counts.";
            }
            continue;
          }
          checkedTraversals.insert(config->traversal);

          if (config->traversal == autopas::TraversalOption::vl_list_iteration_c27) {
            // This traversal reports upper bounds.
            ASSERT_EQ(reportedLoopCounts->size(), measuredLoopCounts.size());
            for (size_t i = 0; i < measuredLoopCounts.size(); ++i) {
              EXPECT_GE(reportedLoopCounts->at(i), measuredLoopCounts[i]) << "Loop " << i;
            }
          } else {
            EXPECT_EQ(*reportedLoopCounts, measuredLoopCounts);
          }
        }
      }
    }
  };

  EmptyPairwiseFunctor<Molecule> pairwiseFunctor;
  checkLoopCounts(pairwiseFunctor);

  EmptyTriwiseFunctor<Molecule> triwiseFunctor;
  checkLoopCounts(triwiseFunctor);

  for (const auto &traversalOption : autopas::TraversalOption::getAllOptions()) {
    if (not static1OnlyTraversals.contains(traversalOption)) {
      EXPECT_TRUE(checkedTraversals.contains(traversalOption))
          << traversalOption.to_string() << " did not report loop counts for any configuration.";
    }
  }
#endif
}