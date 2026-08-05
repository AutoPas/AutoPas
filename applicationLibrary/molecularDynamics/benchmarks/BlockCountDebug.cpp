/**
 * @file BlockCountDebug.cpp
 * @date 2026-08-05
 *
 * Standalone debug tool for the VecPattern/sorting mechanism investigation
 * (debug/investigate-vec-patterns branch, see BachelorThesis project notes).
 *
 * Counts SIMD blocks in the SoA pair kernel for the same three paths the timing ablation compares,
 * using the AUTOPAS_LJ_HWY_BLOCK_COUNTERS instrumentation in LJFunctorHWY.h:
 *
 *   unsorted   original particle order, full j-range
 *   noprune    SortedSoAView layout, full j-range   (setDisablePairPruning(true))
 *   sorted     SortedSoAView layout, pruned window
 *
 * Two counters are reported per configuration. "visited" counts every block the j-loop reaches,
 * which is what fillJRegisters is paid on. "surviving" counts blocks passing the cutoff check,
 * which is what handleNewton3Reduction is paid on.
 *
 * Hypothesis under test: pruning removes blocks that would have exited at the cutoff check, so
 * going from noprune to sorted should cut "visited" sharply while leaving "surviving" roughly
 * flat. Reordering (unsorted to noprune) should instead leave "visited" identical and reduce
 * "surviving", because clustering packs interacting particles into fewer blocks.
 *
 * Counts are exact and deterministic for a given seed, so a single run per configuration is
 * sufficient. No timing is involved and the tool needs no machine isolation, unlike the timing
 * benchmarks. Lane count is the only hardware dependence, so results are comparable across
 * machines with the same SIMD width.
 *
 * useMixing is false here, unlike LJFunctorHWYBench. Mixing only changes how epsilon/sigma are
 * fetched per lane and affects neither loop bounds nor the cutoff check, so block counts are
 * identical either way.
 *
 * Output is CSV on stdout. Not part of the thesis-cited benchmark dataset.
 */

#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/options/SortingDirectionOption.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/utils/SortingThresholdInfoSingle.h"
#include "autopas/utils/generators/TwoCellsInteractionHitrateGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorHWY.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"

#ifndef AUTOPAS_LJ_HWY_BLOCK_COUNTERS
#error "BlockCountDebug requires -DAUTOPAS_LJ_HWY_BLOCK_COUNTERS (see benchmarks/CMakeLists.txt)."
#endif

namespace {

using MoleculeType = mdLib::MoleculeLJ;
using FMCell = autopas::FullParticleCell<MoleculeType>;
using VectorizationPattern = autopas::VectorizationPatternOption::Value;

// Matches LJFunctorHWYBench.cpp so the counts line up with the timing grid.
constexpr double kCutoff = 3.0;
constexpr double kEpsilon = 1.0;
constexpr double kSigma = 1.0;
constexpr double kCellSize = 4.0;
constexpr std::array<double, 3> kLow{0.0, 0.0, 0.0};
constexpr std::array<double, 3> kCell1High{kCellSize, kCellSize, kCellSize};
// Face-adjacent cell2, matching setupCellPair(face) in LJFunctorHWYBench.cpp.
constexpr std::array<double, 3> kCell2Low{kCellSize, 0.0, 0.0};
constexpr std::array<double, 3> kCell2High{2 * kCellSize, kCellSize, kCellSize};
constexpr std::array<double, 3> kSortingDirection{1.0, 0.0, 0.0};

constexpr unsigned kSeed = 42;

using DebugFunctor = mdLib::LJFunctorHWY<MoleculeType, /*shifting=*/true, /*useMixing=*/false,
                                         autopas::FunctorN3Modes::Both, /*calculateGlobals=*/false,
                                         /*countFLOPs=*/false>;

const std::vector<size_t> kNValues = {10, 25, 50, 75, 100, 125, 150};
const std::vector<int> kHitrates = {0, 20, 30, 40, 50, 70, 90};
const std::vector<VectorizationPattern> kVecPatterns = {
    VectorizationPattern::p1xVec, VectorizationPattern::p2xVecDiv2, VectorizationPattern::pVecDiv2x2,
    VectorizationPattern::pVecx1};

const char *patternName(VectorizationPattern p) {
  switch (p) {
    case VectorizationPattern::p1xVec:
      return "p1xVec";
    case VectorizationPattern::p2xVecDiv2:
      return "p2xVecDiv2";
    case VectorizationPattern::pVecDiv2x2:
      return "pVecDiv2x2";
    case VectorizationPattern::pVecx1:
      return "pVecx1";
  }
  return "unknown";
}

/**
 * Runs one cell-pair interaction and reports the block counts it produced.
 * @param path One of "unsorted", "noprune", "sorted".
 * @param n Particles per cell.
 * @param hitrate Interacting-particle fraction handed to the generator.
 * @param newton3 Whether newton3 is enabled.
 * @param pattern Vectorization pattern to run.
 */
void runOne(const std::string &path, size_t n, int hitrate, bool newton3, VectorizationPattern pattern) {
  FMCell cell1, cell2;
  const MoleculeType defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopas::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
      cell1, cell2, kLow, kCell1High, kCell2Low, kCell2High, n, hitrate / 100.0, kCutoff, defaultParticle, kSeed);

  DebugFunctor functor(kCutoff);
  functor.setParticleProperties(kEpsilon * 24.0, kSigma * kSigma);
  functor.setVecPattern(pattern);
  if (path == "noprune") {
    functor.setDisablePairPruning(true);
  }
  functor.initTraversal();

  autopas::internal::CellFunctor<FMCell, DebugFunctor, /*bidirectional=*/true> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::soa, newton3);
  // Threshold 0 forces the sorted path, a huge threshold keeps the unsorted one.
  cellFunctor.setSoASortingThresholds(
      autopas::SortingThresholdInfoSingle(path == "unsorted" ? std::numeric_limits<size_t>::max() : 0));

  functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);

  functor.resetBlockCounters();
  cellFunctor.processCellPair(cell1, cell2, kSortingDirection);
  const auto visited = functor.getVisitedBlocks();
  const auto surviving = functor.getSurvivingBlocks();

  functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
  functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
  functor.endTraversal(newton3);

  std::printf("%s,%zu,%d,%d,%s,%zu,%zu\n", path.c_str(), n, newton3 ? 1 : 0, hitrate, patternName(pattern), visited,
              surviving);
}

}  // namespace

int main() {
  std::printf("path,N,n3,hitrate,vecPat,visited,surviving\n");
  for (const auto &path : {std::string("unsorted"), std::string("noprune"), std::string("sorted")}) {
    for (const auto n : kNValues) {
      for (const auto hitrate : kHitrates) {
        for (const auto newton3 : {false, true}) {
          for (const auto pattern : kVecPatterns) {
            runOne(path, n, hitrate, newton3, pattern);
          }
        }
      }
    }
  }
  return 0;
}
