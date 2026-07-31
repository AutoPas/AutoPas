/**
 * @file PruningRatioDebug.cpp
 * @date 2026-07-31
 *
 * Standalone debug tool (Tier B of the VecPattern/sorting investigation, see BachelorThesis
 * project notes). Does NOT touch any production hot-path code: it only reads out the
 * minIndex/maxIndex bounds that CellFunctor::computeSortingData() already computes (public
 * method) before the sorted kernel runs, and sums them to get the fraction of the full n1*n2
 * space that survives the 1-D cutoff prune.
 *
 * Hypothesis under test: the sorted SoA path's pattern-independence is driven by how much of
 * the pruned window survives (fewer iterations -> less relative opportunity for per-pattern
 * register-fill overhead to show up). If so, this pruning ratio should correlate with the
 * pattern-spread-vs-hitrate trend already found in the benchmark dataset.
 *
 * Not wired into CMake by default; build manually alongside LJFunctorHWYBench (see CMakeLists.txt
 * in this directory). Not part of the thesis-cited benchmark dataset -- this is a debug/triage
 * tool for the debug/investigate-vec-patterns branch only.
 */

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/utils/SortedSoAView.h"
#include "autopas/utils/generators/TwoCellsInteractionHitrateGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorHWY.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"

namespace {

using MoleculeType = mdLib::MoleculeLJ;
using FMCell = autopas::FullParticleCell<MoleculeType>;

constexpr double kCutoff = 3.0;
constexpr double kCellSize = 4.0;
constexpr std::array<double, 3> kLow{0.0, 0.0, 0.0};
constexpr std::array<double, 3> kCell1High{kCellSize, kCellSize, kCellSize};
// Face-adjacent cell2, matching setupCellPair(face) in LJFunctorHWYBench.cpp.
constexpr std::array<double, 3> kCell2Low{kCellSize, 0.0, 0.0};
constexpr std::array<double, 3> kCell2High{2 * kCellSize, kCellSize, kCellSize};
constexpr std::array<double, 3> kSortingDirection{1.0, 0.0, 0.0};

using BenchFunctor = mdLib::LJFunctorHWY<MoleculeType, /*shifting=*/true, /*useMixing=*/false,
                                         autopas::FunctorN3Modes::Both, /*calculateGlobals=*/false,
                                         /*countFLOPs=*/false>;

BenchFunctor makeFunctor() {
  BenchFunctor f(kCutoff);
  f.setParticleProperties(1.0 * 24.0, 1.0 * 1.0);
  return f;
}

/**
 * Generates one (cell1, cell2) sample at the given N/hitrate and returns the pruning ratio:
 * sum_{i=startI}^{endI-1}(maxIndex[i] - minIndex[i]) / (n1 * n2).
 * A ratio of 1.0 means no pruning at all (equivalent to the unsorted path's full scan); a small
 * ratio means most of the n1*n2 space was skipped before the kernel ever runs.
 */
double pruningRatioSample(std::size_t n, double hitrate, unsigned seed) {
  FMCell cell1, cell2;
  const MoleculeType defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopas::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
      cell1, cell2, kLow, kCell1High, kCell2Low, kCell2High, n, hitrate, kCutoff, defaultParticle, seed);

  auto functor = makeFunctor();
  functor.initTraversal();
  functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
  functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);

  // Same CellFunctor construction as the bench file; only used here for its public
  // computeSortingData(), never for processCellPair().
  autopas::internal::CellFunctor<FMCell, BenchFunctor, /*bidirectional=*/true> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::soa, /*useNewton3=*/true);

  autopas::SoA<MoleculeType::SoAArraysType> sortedSoa1, sortedSoa2;
  std::vector<std::pair<double, std::size_t>> projIdx1, projIdx2;
  autopas::SortedSoAView<MoleculeType, BenchFunctor> view1(cell1._particleSoABuffer, kSortingDirection, sortedSoa1,
                                                           projIdx1);
  autopas::SortedSoAView<MoleculeType, BenchFunctor> view2(cell2._particleSoABuffer, kSortingDirection, sortedSoa2,
                                                           projIdx2);

  std::vector<std::size_t> maxIndex, minIndex;
  const auto sd = cellFunctor.computeSortingData(projIdx1, projIdx2, maxIndex, minIndex);

  const std::size_t n1 = projIdx1.size();
  const std::size_t n2 = projIdx2.size();
  if (n1 == 0 or n2 == 0) {
    functor.endTraversal(true);
    return 0.0;
  }

  std::size_t survivingWindow = 0;
  for (std::size_t i = sd.startI; i < sd.endI; ++i) {
    survivingWindow += maxIndex[i] - minIndex[i];
  }

  functor.endTraversal(true);
  return static_cast<double>(survivingWindow) / static_cast<double>(n1 * n2);
}

}  // namespace

int main() {
  const std::vector<std::size_t> nValues = {50, 75, 100};
  const std::vector<double> hitrates = {0, 5, 10, 15, 20, 30, 50};
  constexpr int kRepsPerConfig = 20;

  std::printf("N,hitrate,mean_prune_ratio,min_prune_ratio,max_prune_ratio\n");
  unsigned seed = 42;
  for (auto n : nValues) {
    for (auto hitratePct : hitrates) {
      double sum = 0.0;
      double lo = 1.0;
      double hi = 0.0;
      for (int rep = 0; rep < kRepsPerConfig; ++rep) {
        const double ratio = pruningRatioSample(n, hitratePct / 100.0, seed++);
        sum += ratio;
        lo = std::min(lo, ratio);
        hi = std::max(hi, ratio);
      }
      std::printf("%zu,%.0f,%.4f,%.4f,%.4f\n", n, hitratePct, sum / kRepsPerConfig, lo, hi);
    }
  }
  return 0;
}
