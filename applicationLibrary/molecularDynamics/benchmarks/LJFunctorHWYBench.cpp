/**
 * @file LJFunctorHWYBench.cpp
 * @date 25.6.2026
 * @author hmeyran
 *
 * Micro-benchmarks for mdLib::LJFunctorHWY covering three groups:
 *
 *   1. Face pair — hitrate-controlled (TwoCellsInteractionHitrateGenerator):
 *        a. Hitrate study: vecPattern fixed to p1xVec; hitrate and N (number of particles) swept.
 *        b. VecPattern study: all four vecPatterns × hitrate × N (reduced set of N to reduce number of benchmarks).
 *      Cells are face-adjacent (shared yz-plane at x = kBoundary).
 *
 *   2. Cell Pair with specific Layout,  uniform generator.
 *      Each layout has an unsorted (SoAFunctorPair) and a sorted
 *      (SoAFunctorPairSorted) variant. all four vecPatterns are swept.
 *
 *   3. Baselines: AoSFunctor, SoAFunctorSingle, SoAFunctorVerlet — single cell,
 *      varying N and newton3. Reference points for absolute throughput.
 *
 * Cell size: kCellSize = 5 (≈ 1.67 × cutoff), N ∈ {10…150} per cell.
 */

#include <benchmark/benchmark.h>

#include <array>
#include <vector>

#include "autopas/baseFunctors/CellFunctor.h"
#include "autopas/cells/FullParticleCell.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/utils/AlignedAllocator.h"
#include "autopas/utils/ArrayMath.h"
#include "autopasTools/generators/TwoCellsInteractionHitrateGenerator.h"
#include "autopasTools/generators/UniformGenerator.h"
#include "molecularDynamicsLibrary/LJFunctorHWY.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"

namespace {

using Molecule_Type = mdLib::MoleculeLJ;
using FMCell = autopas::FullParticleCell<Molecule_Type>;
using VectorizationPattern = autopas::VectorizationPatternOption::Value;
/**
 * Cutoff used by the Benchmarks
 */
constexpr double kCutoff = 3.0;
/**
 * Epsilon used by the Benchmarks.
 */
constexpr double kEpsilon = 1.0;
/**
 * Sigma used by the Benchmarks.
 */
constexpr double kSigma = 1.0;
/**
 * Lower Left Corner of the Cells used by the Benchmarks.
 */
constexpr std::array<double, 3> kLow{0.0, 0.0, 0.0};
/**
 * Upper Right Corner of the Cells used by the Benchmarks.
 */
constexpr std::array<double, 3> kHigh{5.0, 5.0, 5.0};
/**
 * X-coordinate of the shared cell boundary used by the Benchmarks.
 */
constexpr double kBoundary = 5.0;
/**
 * X-extent of each pair cell (> 1.5 * kCutoff = 4.5) used by the Benchmarks.
 */
constexpr double kCellSize = 5.0;

/**
 * Number of Particles swept by the Benchmarks in the Hitrate, Layout and Baseline Benchmarks.
 */
const std::vector<int64_t> kNValues = {10, 25, 50, 75, 100, 125, 150};
/**
 * Number of Particles swept by the Benchmarks in the VecPattern study, reduced to reduce number of different Benchmarks
 * run.
 */
const std::vector<int64_t> kNValuesReduced = {50, 75, 100};
/**
 * Hitrates swept by the Benchmarks in the Hitrate study.
 */
const std::vector<int64_t> hitrates = {0, 5, 10, 15, 20, 30, 50};

/**
 * VecPattern used by the Benchmarks in the Hitrate study.
 */
constexpr int kFixedVecPattern = VectorizationPattern::p1xVec;

/**
 * Functor config used to benchmark.
 */
using BenchFunctor = mdLib::LJFunctorHWY<Molecule_Type, /*shifting=*/true, /*useMixing=*/false,
                                         autopas::FunctorN3Modes::Both, /*calculateGlobals=*/false,
                                         /*countFLOPs=*/false>;
/**
 * Fills a Cell with Particles (of type MoleculeLJ) using the AutoPas uniform generator.
 * @param cell Cell to fill.
 * @param low lower left corner of the cell to fill.
 * @param high upper right corner of the cell to fill.
 * @param n number of particles to fill the cell with.
 * @param seed seed used by the Random Number Generator.
 */
void fillCell(FMCell &cell, const std::array<double, 3> &low, const std::array<double, 3> &high, std::size_t n,
              unsigned seed) {
  const Molecule_Type defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell, defaultParticle, low, high, n, seed);
}

/**
 * Builds neighborLists for use in the VerletList Benchmark
 * @param cell Cell to build NeighborLists for
 * @param interactionLen
 * @param newton3 wether newton3 optimization is being used or not
 * @return
 */
std::vector<std::vector<std::size_t, autopas::AlignedAllocator<std::size_t>>> buildNeighborLists(const FMCell &cell,
                                                                                                 double interactionLen,
                                                                                                 bool newton3) {
  using namespace autopas::utils::ArrayMath::literals;
  std::vector<std::vector<std::size_t, autopas::AlignedAllocator<std::size_t>>> lists(cell.size());
  const double interactionLenSq = interactionLen * interactionLen;
  for (std::size_t i = 0; i < cell.size(); ++i) {
    for (std::size_t j = newton3 ? i + 1 : 0; j < cell.size(); ++j) {
      if (i == j) continue;
      const auto dr = cell[i].getR() - cell[j].getR();
      if (autopas::utils::ArrayMath::dot(dr, dr) <= interactionLenSq) {
        lists[i].push_back(j);
      }
    }
  }
  return lists;
}

/**
 * Creates the Functor use for Benchmarking with cutoff, sigma and epsilon set.
 * @return the Functor.
 */
BenchFunctor makeFunctor() {
  BenchFunctor f(kCutoff);
  f.setParticleProperties(kEpsilon * 24.0, kSigma * kSigma);
  return f;
}

}  // namespace

/**
 * Benchmark of the AoSFunctor
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 *
 * The Benchmark uses fillCell() with kLow, kHigh and the seed 42.
 * It manually loops through all particles, applying the AoSFunctor for each needed pair (depending on newton3)
 * and measures the total calculation time.
 */
static void BM_AoSFunctor(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;

  FMCell cell;
  fillCell(cell, kLow, kHigh, n, 42);

  auto functor = makeFunctor();
  functor.initTraversal();

  for (auto _ : state) {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = newton3 ? i + 1 : 0; j < n; ++j) {
        if (i == j) continue;
        functor.AoSFunctor(cell[i], cell[j], newton3);
      }
    }
    benchmark::DoNotOptimize(cell);
  }
  functor.endTraversal(newton3);
}
/**
 * @name BM_AoS
 * @details
 * - Arguments: {Number of Particles, newton3}
 * - Sweeps kNValues with both newton3 on and off
 */
BENCHMARK(BM_AoSFunctor)->ArgsProduct({kNValues, {0, 1}})->ArgNames({"N", "n3"})->Repetitions(5)->Name("BM_AoS");

/**
 * Benchmark of the AoSFunctorPair (sorted path) with Cells in a Face Layout, via CellFunctor.
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 *
 * Uses the same cell layout and UniformGenerator as BM_SoAFunctorSortedPairFaceUniform so that
 * the two benchmarks are directly comparable. Sorting is forced by setting aosSortingThreshold=0 and
 * passing the face-normal direction (x-axis). No SoA load/extract is involved.
 */
static void BM_AoSFunctorSortedPairFaceUniform(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;

  FMCell cell1, cell2;
  const Molecule_Type defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell2, defaultParticle, {kBoundary, kLow[1], kLow[2]},
                                                                {kBoundary + kCellSize, kHigh[1], kHigh[2]}, n, 43);

  auto functor = makeFunctor();
  functor.initTraversal();

  autopas::internal::CellFunctor<FMCell, BenchFunctor, /*bidirectional=*/false> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::aos, newton3);
  cellFunctor.setAoSSortingThreshold(0);  // always take the sorted path

  const std::array<double, 3> sortingDirection{1.0, 0.0, 0.0};

  for (auto _ : state) {
    cellFunctor.processCellPair(cell1, cell2, sortingDirection);
    benchmark::DoNotOptimize(cell1);
    benchmark::DoNotOptimize(cell2);
  }
  functor.endTraversal(newton3);
}
/**
 * @name BM_AoS_PairSorted_Face_Uniform
 * @details
 * - Arguments: {Number of Particles, newton3}
 * - Sweeps kNValues with both newton3 on and off
 * - This is part of the Layout study; directly comparable to BM_SoA_PairSorted_Face_Uniform.
 */
BENCHMARK(BM_AoSFunctorSortedPairFaceUniform)
    ->ArgsProduct({kNValues, {0, 1}})
    ->ArgNames({"N", "n3"})
    ->Repetitions(5)
    ->Name("BM_AoS_PairSorted_Face_Uniform");

/**
 * Benchmark of the SoAFunctorSingle
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 *
 * The Benchmark uses fillCell() with kLow, kHigh and the seed 42.
 * It measures total time of the SoALoader + SoaFunctorSingle + SoAExtractor to get a full picture and capture memory
 * performance.
 */
static void BM_SoAFunctorSingle(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;

  FMCell cell;
  fillCell(cell, kLow, kHigh, n, 42);

  auto functor = makeFunctor();
  functor.initTraversal();

  for (auto _ : state) {
    functor.SoALoader(cell, cell._particleSoABuffer, 0, /*skipSoAResize=*/false);
    functor.SoAFunctorSingle(cell._particleSoABuffer, newton3);
    functor.SoAExtractor(cell, cell._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
/**
 * @name BM_SoA_Single
 * @details
 * - Arguments: {Number of Particles, newton3}
 * - Sweeps kNValues with both newton3 on and off
 */
BENCHMARK(BM_SoAFunctorSingle)
    ->ArgsProduct({kNValues, {0, 1}})
    ->ArgNames({"N", "n3"})
    ->Repetitions(5)
    ->Name("BM_SoA_Single");

/**
 * Benchmark of the SoAFunctorPair (unsorted path) with Cells in a Face Layout, via CellFunctor.
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 * - state.range(2): The VecPattern to use.
 * - state.range(3): The hitrate to use.
 *
 * The Benchmark uses the TwoCellsInteractionHitrateGenerator to generate the cells in a Face Layout with the specified
 * n and hitrate using the seed 42. It measures total time of the SoALoader + SoAFunctorPair + SoAExtractor (Loader
 * and Extractor for both SoAs) to get a full picture and capture memory performance.
 * Sorting is disabled by passing a zero sortingDirection to the CellFunctor.
 */
static void BM_SoAFunctorPairFace(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));
  const double hitrate = state.range(3) / 100.0;

  FMCell cell1, cell2;
  const Molecule_Type defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
      cell1, cell2, kLow, {kBoundary, kHigh[1], kHigh[2]}, {kBoundary, kLow[1], kLow[2]},
      {kBoundary + kCellSize, kHigh[1], kHigh[2]}, n, hitrate, kCutoff, defaultParticle, /*seed=*/42);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  autopas::internal::CellFunctor<FMCell, BenchFunctor, /*bidirectional=*/false> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::soa, newton3);

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    cellFunctor.processCellPair(cell1, cell2);  // {0,0,0} default direction disables sorting
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
/**
 * @name BM_SoA_Pair_Face_Hitrate
 * @details
 * - Arguments: {Number of Particles, newton3, VecPattern, hitrate}
 * - Sweeps kNValues, both newton3 on and off, kFixedVecPattern with all hitrates defined in the global hitrates
 * - This constitutes the Hitrate study.
 */
BENCHMARK(BM_SoAFunctorPairFace)
    ->ArgsProduct({kNValues, {0, 1}, {kFixedVecPattern}, hitrates})
    ->ArgNames({"N", "n3", "vecPat", "hitrate%"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Face_Hitrate");
/**
 * @name BM_SoA_Pair_Face_VecPatterns
 * @details
 * - Arguments: {Number of Particles, newton3, VecPattern, hitrate}
 * - Sweeps kNValuesReduced, both newton3 on and off, all VecPatterns with all hitrates defined in the global hitrates
 * - This constitutes the VecPattern study.
 */
BENCHMARK(BM_SoAFunctorPairFace)
    ->ArgsProduct({kNValuesReduced,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2), static_cast<int>(VectorizationPattern::pVecx1)},
                   hitrates})
    ->ArgNames({"N", "n3", "vecPat", "hitrate%"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Face_VecPatterns");

/**
 * Benchmark of the SoAFunctorPairSorted (sorted path) with Cells in a Face Layout, via CellFunctor.
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 * - state.range(2): The VecPattern to use.
 * - state.range(3): The hitrate to use.
 *
 * The Benchmark uses the TwoCellsInteractionHitrateGenerator to generate the cells in a Face Layout with the specified
 * n and hitrate using the seed 42. It measures total time of the SoALoader + SoAFunctorPairSorted + SoAExtractor
 * (Loader and Extractor for both SoAs) to get a full picture and capture memory performance.
 * Sorting is forced by setting soaSortingThreshold=0 and passing the face-normal direction (x-axis).
 */
static void BM_SoAFunctorPairSortedFace(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));
  const double hitrate = state.range(3) / 100.0;

  FMCell cell1, cell2;
  const Molecule_Type defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::TwoCellsInteractionHitrateGenerator::fillWithParticles(
      cell1, cell2, kLow, {kBoundary, kHigh[1], kHigh[2]}, {kBoundary, kLow[1], kLow[2]},
      {kBoundary + kCellSize, kHigh[1], kHigh[2]}, n, hitrate, kCutoff, defaultParticle, /*seed=*/42);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  autopas::internal::CellFunctor<FMCell, BenchFunctor, /*bidirectional=*/false> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::soa, newton3);
  cellFunctor.setSoASortingThreshold(0);  // always take the sorted path

  const std::array<double, 3> sortingDirection{1.0, 0.0, 0.0};

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    cellFunctor.processCellPair(cell1, cell2, sortingDirection);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
/**
 * @name BM_SoA_PairSorted_Face_Hitrate
 * @details
 * - Arguments: {Number of Particles, newton3, VecPattern, hitrate}
 * - Sweeps kNValues, both newton3 on and off, kFixedVecPattern with all hitrates defined in the global hitrates
 * - This constitutes the Hitrate study.
 */
BENCHMARK(BM_SoAFunctorPairSortedFace)
    ->ArgsProduct({kNValues, {0, 1}, {kFixedVecPattern}, hitrates})
    ->ArgNames({"N", "n3", "vecPat", "hitrate%"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Face_Hitrate");
/**
 * @name BM_SoA_PairSorted_Face_VecPatterns
 * @details
 * - Arguments: {Number of Particles, newton3, VecPattern, hitrate}
 * - Sweeps kNValuesReduced, both newton3 on and off, all VecPatterns with all hitrates defined in the global hitrates
 * - This constitutes the VecPattern study.
 */
BENCHMARK(BM_SoAFunctorPairSortedFace)
    ->ArgsProduct({kNValuesReduced,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2), static_cast<int>(VectorizationPattern::pVecx1)},
                   hitrates})
    ->ArgNames({"N", "n3", "vecPat", "hitrate%"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Face_VecPatterns");

/**
 * Benchmark of the SoAFunctorPair (unsorted path) with Cells in a Face Layout, via CellFunctor.
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 * - state.range(2): The VecPattern to use.
 *
 * The Benchmark uses the UniformGenerator to generate the cells in a Face Layout, without a specified hitrate.
 * It measures total time of the SoALoader + SoAFunctorPair + SoAExtractor (Loader
 * and Extractor for both SoAs) to get a full picture and capture memory performance.
 * Sorting is disabled by passing a zero sortingDirection to the CellFunctor.
 */
static void BM_SoAFunctorPairFaceUniform(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule_Type defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell2, defaultParticle, {kBoundary, kLow[1], kLow[2]},
                                                                {kBoundary + kCellSize, kHigh[1], kHigh[2]}, n, 43);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  autopas::internal::CellFunctor<FMCell, BenchFunctor, /*bidirectional=*/false> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::soa, newton3);

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    cellFunctor.processCellPair(cell1, cell2);  // {0,0,0} default direction disables sorting
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
/**
 * @name BM_SoA_Pair_Face_Uniform
 * @details
 * - Arguments: {Number of Particles, newton3, VecPattern}
 * - Sweeps kNValues, both newton3 on and off with all VecPatterns.
 * - This is part of the Layout study.
 */
BENCHMARK(BM_SoAFunctorPairFaceUniform)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Face_Uniform");

/**
 * Benchmark of the SoAFunctorPairSorted (sorted path) with Cells in a Face Layout, via CellFunctor.
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 * - state.range(2): The VecPattern to use.
 *
 * The Benchmark uses the UniformGenerator to generate the cells in a Face Layout, without a specified hitrate.
 * It measures total time of the SoALoader + SoAFunctorPairSorted + SoAExtractor (Loader
 * and Extractor for both SoAs) to get a full picture and capture memory performance.
 * Sorting is forced by setting soaSortingThreshold=0 and passing the face-normal direction (x-axis).
 */
static void BM_SoAFunctorSortedPairFaceUniform(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule_Type defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell2, defaultParticle, {kBoundary, kLow[1], kLow[2]},
                                                                {kBoundary + kCellSize, kHigh[1], kHigh[2]}, n, 43);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  autopas::internal::CellFunctor<FMCell, BenchFunctor, /*bidirectional=*/false> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::soa, newton3);
  cellFunctor.setSoASortingThreshold(0);  // always take the sorted path

  const std::array<double, 3> sortingDirection{1.0, 0.0, 0.0};

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    cellFunctor.processCellPair(cell1, cell2, sortingDirection);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
/**
 * @name BM_SoA_PairSorted_Face_Uniform
 * @details
 * - Arguments: {Number of Particles, newton3, VecPattern}
 * - Sweeps kNValues, both newton3 on and off with all VecPatterns.
 * - This is part of the Layout study.
 */
BENCHMARK(BM_SoAFunctorSortedPairFaceUniform)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Face_Uniform");

/**
 * Benchmark of the SoAFunctorPair (unsorted path) with Cells in an Edge Layout, via CellFunctor.
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 * - state.range(2): The VecPattern to use.
 *
 * The Benchmark uses the UniformGenerator to generate the cells in an Edge Layout, without a specified hitrate.
 * It measures total time of the SoALoader + SoAFunctorPair + SoAExtractor (Loader
 * and Extractor for both SoAs) to get a full picture and capture memory performance.
 * Sorting is disabled by passing a zero sortingDirection to the CellFunctor.
 */
static void BM_SoAFunctorPairEdge(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule_Type defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2, defaultParticle, {kLow[0] + kCellSize, kLow[1] + kCellSize, kLow[2]},
      {kHigh[0] + kCellSize, kHigh[1] + kCellSize, kHigh[2]}, n, 43);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  autopas::internal::CellFunctor<FMCell, BenchFunctor, /*bidirectional=*/false> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::soa, newton3);

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    cellFunctor.processCellPair(cell1, cell2);  // {0,0,0} default direction disables sorting
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}

/**
 * @name BM_SoA_Pair_Edge
 * @details
 * - Arguments: {Number of Particles, newton3, VecPattern}
 * - Sweeps kNValues, both newton3 on and off with all VecPatterns.
 * - This is part of the Layout study.
 */
BENCHMARK(BM_SoAFunctorPairEdge)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Edge");

/**
 * Benchmark of the SoAFunctorPair (unsorted path) with Cells in a Corner Layout, via CellFunctor.
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 * - state.range(2): The VecPattern to use.
 *
 * The Benchmark uses the UniformGenerator to generate the cells in a Corner Layout, without a specified hitrate.
 * It measures total time of the SoALoader + SoAFunctorPair + SoAExtractor (Loader
 * and Extractor for both SoAs) to get a full picture and capture memory performance.
 * Sorting is disabled by passing a zero sortingDirection to the CellFunctor.
 */
static void BM_SoAFunctorPairCorner(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule_Type defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2, defaultParticle, {kLow[0] + kCellSize, kLow[1] + kCellSize, kLow[2] + kCellSize},
      {kHigh[0] + kCellSize, kHigh[1] + kCellSize, kHigh[2] + kCellSize}, n, 43);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  autopas::internal::CellFunctor<FMCell, BenchFunctor, /*bidirectional=*/false> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::soa, newton3);

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    cellFunctor.processCellPair(cell1, cell2);  // {0,0,0} default direction disables sorting
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}

/**
 * @name BM_SoA_Pair_Corner
 * @details
 * - Arguments: {Number of Particles, newton3, VecPattern}
 * - Sweeps kNValues, both newton3 on and off with all VecPatterns.
 * - This is part of the Layout study.
 */
BENCHMARK(BM_SoAFunctorPairCorner)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_Pair_Corner");

/**
 * Benchmark of the SoAFunctorPairSorted (sorted path) with Cells in an Edge Layout, via CellFunctor.
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 * - state.range(2): The VecPattern to use.
 *
 * The Benchmark uses the UniformGenerator to generate the cells in an Edge Layout, without a specified hitrate.
 * It measures total time of the SoALoader + SoAFunctorPairSorted + SoAExtractor (Loader
 * and Extractor for both SoAs) to get a full picture and capture memory performance.
 * Sorting is forced by setting soaSortingThreshold=0 and passing the edge diagonal direction (x+y axis).
 */
static void BM_SoAFunctorSortedPairEdge(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule_Type defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2, defaultParticle, {kLow[0] + kCellSize, kLow[1] + kCellSize, kLow[2]},
      {kHigh[0] + kCellSize, kHigh[1] + kCellSize, kHigh[2]}, n, 43);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  autopas::internal::CellFunctor<FMCell, BenchFunctor, /*bidirectional=*/false> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::soa, newton3);
  cellFunctor.setSoASortingThreshold(0);  // always take the sorted path

  const double normalized = 1.0 / sqrt(2.0);
  const std::array<double, 3> sortingDirection{normalized, normalized, 0.0};

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    cellFunctor.processCellPair(cell1, cell2, sortingDirection);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}

/**
 * @name BM_SoA_PairSorted_Edge
 * @details
 * - Arguments: {Number of Particles, newton3, VecPattern}
 * - Sweeps kNValues, both newton3 on and off with all VecPatterns.
 * - This is part of the Layout study.
 */
BENCHMARK(BM_SoAFunctorSortedPairEdge)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Edge");

/**
 * Benchmark of the SoAFunctorPairSorted (sorted path) with Cells in a Corner Layout, via CellFunctor.
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 * - state.range(2): The VecPattern to use.
 *
 * The Benchmark uses the UniformGenerator to generate the cells in a Corner Layout, without a specified hitrate.
 * It measures total time of the SoALoader + SoAFunctorPairSorted + SoAExtractor (Loader
 * and Extractor for both SoAs) to get a full picture and capture memory performance.
 * Sorting is forced by setting soaSortingThreshold=0 and passing the corner diagonal direction (x+y+z axis).
 */
static void BM_SoAFunctorSortedPairCorner(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;
  const auto pattern = static_cast<VectorizationPattern>(state.range(2));

  FMCell cell1, cell2;
  const Molecule_Type defaultParticle({0, 0, 0}, {0, 0, 0}, 0, 0);
  autopasTools::generators::UniformGenerator::fillWithParticles(cell1, defaultParticle, kLow, kHigh, n, 42);
  autopasTools::generators::UniformGenerator::fillWithParticles(
      cell2, defaultParticle, {kLow[0] + kCellSize, kLow[1] + kCellSize, kLow[2] + kCellSize},
      {kHigh[0] + kCellSize, kHigh[1] + kCellSize, kHigh[2] + kCellSize}, n, 43);

  auto functor = makeFunctor();
  functor.setVecPattern(pattern);
  functor.initTraversal();

  autopas::internal::CellFunctor<FMCell, BenchFunctor, /*bidirectional=*/false> cellFunctor(
      functor, kCutoff, autopas::DataLayoutOption::soa, newton3);
  cellFunctor.setSoASortingThreshold(0);  // always take the sorted path

  const double invSqrt3 = 1.0 / sqrt(3.0);
  const std::array<double, 3> sortingDirection{invSqrt3, invSqrt3, invSqrt3};

  for (auto _ : state) {
    functor.SoALoader(cell1, cell1._particleSoABuffer, 0, false);
    functor.SoALoader(cell2, cell2._particleSoABuffer, 0, false);
    cellFunctor.processCellPair(cell1, cell2, sortingDirection);
    functor.SoAExtractor(cell1, cell1._particleSoABuffer, 0);
    functor.SoAExtractor(cell2, cell2._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell1._particleSoABuffer);
    benchmark::DoNotOptimize(cell2._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}

/**
 * @name BM_SoA_PairSorted_Corner
 * @details
 * - Arguments: {Number of Particles, newton3, VecPattern}
 * - Sweeps kNValues, both newton3 on and off with all VecPatterns.
 * - This is part of the Layout study.
 */
BENCHMARK(BM_SoAFunctorSortedPairCorner)
    ->ArgsProduct({kNValues,
                   {0, 1},
                   {static_cast<int>(VectorizationPattern::p1xVec), static_cast<int>(VectorizationPattern::p2xVecDiv2),
                    static_cast<int>(VectorizationPattern::pVecDiv2x2),
                    static_cast<int>(VectorizationPattern::pVecx1)}})
    ->ArgNames({"N", "n3", "vecPat"})
    ->Repetitions(5)
    ->Name("BM_SoA_PairSorted_Corner");

/**
 * Benchmark of the SoAFunctorVerlet
 * @param state The Benchmark state.
 * - state.range(0): The number of Particles in the Cell
 * - state.range(1): Wether to use newton3 optimization or not.
 *
 * The Benchmark uses fillCell() with kLow, kHigh and the seed 42.
 * It measures total time of the SoALoader + SoAFunctorVerlet + SoAExtractor to get a full picture and capture memory
 * performance. The NeighborList generation is NOT timed.
 */
static void BM_SoAFunctorVerlet(benchmark::State &state) {
  const auto n = static_cast<std::size_t>(state.range(0));
  const bool newton3 = state.range(1) != 0;

  FMCell cell;
  fillCell(cell, kLow, kHigh, n, 42);

  const double skin = 0.5;
  const auto lists = buildNeighborLists(cell, kCutoff + skin, newton3);

  auto functor = makeFunctor();
  functor.initTraversal();

  for (auto _ : state) {
    functor.SoALoader(cell, cell._particleSoABuffer, 0, false);
    for (std::size_t i = 0; i < n; ++i) {
      functor.SoAFunctorVerlet(cell._particleSoABuffer, i, lists[i], newton3);
    }
    functor.SoAExtractor(cell, cell._particleSoABuffer, 0);
    benchmark::DoNotOptimize(cell._particleSoABuffer);
  }
  functor.endTraversal(newton3);
}
/**
 * @name BM_SoA_Verlet
 * @details
 * - Arguments: {Number of Particles, newton3}
 * - Sweeps kNValues with both newton3 on and off
 * - This is a Baseline Benchmark.
 */
BENCHMARK(BM_SoAFunctorVerlet)
    ->ArgsProduct({kNValues, {0, 1}})
    ->ArgNames({"N", "n3"})
    ->Repetitions(5)
    ->Name("BM_SoA_Verlet");

BENCHMARK_MAIN();
