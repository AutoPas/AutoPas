/**
 * @file PatternBenchmark.h
 * @author J. Rief
 * @date 29.08.25
 */

#pragma once

#include <filesystem>
#include <fstream>
#include <vector>

#include "autopas/cells/FullParticleCell.h"
#include "autopas/options/VectorizationPatternOption.h"
#include "autopas/utils/Timer.h"

namespace autopas {

/* PatternBenchmark class contains all functionalities of the benchmark vectorization pattern selection.
 * It calculates and store the benchmark results
 */
class PatternBenchmark {
 public:
  /**
   * maximum particles per cell for vectorization pattern selection benchmarking.
   */
  static constexpr size_t _benchmarkSize = 30;


  /**
   * Helper method to initialize particle cells used in pattern benchmark by adding particles to each cell to match
   * given numbers of particles and hit rates.
   *
   * @tparam Functor_T Functor type for benchmark.
   * @tparam Cell_T Cell type for benchmark.
   * @tparam Particle_T Particle type in benchmark.
   * @param functor functor that is used in benchmark.
   * @param cells list of cells included in benchmark
   * @param numParticlesPerCell list of number of particles per cell
   * @param cutoff cutoff used in benchmark.
   * @param hitRate average hitRate that a particle pair is inside cutoff range in benchmark.
   */
  template <class Functor_T, class Cell_T, class Particle_T>
  void initialization(Functor_T &functor, std::pair<Cell_T, Cell_T> &cells,
                      const std::pair<size_t, size_t> &numParticlesPerCell, double cutoff, double hitRate) {
    // initialize cells with randomly distributed particles

    // this is a formula determined by regression (based on a mapping from hitrate to div with random sample values)
    const double div = 2.86 * (hitRate * hitRate * hitRate) - 4.13 * (hitRate * hitRate) + 2.81 * hitRate + 0.42;
    const double cellLength = cutoff / div;

    cells.first.reserve(numParticlesPerCell.first);
    cells.second.reserve(numParticlesPerCell.second);
    for (size_t cellId = 0; cellId < 2; ++cellId) {
      for (size_t particleId = 0; particleId < (cellId == 0 ? numParticlesPerCell.first : numParticlesPerCell.second);
           ++particleId) {
        Particle_T p{
            {
                // particles are next to each other in X direction
                std::rand() / static_cast<double>(RAND_MAX) * cellLength + cellLength * static_cast<double>(cellId),
                std::rand() / static_cast<double>(RAND_MAX) * cellLength,
                std::rand() / static_cast<double>(RAND_MAX) * cellLength,
            },
            {
                0.,
                0.,
                0.,
            },
            // every cell gets its own id space
            particleId + ((std::numeric_limits<size_t>::max() / 2) * cellId),
            particleId % 5};
        (cellId == 0 ? cells.first : cells.second).addParticle(p);
      }
      functor.SoALoader((cellId == 0 ? cells.first : cells.second),
                        (cellId == 0 ? cells.first : cells.second)._particleSoABuffer, 0, false);
    }
  }

  /**
   * Helper method to apply functor in pattern benchmark.
   *
   * @tparam Functor_T Functor type for benchmark.
   * @tparam Cell_T Cell type for benchmark.
   * @param functor functor that is used in benchmark.
   * @param cells list of cells included in benchmark
   * @param timer timer to track used time in benchmark.
   * @param newton3 true if newton3 optimization is used.
   */
  template <class Functor_T, class Cell_T>
  void applyFunctorPair(Functor_T &functor, std::pair<Cell_T, Cell_T> &cells,
                        std::map<std::string, autopas::utils::Timer> &timer, bool newton3) {
    timer.at("Functor").start();

    functor.SoAFunctorPair(cells.first._particleSoABuffer, cells.second._particleSoABuffer, newton3);

    timer.at("Functor").stop();
  }
  /**
   * Helper method to run a single pattern benchmark for one combination of the number of particles in the first and
   * second SoA buffer.
   *
   * @tparam Functor_T Functor type for benchmark.
   * @tparam Particle_T Particle type in benchmark.
   * @param functor functor that is used in benchmark.
   * @param repetitions amount of repetitions for this combination; particle positions are newly generated in each
   * repetition
   * @param iterations amount of iterations per repetition; particle positions remain the same for all iterations; per
   * iteration the pairwise forces between particles in the first and second SoA buffer are calculated and the execution
   * time is measured
   * @param firstnumParticles amount of particles in first SoA buffer
   * @param secondnumParticles amount of particles in second SoA buffer
   * @param hitRate average hitRate that a particle pair is inside cutoff range in benchmark.
   * @param vecPattern vectorization pattern applied in this benchmark
   * @param timer timer to measure time in benchmark
   * @param newton3 true if newton3 optimization is applied
   * @return amount of time the benchmark took in nanoseconds
   */
  template <class Functor_T, class Particle_T>
  unsigned long runSingleBenchmark(Functor_T &functor, size_t repetitions, size_t iterations, size_t firstnumParticles,
                                   size_t secondnumParticles, double hitRate,
                                   autopas::VectorizationPatternOption::Value vecPattern,
                                   std::map<std::string, autopas::utils::Timer> &timer, bool newton3) {
    using ParticleType = Particle_T;
    using CellType = autopas::FullParticleCell<ParticleType>;

    const double cutoff{functor.getCutoff()};  // is also the cell size

    std::vector<long> times{};
    times.reserve(repetitions);
    functor.setVecPattern(vecPattern);

    for (int n = 0; n < repetitions; ++n) {
      // define scenario
      const std::pair<size_t, size_t> numParticlesPerCell{firstnumParticles, secondnumParticles};

      // repeat the whole experiment multiple times and average results
      std::pair<CellType, CellType> cells{};

      initialization<Functor_T, CellType, ParticleType>(functor, cells, numParticlesPerCell, cutoff, hitRate);
      for (size_t iteration = 0; iteration < iterations; ++iteration) {
        // actual benchmark
        applyFunctorPair(functor, cells, timer, newton3);
      }

      times.push_back(timer["Functor"].getTotalTime());
      for (auto &[_, t] : timer) {
        t.reset();
      }
    }
    unsigned long sum = std::accumulate(times.begin(), times.end(), static_cast<long>(0));
    return sum;
  }

  /**
   * This method calculates the benchmarkSize x benchmarkSize pattern benchmark to determine the optimal vector pattern
   * for some given combination of SoA sizes.
   * @tparam Functor_T Functor type for benchmark.
   * @tparam Particle_T Particle type in benchmark.
   * @param functor functor that is used in benchmark.
   * @param newton3 true if newton3 optimization is applied
   * @return array with optimal vector pattern for different combinations of SoA sizes.
   */
  template <class Functor_T, class Particle_T>
  std::array<autopas::VectorizationPatternOption::Value, _benchmarkSize * _benchmarkSize> calculateVecPatternMap(
      Functor_T &functor, bool newton3, double hitRate) {
    using PatternType = autopas::VectorizationPatternOption::Value;
    std::set<VectorizationPatternOption> patternSet = VectorizationPatternOption::getAllOptions();
    std::vector<PatternType> patterns{};
    for (auto vecOption : patternSet) {
      patterns.push_back(vecOption);
    }
    constexpr size_t numPatterns = 4;
    if (patterns.size() != numPatterns) {
      utils::ExceptionHandler::exception("PatternBenchmark assumes there are only {} vectorization patterns, but"
      "there is actually {}!", numPatterns, patterns.size());
    }
    std::array<std::array<unsigned long, 4>, _benchmarkSize * _benchmarkSize> allResults{};
    std::map<std::string, autopas::utils::Timer> timer{
        {"Initialization", autopas::utils::Timer()},
        {"Functor", autopas::utils::Timer()},
        {"Output", autopas::utils::Timer()},
        {"InteractionCounter", autopas::utils::Timer()},
    };
    for (const auto pattern : patterns) {

      for (size_t numberParticlesSoA1 = 1; numberParticlesSoA1 <= _benchmarkSize; numberParticlesSoA1++) {
        for (size_t numberParticlesSoA2 = 1; numberParticlesSoA2 <= _benchmarkSize; numberParticlesSoA2++) {
          // we set the hit rate to 16 percent, as this is the average for linked cells

          allResults[numberParticlesSoA1 - 1 + _benchmarkSize * (numberParticlesSoA2 - 1)][pattern] =
              runSingleBenchmark<Functor_T, Particle_T>(functor, 100, 100, numberParticlesSoA1, numberParticlesSoA2,
                                                        hitRate, pattern, timer, newton3);
        }
      }
    }

    std::array<PatternType, _benchmarkSize * _benchmarkSize> optimalPatterns{};
    for (size_t i = 0; i < _benchmarkSize * _benchmarkSize; i++) {
      long min = allResults[i][PatternType::p1xVec];
      PatternType minIndex = PatternType::p1xVec;
      for (const auto pattern : patterns) {
        if (allResults[i][pattern] < min) {
          min = allResults[i][pattern];
          minIndex = pattern;
        }
      }

      optimalPatterns[i] = minIndex;
    }

    return optimalPatterns;
  }

  /**
   * This method initializes the PatternBenchmark class by running the benchmark with and without newton3 applied and
   * storing the results
   * @tparam Functor_T Functor type for benchmark.
   * @tparam Particle_T Particle type in benchmark.
   * @param functor functor that is used in benchmark.
   * @param printPatternResults whether csv output for the benchmark results should be created
   */
  template <class Functor_T, typename Particle_T>
  void runBenchmark(Functor_T &functor, bool printPatternResults) {
    _optimalPatternsNewton3OnSide = calculateVecPatternMap<Functor_T, Particle_T>(functor, true, 0.3);
    _optimalPatternsNewton3OffSide = calculateVecPatternMap<Functor_T, Particle_T>(functor, false, 0.3);
    _optimalPatternsNewton3OnEdge = calculateVecPatternMap<Functor_T, Particle_T>(functor, true, 0.08);
    _optimalPatternsNewton3OffEdge = calculateVecPatternMap<Functor_T, Particle_T>(functor, false, 0.08);
    _optimalPatternsNewton3OnCorner = calculateVecPatternMap<Functor_T, Particle_T>(functor, true, 0.02);
    _optimalPatternsNewton3OffCorner = calculateVecPatternMap<Functor_T, Particle_T>(functor, false, 0.02);
    _patternsCalculated = true;

    if (printPatternResults) {
      const auto *fillerAfterSuffix = outputSuffix.empty() or outputSuffix.back() == '_' ? "" : "_";
      std::ofstream newton3SideFile(
          "newton3on_side_patterns" + outputSuffix + fillerAfterSuffix + utils::Timer::getDateStamp() + ".csv",
          std::ios::out);
      if (newton3SideFile.is_open()) {
        newton3SideFile << "pattern,firstCellSize,SecondCellSize\n";
        for (size_t secondSize = 1; secondSize <= _benchmarkSize; secondSize++) {
          for (size_t firstSize = 1; firstSize <= _benchmarkSize; firstSize++) {
            newton3SideFile << VectorizationPatternOption{_optimalPatternsNewton3OnSide[(secondSize - 1) * _benchmarkSize +
                                                                                (firstSize - 1)]}
                               .to_string()
                        << "," << firstSize << "," << secondSize << "\n";
          }
        }
        newton3SideFile.close();
      }
      std::ofstream newton3offSideFile(
          "newton3off_side_patterns" + outputSuffix + fillerAfterSuffix + utils::Timer::getDateStamp() + ".csv",
          std::ios::out);
      if (newton3offSideFile.is_open()) {
        newton3offSideFile << "pattern,firstCellSize,SecondCellSize\n";
        for (size_t secondSize = 1; secondSize <= _benchmarkSize; secondSize++) {
          for (size_t firstSize = 1; firstSize <= _benchmarkSize; firstSize++) {
            newton3offSideFile << VectorizationPatternOption{_optimalPatternsNewton3OffSide[(secondSize - 1) * _benchmarkSize +
                                                                                    (firstSize - 1)]}
                                  .to_string()
                           << "," << firstSize << "," << secondSize << "\n";
          }
        }
        newton3offSideFile.close();
      }
      std::ofstream newton3EdgeFile(
          "newton3on_edge_patterns" + outputSuffix + fillerAfterSuffix + utils::Timer::getDateStamp() + ".csv",
          std::ios::out);
      if (newton3EdgeFile.is_open()) {
        newton3EdgeFile << "pattern,firstCellSize,SecondCellSize\n";
        for (size_t secondSize = 1; secondSize <= _benchmarkSize; secondSize++) {
          for (size_t firstSize = 1; firstSize <= _benchmarkSize; firstSize++) {
            newton3EdgeFile << VectorizationPatternOption{_optimalPatternsNewton3OnSide[(secondSize - 1) * _benchmarkSize +
                                                                                (firstSize - 1)]}
            .to_string()
     << "," << firstSize << "," << secondSize << "\n";
          }
        }
        newton3EdgeFile.close();
      }
      std::ofstream newton3offEdgeFile(
          "newton3off_edge_patterns" + outputSuffix + fillerAfterSuffix + utils::Timer::getDateStamp() + ".csv",
          std::ios::out);
      if (newton3offEdgeFile.is_open()) {
        newton3offEdgeFile << "pattern,firstCellSize,SecondCellSize\n";
        for (size_t secondSize = 1; secondSize <= _benchmarkSize; secondSize++) {
          for (size_t firstSize = 1; firstSize <= _benchmarkSize; firstSize++) {
            newton3offEdgeFile << VectorizationPatternOption{_optimalPatternsNewton3OffSide[(secondSize - 1) * _benchmarkSize +
                                                                                    (firstSize - 1)]}
            .to_string()
     << "," << firstSize << "," << secondSize << "\n";
          }
        }
        newton3offEdgeFile.close();
      }
      std::ofstream newton3CornerFile(
          "newton3on_corner_patterns" + outputSuffix + fillerAfterSuffix + utils::Timer::getDateStamp() + ".csv",
          std::ios::out);
      if (newton3CornerFile.is_open()) {
        newton3CornerFile << "pattern,firstCellSize,SecondCellSize\n";
        for (size_t secondSize = 1; secondSize <= _benchmarkSize; secondSize++) {
          for (size_t firstSize = 1; firstSize <= _benchmarkSize; firstSize++) {
            newton3CornerFile << VectorizationPatternOption{_optimalPatternsNewton3OnSide[(secondSize - 1) * _benchmarkSize +
                                                                                (firstSize - 1)]}
            .to_string()
     << "," << firstSize << "," << secondSize << "\n";
          }
        }
        newton3EdgeFile.close();
      }
      std::ofstream newton3offCornerFile(
          "newton3off_corner_patterns" + outputSuffix + fillerAfterSuffix + utils::Timer::getDateStamp() + ".csv",
          std::ios::out);
      if (newton3offCornerFile.is_open()) {
        newton3offCornerFile << "pattern,firstCellSize,SecondCellSize\n";
        for (size_t secondSize = 1; secondSize <= _benchmarkSize; secondSize++) {
          for (size_t firstSize = 1; firstSize <= _benchmarkSize; firstSize++) {
            newton3offCornerFile << VectorizationPatternOption{_optimalPatternsNewton3OffSide[(secondSize - 1) * _benchmarkSize +
                                                                                    (firstSize - 1)]}
            .to_string()
     << "," << firstSize << "," << secondSize << "\n";
          }
        }
        newton3offCornerFile.close();
      }
    }
  }
  /**
   * This method gets an optimal vectorization pattern from the pattern benchmark results for a buffer-size-pair also
   * depending on the newton3 optimization. If at least one buffer has 0 particles it should default to p1xVec, as in
   * this case there the optimal vectorization pattern doesn't exist as there are no pairwise force calculations between
   * the particles in the two buffers. As we want a standardized behavior in this edge case, we always default to
   * p1xVec.
   * @param firstBufferSize number of particles in the first buffer
   * @param secondBufferSize number of particles in the second buffer
   * @param newton3 whether newton3 optimization is applied.
   */
  inline VectorizationPatternOption::Value getBenchmarkResult(size_t firstBufferSize, size_t secondBufferSize,
                                                              bool newton3, size_t type) { // type is just a temp name and should be an enum: 0 = side, 1 = edge, 2 =  corner
    autopas::VectorizationPatternOption::Value vectorizationPattern =
        autopas::VectorizationPatternOption::Value::p1xVec;
    // if at least one buffer has 0 particles it should just return p1xVec
    if (firstBufferSize != 0 and secondBufferSize != 0) {
      if (newton3) {
        if (type == 0) {
          vectorizationPattern =
              _optimalPatternsNewton3OnSide[(std::min(firstBufferSize, _benchmarkSize) - 1) +
                                          _benchmarkSize * (std::min(secondBufferSize, _benchmarkSize) - 1)];
        } else if (type == 1) {
          vectorizationPattern =
              _optimalPatternsNewton3OnEdge[(std::min(firstBufferSize, _benchmarkSize) - 1)];
        } else if (type == 2){
           vectorizationPattern =
             _optimalPatternsNewton3OnCorner[(std::min(firstBufferSize, _benchmarkSize) - 1)];
        }
      } else {
        if (type == 0) {
          vectorizationPattern =
              _optimalPatternsNewton3OffSide[(std::min(firstBufferSize, _benchmarkSize) - 1) +
                                           _benchmarkSize * (std::min(secondBufferSize, _benchmarkSize) - 1)];
        } else if (type == 1) {
          vectorizationPattern =
              _optimalPatternsNewton3OffEdge[(std::min(firstBufferSize, _benchmarkSize) - 1)];
        } else if (type == 2){
           vectorizationPattern =
             _optimalPatternsNewton3OffCorner[(std::min(firstBufferSize, _benchmarkSize) - 1)];
        }
      }
    }
    return vectorizationPattern;
  }

  /**
   * both variables represent the optimal pattern results for vectorization pattern selection benchmarking. They are
   * only set once at the beginning of the simulation.
   */
  std::array<autopas::VectorizationPatternOption::Value, _benchmarkSize * _benchmarkSize> _optimalPatternsNewton3OnSide;
  std::array<autopas::VectorizationPatternOption::Value, _benchmarkSize * _benchmarkSize> _optimalPatternsNewton3OffSide;
  std::array<autopas::VectorizationPatternOption::Value, _benchmarkSize * _benchmarkSize> _optimalPatternsNewton3OnEdge;
  std::array<autopas::VectorizationPatternOption::Value, _benchmarkSize * _benchmarkSize> _optimalPatternsNewton3OffEdge;
  std::array<autopas::VectorizationPatternOption::Value, _benchmarkSize * _benchmarkSize> _optimalPatternsNewton3OnCorner;
  std::array<autopas::VectorizationPatternOption::Value, _benchmarkSize * _benchmarkSize> _optimalPatternsNewton3OffCorner;

  /**
   * boolean to determine if pattern benchmark has already been executed
   */
  bool _patternsCalculated{false};
  // string which store the outputSuffix for the filenames of the pattern benchmark results
  std::string outputSuffix;
};
}  // namespace autopas