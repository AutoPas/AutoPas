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
   * Helper method to run benchmark for optimal patterns selection
   *
   * @tparam Functor Functor type for benchmark.
   * @tparam Cell Cell type for benchmark.
   * @param functor functor that is used in benchmark.
   * @param cells list of cells included in benchmark
   */
  template <class Functor, class Cell>
  void csvOutput(Functor &functor, std::vector<Cell> &cells) {
    std::ofstream csvFile("particles.csv");
    if (not csvFile.is_open()) {
      throw std::runtime_error("FILE NOT OPEN!");
    }
    csvFile << "CellId,ParticleId,rX,rY,rZ,fX,fY,fZ\n";
    for (size_t cellId = 0; cellId < cells.size(); ++cellId) {
      functor.SoAExtractor(cells[cellId], cells[cellId]._particleSoABuffer, 0);
      for (size_t particleId = 0; particleId < cells[cellId].getNumberOfParticles(autopas::IteratorBehavior::owned);
           ++particleId) {
        const auto &p = cells[cellId][particleId];
        using autopas::utils::ArrayUtils::to_string;
        csvFile << cellId << "," << p.getID() << "," << to_string(p.getR(), ",", {"", ""}) << ","
                << to_string(p.getF(), ",", {"", ""}) << "\n";
      }
    }
    csvFile.close();
  }
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
  void initialization(Functor_T &functor, std::vector<Cell_T> &cells, const std::vector<size_t> &numParticlesPerCell,
                      double cutoff, double hitRate) {
    // initialize cells with randomly distributed particles

    // this is a formula determined by regression (based on a mapping from hitrate to div with random sample values)
    const double div = 2.86 * (hitRate * hitRate * hitRate) - 4.13 * (hitRate * hitRate) + 2.81 * hitRate + 0.42;
    const double cellLength = cutoff / div;

    cells[0].reserve(numParticlesPerCell[0]);
    cells[1].reserve(numParticlesPerCell[1]);
    for (size_t cellId = 0; cellId < numParticlesPerCell.size(); ++cellId) {
      for (size_t particleId = 0; particleId < numParticlesPerCell[cellId]; ++particleId) {
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
            particleId + ((std::numeric_limits<size_t>::max() / numParticlesPerCell.size()) * cellId),
            particleId % 5};
        cells[cellId].addParticle(p);
      }
      functor.SoALoader(cells[cellId], cells[cellId]._particleSoABuffer, 0, false);
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
  void applyFunctorPair(Functor_T &functor, std::vector<Cell_T> &cells,
                        std::map<std::string, autopas::utils::Timer> &timer, bool newton3) {
    timer.at("Functor").start();

    functor.SoAFunctorPair(cells[0]._particleSoABuffer, cells[1]._particleSoABuffer, newton3);

    timer.at("Functor").stop();
  }
  /**
   * Helper method to run a single pattern benchmark for one combination of the number of particles in the first and
   * second SoA buffer.
   *
   * @tparam Functor_T Functor type for benchmark.
   * @tparam Particle_T Particle type in benchmark.
   * @param functor functor that is used in benchmark.
   * @param repetitions amount of repetitions for this constellation; particle positions are newly generated in each
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
      const std::vector<size_t> numParticlesPerCell{firstnumParticles, secondnumParticles};

      // repeat the whole experiment multiple times and average results
      std::vector<CellType> cells{2};

      initialization<Functor_T, CellType, ParticleType>(functor, cells, numParticlesPerCell, cutoff, hitRate);
      for (size_t iteration = 0; iteration < iterations; ++iteration) {
        // actual benchmark
        applyFunctorPair(functor, cells, timer, newton3);
      }
      // print particles to CSV for checking and prevent compiler from optimizing everything away.
      csvOutput(functor, cells);
      times.push_back(timer["Functor"].getTotalTime());
      for (auto &[_, t] : timer) {
        t.reset();
      }
    }
    unsigned long sum = std::accumulate(times.begin(), times.end(), static_cast<long>(0));
    return sum;
  }
  /**
   * helper method to turn vectorization pattern enum into string
   *
   * @param vecPattern vectorization pattern
   * @return string
   */
  inline std::string checkVecPattern(const autopas::VectorizationPatternOption::Value vecPattern) {
    if (vecPattern == autopas::VectorizationPatternOption::Value::p1xVec) {
      return "p1xVec";
    } else if (vecPattern == autopas::VectorizationPatternOption::Value::p2xVecDiv2) {
      return "p2xVecDiv2";
    } else if (vecPattern == autopas::VectorizationPatternOption::Value::pVecDiv2x2) {
      return "pVecDiv2x2";
    } else {
      return "pVecx1";
    }
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
      Functor_T &functor, bool newton3) {
    using patterntype = autopas::VectorizationPatternOption::Value;
    std::vector<patterntype> patterns = {patterntype::p1xVec, patterntype::p2xVecDiv2, patterntype::pVecDiv2x2,
                                         patterntype::pVecx1};
    std::array<std::array<unsigned long, 4>, _benchmarkSize * _benchmarkSize> all_results{};
    std::map<std::string, autopas::utils::Timer> timer{
        {"Initialization", autopas::utils::Timer()},
        {"Functor", autopas::utils::Timer()},
        {"Output", autopas::utils::Timer()},
        {"InteractionCounter", autopas::utils::Timer()},
    };
    for (size_t i = 0; i < patterns.size(); i++) {
      patterntype current_pattern = patterns[i];
      for (size_t numberParticlesSoA1 = 1; numberParticlesSoA1 <= _benchmarkSize; numberParticlesSoA1++) {
        for (size_t numberParticlesSoA2 = 1; numberParticlesSoA2 <= _benchmarkSize; numberParticlesSoA2++) {
          // we set the hit rate to 16 percent, as this is the average for linked cells

          all_results[numberParticlesSoA1 - 1 + _benchmarkSize * (numberParticlesSoA2 - 1)][i] =
              runSingleBenchmark<Functor_T, Particle_T>(functor, 100, 100, numberParticlesSoA1, numberParticlesSoA2,
                                                        0.16, current_pattern, timer, newton3);
        }
      }
    }

    std::array<patterntype, _benchmarkSize * _benchmarkSize> optimalPatterns{};
    for (size_t i = 0; i < _benchmarkSize * _benchmarkSize; i++) {
      long min = all_results[i][0];
      long min_index = 0;
      for (int j = 1; j < patterns.size(); j++) {
        if (all_results[i][j] < min) {
          min = all_results[i][j];
          min_index = j;
        }
      }

      optimalPatterns[i] = patterns[min_index];
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
  void runBenchmark(Functor_T functor, bool printPatternResults) {
    _optimalPatternsNewton3On = calculateVecPatternMap<Functor_T, Particle_T>(functor, true);
    _optimalPatternsNewton3Off = calculateVecPatternMap<Functor_T, Particle_T>(functor, false);
    _patternsCalculated = true;

    if (printPatternResults) {
      const auto *fillerAfterSuffix = outputSuffix.empty() or outputSuffix.back() == '_' ? "" : "_";
      std::ofstream newton3File(
          "newton3on_patterns" + outputSuffix + fillerAfterSuffix + utils::Timer::getDateStamp() + ".csv",
          std::ios::out);
      if (newton3File.is_open()) {
        newton3File << "pattern,fcs,scs\n";
        for (size_t second_size = 1; second_size <= _benchmarkSize; second_size++) {
          for (size_t first_size = 1; first_size <= _benchmarkSize; first_size++) {
            newton3File << checkVecPattern(
                               _optimalPatternsNewton3On[(second_size - 1) * _benchmarkSize + (first_size - 1)])
                        << "," << first_size << "," << second_size << "\n";
          }
        }
        newton3File.close();
      }
      std::ofstream newton3offFile(
          "newton3off_patterns" + outputSuffix + fillerAfterSuffix + utils::Timer::getDateStamp() + ".csv",
          std::ios::out);
      if (newton3offFile.is_open()) {
        newton3offFile << "pattern,fcs,scs\n";
        for (size_t second_size = 1; second_size <= _benchmarkSize; second_size++) {
          for (size_t first_size = 1; first_size <= _benchmarkSize; first_size++) {
            newton3offFile << checkVecPattern(
                                  _optimalPatternsNewton3Off[(second_size - 1) * _benchmarkSize + (first_size - 1)])
                           << "," << first_size << "," << second_size << "\n";
          }
        }
        newton3offFile.close();
      }
    }
  }
  /**
   * This method initializes the PatternBenchmark class by running the benchmark with and without newton3 applied and
   * storing the results.
   * @param firstBufferSize number of particles in the first buffer
   * @param secondBufferSize number of particles in the second buffer
   * @param newton3 whether newton3 optimization is applied.
   */
  inline VectorizationPatternOption::Value getBenchmarkResult(size_t firstBufferSize, size_t secondBufferSize,
                                                              bool newton3) {
    autopas::VectorizationPatternOption::Value vectorizationPattern =
        autopas::VectorizationPatternOption::Value::p1xVec;
    if (firstBufferSize != 0 and secondBufferSize != 0) {
      if (newton3) {
        vectorizationPattern =
            _optimalPatternsNewton3On[(std::min(firstBufferSize, _benchmarkSize) - 1) +
                                      _benchmarkSize * (std::min(secondBufferSize, _benchmarkSize) - 1)];
      } else {
        vectorizationPattern =
            _optimalPatternsNewton3Off[(std::min(firstBufferSize, _benchmarkSize) - 1) +
                                       _benchmarkSize * (std::min(secondBufferSize, _benchmarkSize) - 1)];
      }
    }
    return vectorizationPattern;
  }

  /**
   * both variables represent the optimal pattern results for vectorization pattern selection benchmarking. They are
   * only set once at the beginning of the simulation.
   */
  std::array<autopas::VectorizationPatternOption::Value, _benchmarkSize * _benchmarkSize> _optimalPatternsNewton3On;
  std::array<autopas::VectorizationPatternOption::Value, _benchmarkSize * _benchmarkSize> _optimalPatternsNewton3Off;

  /**
   * boolean to determine if pattern benchmark has already been executed
   */
  bool _patternsCalculated{false};
  // string which store the outputSuffix for the filenames of the pattern benchmark results
  std::string outputSuffix;
};
}  // namespace autopas