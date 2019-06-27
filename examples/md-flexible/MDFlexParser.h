/**
 * @file MDFlexParser.h
 * @date 23.02.2018
 * @author F. Gratl
 */

#pragma once
#include "MDFlexParserInclude.h"


using namespace std;

class MDFlexParser {
 public:
  enum FunctorOption { lj12_6, lj12_6_AVX };
  enum GeneratorOption { grid, uniform, gaussian };

  MDFlexParser() = default;

  double getBoxLength();
  std::set<autopas::ContainerOption> getContainerOptions() const;
  autopas::SelectorStrategyOption getSelectorStrategy() const;
  double getCutoff() const;
  const autopas::NumberSet<double> &getCellSizeFactors() const;
  std::set<autopas::DataLayoutOption> getDataLayoutOptions() const;
  double getDistributionMean() const;
  double getDistributionStdDev() const;
  FunctorOption getFunctorOption() const;
  GeneratorOption getGeneratorOption() const;
  size_t getIterations() const;
  bool getMeasureFlops() const;
  std::set<autopas::Newton3Option> getNewton3Options() const;
  const std::string &getLogFileName() const;
  spdlog::level::level_enum getLogLevel() const;
  double getParticleSpacing() const;
  size_t getParticlesTotal() const;
  size_t getParticlesPerDim() const;
  double getEpsilon() const;
  double getSigma() const;
  double getDeltaT() const;

  unsigned int getTuningInterval() const;
  unsigned int getTuningSamples() const;
  autopas::TuningStrategyOption getTuningStrategyOption() const;
  std::string getWriteVTK() const;
  const std::set<autopas::TraversalOption> &getTraversalOptions() const;
  unsigned int getVerletRebuildFrequency() const;
  double getVerletSkinRadius() const;
  bool parseInput(int argc, char **argv);
  void printConfig();

  double getMass() const;

 private:
  static constexpr size_t valueOffset = 32;

  // defaults:
  // AutoPas options:
  std::set<autopas::ContainerOption> containerOptions = autopas::allContainerOptions;
  std::set<autopas::DataLayoutOption> dataLayoutOptions = autopas::allDataLayoutOptions;
  autopas::SelectorStrategyOption selectorStrategy = autopas::SelectorStrategyOption::fastestAbs;
  std::set<autopas::TraversalOption> traversalOptions = autopas::allTraversalOptions;
  autopas::TuningStrategyOption tuningStrategyOption = autopas::TuningStrategyOption::fullSearch;
  std::set<autopas::Newton3Option> newton3Options = autopas::allNewton3Options;
  std::unique_ptr<autopas::NumberSet<double>> cellSizeFactors =
      std::make_unique<autopas::NumberSetFinite<double>>(std::set<double>{1.});

  // Simulation Options:
  double boxLength = -1;
  double cutoff = 1.;
  double distributionMean = 5.;
  double distributionStdDev = 2.;
  FunctorOption functorOption = FunctorOption::lj12_6;
  GeneratorOption generatorOption = GeneratorOption::grid;
  size_t iterations = 10;
  spdlog::level::level_enum logLevel = spdlog::level::info;
  bool measureFlops = true;
  size_t particlesPerDim = 20;
  size_t particlesTotal = 1000;
  double particleSpacing = .4;
  unsigned int tuningInterval = 100;
  unsigned int tuningSamples = 3;
  std::string writeVTK = "";
  std::string logFileName = "";
  unsigned int verletRebuildFrequency = 5;
  double verletSkinRadius = .2;
  double epsilon = 5.0;
  double sigma = 1.0;
  double delta_t = 0.001;
  double mass = 1.0;
};
