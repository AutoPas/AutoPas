/**
 * @file MDFlexParser.h
 * @date 23.02.2018
 * @author F. Gratl
 */

#pragma once

#include <getopt.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include "autopas/AutoPas.h"

using namespace std;

class MDFlexParser {
 public:
  enum FunctorOption { lj12_6, lj12_6_AVX, lj12_6_Globals };
  enum GeneratorOption { grid, uniform, gaussian };
  enum PrecisionOption { FP32, FP64 };

  MDFlexParser() = default;

  double getBoxLength();
  std::vector<autopas::ContainerOption> getContainerOptions() const;
  autopas::SelectorStrategy getSelectorStrategy() const;
  double getCutoff() const;
  double getCellSizeFactor() const;
  vector<autopas::DataLayoutOption> getDataLayoutOptions() const;
  double getDistributionMean() const;
  double getDistributionStdDev() const;
  FunctorOption getFunctorOption() const;
  GeneratorOption getGeneratorOption() const;
  PrecisionOption getPrecisionOption() const;
  size_t getIterations() const;
  bool getMeasureFlops() const;
  std::vector<autopas::Newton3Option> getNewton3Options() const;
  const string &getLogFileName() const;
  spdlog::level::level_enum getLogLevel() const;
  double getParticleSpacing() const;
  size_t getParticlesTotal() const;
  size_t getParticlesPerDim() const;
  unsigned int getTuningInterval() const;
  unsigned int getTuningSamples() const;
  string getWriteVTK() const;
  const vector<autopas::TraversalOption> &getTraversalOptions() const;
  unsigned int getVerletRebuildFrequency() const;
  unsigned int getVerletClusterSize() const;
  double getVerletSkinRadius() const;
  bool parseInput(int argc, char **argv);
  void printConfig();

 private:
  static constexpr size_t valueOffset = 32;

  // defaults:
  std::vector<autopas::ContainerOption> containerOptions = autopas::allContainerOptions;
  autopas::SelectorStrategy selectorStrategy = autopas::SelectorStrategy::fastestAbs;
  std::vector<autopas::DataLayoutOption> dataLayoutOptions = autopas::allDataLayoutOptions;
  std::vector<autopas::TraversalOption> traversalOptions = autopas::allTraversalOptions;
  std::vector<autopas::Newton3Option> newton3Options = autopas::allNewton3Options;

 private:
  double boxLength = -1;
  double cutoff = 1.;
  double cellSizeFactor = 1.;
  double distributionMean = 5.;
  double distributionStdDev = 2.;
  FunctorOption functorOption = FunctorOption::lj12_6;
  GeneratorOption generatorOption = GeneratorOption::grid;
  PrecisionOption precisionOption = PrecisionOption::FP64;
  size_t iterations = 10;
  spdlog::level::level_enum logLevel = spdlog::level::info;
  bool measureFlops = true;
  size_t particlesPerDim = 20;
  size_t particlesTotal = 1000;
  double particleSpacing = .4;
  unsigned int tuningInterval = 100;
  unsigned int tuningSamples = 3;
  string writeVTK = "";
  string logFileName = "";
  unsigned int verletRebuildFrequency = 5;
  unsigned int verletClusterSize = 64;
  double verletSkinRadius = .2;
};
