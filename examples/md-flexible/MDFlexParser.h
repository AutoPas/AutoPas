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
  enum FunctorOption { lj12_6 };
  enum GeneratorOption { grid, gaussian };

  MDFlexParser() = default;

  double getBoxLength() const;
  std::vector<autopas::ContainerOptions> getContainerOptions() const;
  double getCutoff() const;
  autopas::DataLayoutOption getDataLayoutOption() const;
  double getDistributionMean() const;
  double getDistributionStdDev() const;
  FunctorOption getFunctorOption() const;
  GeneratorOption getGeneratorOption() const;
  size_t getIterations() const;
  bool getMeasureFlops() const;
  spdlog::level::level_enum getLogLevel() const;
  double getParticleSpacing() const;
  size_t getParticlesPerDim() const;
  unsigned int getTuningInterval() const;
  string getWriteVTK() const;
  const vector<autopas::TraversalOptions> &getTraversalOptions() const;
  unsigned int getVerletRebuildFrequency() const;
  double getVerletSkinRadius() const;
  bool parseInput(int argc, char **argv);
  void printConfig();

 private:
  static constexpr size_t valueOffset = 32;

  // defaults:
  std::vector<autopas::ContainerOptions> containerOptions = {autopas::ContainerOptions::verletLists};
  autopas::DataLayoutOption dataLayoutOption = autopas::DataLayoutOption::soa;
  std::vector<autopas::TraversalOptions> traversalOptions;

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
  double particleSpacing = .4;
  unsigned int tuningInterval = 100;
  string writeVTK = "";
  unsigned int verletRebuildFrequency = 5;
  double verletSkinRadius = .2;
};
