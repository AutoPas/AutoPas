/**
 * @file MDFlexParser.h
 * @date 23.02.2018
 * @author F. Gratl
 */

#pragma once
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/SelectorStrategyOption.h"
#include "autopas/options/TraversalOption.h"
#include "autopas/options/TuningStrategyOption.h"
#include "autopas/utils/Logger.h"
#include "autopas/utils/NumberSet.h"

#include <set>

class MDFlexParser {
 public:
  // options specific for the md-flex example

  /**
   * Choice of the functor
   */
  enum class FunctorOption { lj12_6, lj12_6_AVX };

  /**
   * Choice of the particle generator
   */
  enum class GeneratorOption { grid, uniform, gaussian };

  MDFlexParser() = default;

  MDFlexParser(const MDFlexParser &mdparser) = default;

  MDFlexParser &operator=(const MDFlexParser &mdfparser) = default;

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
  autopas::Logger::LogLevel getLogLevel() const;
  double getParticleSpacing() const;
  size_t getParticlesTotal() const;
  size_t getParticlesPerDim() const;
  double getEpsilon() const;
  double getSigma() const;
  double getDeltaT() const;

  unsigned int getTuningInterval() const;
  unsigned int getTuningSamples() const;
  unsigned int getTuningMaxEvidence() const;
  autopas::TuningStrategyOption getTuningStrategyOption() const;
  std::string getVTKFileName() const;
  size_t getVtkWriteFrequency() const;
  const std::set<autopas::TraversalOption> &getTraversalOptions() const;
  unsigned int getVerletRebuildFrequency() const;
  double getVerletSkinRadius() const;
  bool parseInput(int argc, char **argv);
  void printConfig() const;
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
  std::shared_ptr<autopas::NumberSet<double>> cellSizeFactors =
      std::make_shared<autopas::NumberSetFinite<double>>(std::set<double>{1.});

  // Simulation Options:
  double boxLength = -1;
  double cutoff = 1.;
  double delta_t = 0.001;
  double distributionMean = 5.;
  double distributionStdDev = 2.;
  double epsilon = 5.0;
  FunctorOption functorOption = FunctorOption::lj12_6;
  GeneratorOption generatorOption = GeneratorOption::grid;
  size_t iterations = 10;
  autopas::Logger::LogLevel logLevel = autopas::Logger::LogLevel::info;
  std::string logFileName = "";
  double mass = 1.0;
  bool measureFlops = true;
  size_t particlesPerDim = 20;
  size_t particlesTotal = 1000;
  double particleSpacing = .4;
  double sigma = 1.0;
  unsigned int tuningInterval = 100;
  unsigned int tuningSamples = 3;
  unsigned int tuningMaxEvidence = 10;
  unsigned int verletRebuildFrequency = 5;
  std::string vtkFileName = "";
  size_t vtkWriteFrequency = 100;
  double verletSkinRadius = .2;
};
