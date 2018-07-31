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

  autopas::ContainerOptions getContainerOption() const;
  double getCutoff() const;
  autopas::DataLayoutOption getDataLayoutOption() const;
  double getDistributionMean() const;
  double getDistributionStdDev() const;
  FunctorOption getFunctorOption() const;
  GeneratorOption getGeneratorOption() const;
  size_t getIterations() const;
  double getParticleSpacing() const;
  size_t getParticlesPerDim() const;
  const vector<autopas::TraversalOptions> &getTraversalOptions() const;
  size_t getVerletRebuildFrequency() const;
  double getVerletSkinRadius() const;
  bool parseInput(int argc, char **argv);

 private:
  autopas::ContainerOptions containerOption;
  autopas::DataLayoutOption dataLayoutOption;
  std::vector<autopas::TraversalOptions> traversalOptions;

  double cutoff;
  double distributionMean;
  double distributionStdDev;
  FunctorOption functorOption;
  GeneratorOption generatorOption;
  size_t iterations;
  size_t particlesPerDim;
  double particleSpacing;
  size_t verletRebuildFrequency;
  double verletSkinRadius;
};
