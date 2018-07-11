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

  MDFlexParser() = default;

  bool parseInput(int argc, char **argv);
  autopas::ContainerOptions getContainerOption() const;
  double getCutoff() const;
  autopas::DataLayoutOption getDataLayoutOption() const;
  FunctorOption getFunctorOption() const;
  size_t getIterations() const;
  size_t getParticlesPerDim() const;
  double getParticleSpacing() const;
  const vector<autopas::TraversalOptions> &getTraversalOptions() const;

 private:
  autopas::ContainerOptions containerOption;
  autopas::DataLayoutOption dataLayoutOption;
  std::vector<autopas::TraversalOptions> traversalOptions;

 private:
  FunctorOption functorOption;
  size_t particlesPerDim;
  double particleSpacing;
  size_t iterations;
  double cutoff;
};
