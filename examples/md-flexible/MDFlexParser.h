#ifndef AUTOPAS_MDFLEXPARSER_H
#define AUTOPAS_MDFLEXPARSER_H

#include <AutoPas.h>
#include <getopt.h>
#include <algorithm>
#include <iomanip>
#include <iostream>

using namespace std;

class MDFlexParser {
 public:
  enum FunctorOption { lj12_6 };

  MDFlexParser() = default;

  bool parseInput(int argc, char **argv);
  autopas::ContainerOption getContainerOption() const;
  double getCutoff() const;
  autopas::DataLayoutOption getDataLayoutOption() const;
  FunctorOption getFunctorOption() const;
  size_t getIterations() const;
  size_t getParticlesPerDim() const;
  double getParticlesSpacing() const;

 private:
  autopas::ContainerOption containerOption;
  autopas::DataLayoutOption dataLayoutOption;
  FunctorOption functorOption;
  size_t particlesPerDim;
  double particleSpacing;
  size_t iterations;
  double cutoff;
};

#endif  // AUTOPAS_MDFLEXPARSER_H
