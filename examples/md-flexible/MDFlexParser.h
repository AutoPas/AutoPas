#ifndef AUTOPAS_MDFLEXPARSER_H
#define AUTOPAS_MDFLEXPARSER_H

#include <getopt.h>
#include <algorithm>
#include <iomanip>
#include <iostream>

using namespace std;

class MDFlexParser {
 public:
  enum ContainerOption { directSum, linkedCells };
  enum DataLayoutOption { aos, soa };
  enum FunctorOption { lj12_6 };

  MDFlexParser() = default;

  bool parseInput(int argc, char **argv);
  ContainerOption getContainerOption() const;
  DataLayoutOption getDataLayoutOption() const;
  size_t getParticlesPerDim() const;

 private:
  ContainerOption containerOption;
  DataLayoutOption dataLayoutOption;
  FunctorOption functorOption;
  size_t particlesPerDim;
};

#endif  // AUTOPAS_MDFLEXPARSER_H
