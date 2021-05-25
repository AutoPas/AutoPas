/**
 * @file MDFlexSingleRank.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"
#include "src/domainDecomposition/SingleDomain.h"

/**
 * Runs the MD-Flex simulation on a single node.
 * This is the default demonstration of AutoPas.
 */
class MDFlexSingleRank : public MDFlexSimulation {
 public:
  /**
   * Constructor.
   * @param dimensionCount The number of dimensions in the simulation.
   * @param argc The argument count passed to the main function.
   * @param argv The argument vector passed to the main function.
   */
  MDFlexSingleRank(int dimensionCount, int argc, char **argv);

  /**
   * Destructor.
   */
  ~MDFlexSingleRank() = default;

  /**
   * Runs the simulation
   */
  void run() override;

  /**
   * Initializes the domain decomposition for this simulation.
   */
  void initializeDomainDecomposition(int &dimensionCount) override;

  /**
   * Returns the domain decomposition for this simulation.
   */
  DomainDecomposition *getDomainDecomposition() override {
    return static_cast<DomainDecomposition *>(&(*_domainDecomposition));
  }

 private:
  /**
   * Stores the simulations domain decomposition.
   */
  std::shared_ptr<SingleDomain> _domainDecomposition;

  /**
   * Prints the statistics of the simulation.
   */
  void printStatistics();
};
