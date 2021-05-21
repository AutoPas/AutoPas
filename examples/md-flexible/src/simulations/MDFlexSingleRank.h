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
	MDFlexSingleRank(int dimensionCount, int argc, char **argv);
  ~MDFlexSingleRank() = default;

  /**
   * Runs the simulation
   */
  void run() override;

	void initializeDomainDecomposition(int &dimensionCount) override;

	DomainDecomposition* getDomainDecomposition() override {
			return static_cast<DomainDecomposition*>(&(*_domainDecomposition));
	}

	private:
		std::shared_ptr<SingleDomain> _domainDecomposition;

		void printStatistics();
};
