/**
 * @file MDFlexSingleNode.h
 * @author J. Körner
 * @date 07.04.2021
 */
#pragma once

#include "MDFlexSimulation.h"

#include "../domainDecomposition/SingleDomain.h"

/**
 * Runs the MD-Flex simulation on a single node.
 * This is the default demonstration of AutoPas.
 */
class MDFlexSingleNode : protected MDFlexSimulation {
 public:
	MDFlexSingleNode(int dimensionCount, int argc, char **argv);
  ~MDFlexSingleNode() = default;

  /**
   * Runs the simulation
   */
  void run() override;
	void initializeDomainDecomposition(int &dimensionCount) override;

	private:
		std::shared_ptr<SingleDomain> _domainDecomposition;

		template <class FunctorType>
		void calculateForces();
};
